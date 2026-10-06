"""Automatic, read-only PyMOL loading and styling of RFD4 annotated CIFs.

Only the subprocess uses RFD4/AtomWorks; normal PyMOL stays in its own Python.
Load with ``run ~/pymolrc/rfd4_viewer.py``; enter ``help_rfd4`` for usage.
"""
import atexit
from collections import OrderedDict
from copy import deepcopy
import hashlib
import importlib.util
import json
import inspect
import os
from pathlib import Path
import re
import subprocess
import select
import tempfile
import threading
import time

from chempy import Atom, Bond
from chempy.models import Indexed
from pymol import cmd

# Store session state on cmd so re-running this script preserves loaded metadata
# and never nests loader wrappers or spawns another parser unnecessarily.
if not hasattr(cmd, "_rfd4_viewer_state"):
    cmd._rfd4_viewer_state = {"objects": {}, "cache": OrderedDict(), "worker": None,
                              "worker_key": None, "lock": threading.Lock()}
_RFD4_STATE = cmd._rfd4_viewer_state
_RFD4_OBJECTS = _RFD4_STATE["objects"]
_RFD4_DIR = Path(inspect.currentframe().f_code.co_filename).resolve().parent
_cache_spec = importlib.util.spec_from_file_location("_rfd4_model_cache", _RFD4_DIR / "scripts/rfd4_model_cache.py")
_RFD4_CACHE = importlib.util.module_from_spec(_cache_spec)
_cache_spec.loader.exec_module(_RFD4_CACHE)
_RFD4_STATE.setdefault("cache_stats", {"memory_hits": 0, "disk_hits": 0, "native_parses": 0,
                                      "last_seconds": 0.0, "last_route": "none"})


def _selection(obj, suffix):
    return obj + "_" + suffix


def _select_atoms(obj, suffix, predicate, data, _self):
    identifiers = [str(atom["pymol_index"]) for atom in data["atoms"] if predicate(atom)]
    expression = "model " + obj + " and index " + "+".join(identifiers) if identifiers else "none"
    _self.select(_selection(obj, suffix), expression, quiet=1)
    _self.disable(_selection(obj, suffix))


def _stop_native_worker():
    worker = _RFD4_STATE.get("worker")
    if worker is not None:
        if worker.poll() is None:
            worker.terminate()
            try:
                worker.wait(timeout=3)
            except subprocess.TimeoutExpired:
                worker.kill()
                worker.wait(timeout=3)
        worker.stdin.close()
        worker.stdout.close()
        worker._rfd4_stderr.close()
    _RFD4_STATE["worker"] = None


def _native_model(source, repo="", python=""):
    started = time.monotonic()
    config = _RFD4_CACHE.reader_config(repo, python)
    reader = _RFD4_CACHE.provenance(config)
    repo_path, interpreter, exporter = config["repo"], config["python"], config["exporter"]
    worker_key = reader["sha256"]
    source_hash = hashlib.sha256(source.read_bytes()).hexdigest()
    key = (source_hash, worker_key)
    stats = _RFD4_STATE["cache_stats"]

    def record(data, route):
        stats[route] += 1
        stats.update(last_seconds=round(time.monotonic() - started, 4), last_route=route)
        data = deepcopy(data)
        data["source"] = str(source)
        return data

    with _RFD4_STATE["lock"]:
        if key in _RFD4_STATE["cache"]:
            _RFD4_STATE["cache"].move_to_end(key)
            return record(_RFD4_STATE["cache"][key], "memory_hits")
        data = _RFD4_CACHE.read_cached(source, source_hash, worker_key)
        if data is not None:
            _RFD4_STATE["cache"][key] = data
            while len(_RFD4_STATE["cache"]) > 32:
                _RFD4_STATE["cache"].popitem(last=False)
            return record(data, "disk_hits")
        worker = _RFD4_STATE.get("worker")
        if worker is None or worker.poll() is not None or _RFD4_STATE["worker_key"] != worker_key:
            _stop_native_worker()
            environment = dict(os.environ)
            for name in ("PYTHONHOME", "PYTHONEXECUTABLE", "__PYVENV_LAUNCHER__"):
                environment.pop(name, None)
            environment.update(CLUSTER=config["cluster"], CUDA_VISIBLE_DEVICES="", OMP_NUM_THREADS="1",
                               OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
            error_log = tempfile.TemporaryFile(mode="w+t")
            worker = subprocess.Popen([str(interpreter), str(exporter), "--serve"], cwd=str(repo_path),
                                      env=environment, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                      stderr=error_log, text=True, bufsize=1)
            worker._rfd4_stderr = error_log
            _RFD4_STATE.update(worker=worker, worker_key=worker_key)
            print("[RFD4] Preparing a new display cache on CPU (reused across PyMOL sessions).")
        try:
            worker.stdin.write(json.dumps({"source": str(source)}) + "\n")
            worker.stdin.flush()
            if not select.select([worker.stdout], [], [], 180)[0]:
                raise RuntimeError("Native RFD4 reader timed out after 180 seconds.")
            line = worker.stdout.readline()
            if not line:
                worker._rfd4_stderr.seek(0)
                raise RuntimeError("Native reader exited:\n" + worker._rfd4_stderr.read()[-6000:])
            response = json.loads(line)
        except Exception:
            _stop_native_worker()
            raise
        if not response["ok"]:
            raise RuntimeError("Native RFD4 parsing failed:\n" + response["error"][-6000:])
        data = response["model"]
        _RFD4_CACHE.validate_model(data)
        if hashlib.sha256(source.read_bytes()).hexdigest() != source_hash:
            raise RuntimeError("Source CIF changed during native parsing; reload its current contents.")
        if _RFD4_CACHE.provenance(config)["sha256"] != worker_key:
            _stop_native_worker()
            raise RuntimeError("RFD4 reader installation changed during parsing; reload to use the new reader.")
        try:
            _RFD4_CACHE.write_cached(source_hash, reader, data)
        except OSError as error:
            print("[RFD4] Could not save the disposable display cache: " + str(error))
        _RFD4_STATE["cache"][key] = data
        while len(_RFD4_STATE["cache"]) > 32:
            _RFD4_STATE["cache"].popitem(last=False)
        return record(data, "native_parses")


def rfd4_cache_info():
    """Report disk/memory cache hits and reader timing for this PyMOL session."""
    print("[RFD4] Display cache: " + str(_RFD4_CACHE.cache_root()))
    print("[RFD4] " + json.dumps(_RFD4_STATE["cache_stats"], sort_keys=True))


def _is_rfd4_cif(filename):
    try:
        source = Path(str(filename).strip('"\'')).expanduser()
        if source.suffix.lower() not in {".cif", ".mmcif"} or not source.is_file():
            return False
        # A normal CIF follows the original loader unchanged. Test actual CIF
        # field declarations, not filenames, so descriptive renames are harmless.
        with source.open(errors="replace") as handle:
            return any(line.lstrip().lower().startswith(("_condition_", "_atom_site.condition_",
                                        "_atom_site.annotation_expseg")) for line in handle)
    except (OSError, ValueError):
        return False


def rfd4_load(filename, object="", mode="coordinates", labels=1, repo="", python="", _self=cmd,
              _automatic=False, _zoom=-1):
    """Load an annotated CIF; normal ``load`` and command-line opens also work.

    rfd4_load /path/input.cif [, object=motif] [, mode=coordinates] [, labels=1]
    """
    source = Path(str(filename).strip('"\'')).expanduser().resolve(strict=True)
    if mode not in {"coordinates", "sequence", "unindexed", "hbond", "burial", "elements"}:
        raise ValueError("Mode must be coordinates, sequence, unindexed, hbond, burial, or elements.")
    stem = re.sub(r"[^A-Za-z0-9_]", "_", source.stem)
    proposed = str(object) or (stem if _automatic else "rfd4_" + stem)
    if not re.fullmatch(r"[A-Za-z0-9_]+", proposed):
        raise ValueError("Choose an object name using letters, numbers and underscores.")
    # Never append an incompatible motif to an existing object, or overwrite it.
    obj = _self.get_unused_name(proposed) if proposed in _self.get_names("all") else proposed
    if any(name.startswith(obj + "_") for name in _self.get_names("all")):
        obj = _self.get_unused_name(obj + "_view")
    data = _native_model(source, repo=repo, python=python)
    live_objects = set(_self.get_object_list())
    for stale in set(_RFD4_OBJECTS) - live_objects:
        del _RFD4_OBJECTS[stale]
    first_object = not live_objects
    model = Indexed()
    for i, record in enumerate(data["atoms"], 1):
        atom = Atom()
        atom.id, atom.index = i, i
        atom.name, atom.symbol, atom.chain = record["name"], record["element"], record["chain"]
        atom.resi, atom.resn, atom.coord = record["resi"], record["resn"], record["coord"]
        atom.hetatm, atom.q, atom.formal_charge = int(record["hetero"]), record["occupancy"], record["charge"]
        model.atom.append(atom)
    for i, j, order in data["bonds"]:
        bond = Bond()
        bond.index, bond.order = [i, j], order
        model.bond.append(bond)
    _self.load_model(model, obj, zoom=0, discrete=0)
    # PyMOL may sort atoms and replace serial IDs during load. Match identities
    # before building masks so native annotation indices can never shift atoms.
    observed = _self.get_model(obj).atom
    identities = {(atom.chain, atom.resi, atom.resn, atom.name): atom for atom in observed}
    expected_keys = [(atom["chain"], atom["resi"], atom["resn"], atom["name"]) for atom in data["atoms"]]
    if len(identities) != len(observed) or len(set(expected_keys)) != len(expected_keys) or set(identities) != set(expected_keys):
        _self.delete(obj)
        raise RuntimeError("PyMOL did not preserve unique native atom identities; no annotation masks were applied.")
    for atom, key in zip(data["atoms"], expected_keys):
        atom["pymol_index"] = identities[key].index
    _RFD4_OBJECTS[obj] = data
    if not _automatic:
        _ensure_selections(obj, _self)
    rfd4_style(obj, mode, labels=labels, _self=_self)
    if int(_zoom) == 1 or (int(_zoom) < 0 and first_object):
        _self.zoom(obj)
    if not _automatic:
        rfd4_info(obj, _self=_self)
    else:
        print("[RFD4] %s: %d atoms, %d bonds; %d scaffold/unoccupied atoms hidden." %
              (obj, len(data["atoms"]), len(data["bonds"]), data["excluded_atoms"]))
    return obj


def _ensure_selections(obj, _self=cmd):
    data = _RFD4_OBJECTS[obj]
    if data.get("selections_created") and _selection(obj, "fixed") in _self.get_names("all"):
        return
    selections = {
        "fixed": lambda atom: atom["fixed"], "sequence_fixed": lambda atom: atom["sequence_fixed"],
        "sequence_free": lambda atom: not atom["sequence_fixed"], "unindexed": lambda atom: atom["unindexed"],
        "hb_donor": lambda atom: atom["hbond"] in {"DONOR", "BOTH"},
        "hb_acceptor": lambda atom: atom["hbond"] in {"ACCEPTOR", "BOTH"},
        "hb_both": lambda atom: atom["hbond"] == "BOTH",
        "buried": lambda atom: atom["rasa"] == "BURIED", "partial": lambda atom: atom["rasa"] == "INTERMEDIATE",
        "exposed": lambda atom: atom["rasa"] == "EXPOSED",
    }
    for suffix, predicate in selections.items():
        _select_atoms(obj, suffix, predicate, data, _self)
    _self.group(obj + "_annotations", " ".join(_selection(obj, suffix) for suffix in selections), action="add")
    data["selections_created"] = True


def rfd4_style(object="", mode="coordinates", labels=1, _self=cmd):
    """Switch annotation colors: coordinates | sequence | unindexed | hbond | burial | elements."""
    obj = str(object)
    if not obj:
        candidates = [name for name in _self.get_object_list() if name in _RFD4_OBJECTS]
        visible = set(_self.get_names("objects", enabled_only=1))
        active = [name for name in candidates if name in visible]
        if len(active) == 1:
            obj = active[0]
        elif len(candidates) == 1:
            obj = candidates[0]
    if obj not in _RFD4_OBJECTS or obj not in _self.get_names("objects"):
        available = [name for name in _self.get_object_list() if name in _RFD4_OBJECTS]
        raise ValueError("Use the actual object name, or omit it when only one annotated object is visible. "
                         "Annotated objects: " + (", ".join(available) or "none; reopen the source CIF"))
    if mode != "elements":
        _ensure_selections(obj, _self)
    palettes = {
        "coordinates": [("fixed", "cyan")],
        "sequence": [("sequence_fixed", "teal"), ("sequence_free", "gray70")],
        "unindexed": [("unindexed", "violet")],
        "hbond": [("hb_donor", "marine"), ("hb_acceptor", "salmon"), ("hb_both", "magenta")],
        "burial": [("buried", "forest"), ("partial", "orange"), ("exposed", "skyblue")],
        "elements": [],
    }
    if mode not in palettes:
        raise ValueError("Mode must be " + ", ".join(palettes))
    _self.hide("everything", obj)
    _self.show("sticks", obj)
    _self.show("nonbonded", obj)
    _self.set("stick_radius", 0.13, obj)
    _self.set("sphere_scale", 0.25, obj)
    _self.color("gray60", obj)
    if mode == "elements":
        for element, color in [("N", "blue"), ("O", "red"), ("S", "yellow"), ("P", "orange"), ("H", "white")]:
            _self.color(color, "model %s and elem %s" % (obj, element))
    for suffix, color in palettes[mode]:
        _self.color(color, _selection(obj, suffix))
    if mode == "coordinates":
        _self.show("spheres", _selection(obj, "fixed"))
    elif mode == "hbond":
        _self.show("spheres", _selection(obj, "hb_donor") + " or " + _selection(obj, "hb_acceptor"))
    _self.label(obj, "")
    if int(labels):
        _self.label("model %s and name CA" % obj, '"%s/%s%s" % (chain,resn,resi)')
    print("[RFD4] " + mode + ": " + (", ".join(suffix + "=" + color for suffix, color in palettes[mode]) or "element colors"))
    if mode == "coordinates":
        print("[RFD4] Gray atoms are reference geometry without fixed-coordinate targets.")


def rfd4_info(object, _self=cmd):
    """Print condition counts, retained source path, and expandable scaffold instructions."""
    data = _RFD4_OBJECTS.get(str(object))
    if data is None:
        raise ValueError("Load this object with rfd4_load first.")
    atoms = data["atoms"]
    print("[RFD4] Source: " + data["source"])
    print("[RFD4] %d physical atoms; %d native bonds; %d placeholders/unoccupied atoms omitted." %
          (len(atoms), len(data["bonds"]), data["excluded_atoms"]))
    for field in ("fixed", "sequence_fixed", "unindexed"):
        print("  %-15s %d atoms%s" % (field, sum(atom[field] for atom in atoms),
              "" if data["present"][field] else " (annotation absent; native default)"))
    print("  H-bond roles: %d donor atoms; %d acceptor atoms." %
          (sum(atom["hbond"] in {"DONOR", "BOTH"} for atom in atoms),
           sum(atom["hbond"] in {"ACCEPTOR", "BOTH"} for atom in atoms)))
    print("  Burial: %d buried; %d partially buried; %d exposed atoms." %
          tuple(sum(atom["rasa"] == state for atom in atoms) for state in ("BURIED", "INTERMEDIATE", "EXPOSED")))
    for segment in data["segments"]:
        print("  Expandable segment %s/%s: %s-%s residues (not displayed)." %
              (segment["chain"], segment["resi"], segment["min"], segment["max"]))
    print("  H-bond roles and burial labels are requested conditions, not measured contacts/accessibility.")


def help_rfd4():
    print("""
RFD4 annotated-CIF inspection (read-only; no generation)
  pymol /path/input.cif             from a terminal: automatic ordinary-looking view
  load /path/input.cif              inside PyMOL: automatic native parsing
  rfd4_style , hbond                color the single visible annotated object
  rfd4_load /absolute/path/input.cif, motif
  rfd4_style motif, coordinates    cyan fixed atoms; gray reference-only atoms
  rfd4_style motif, sequence       teal fixed identities; gray redesignable identities
  rfd4_style motif, unindexed      violet freely placed sequence groups
  rfd4_style motif, hbond          blue donor / salmon acceptor / magenta both
  rfd4_style motif, burial         green buried / orange partial / blue exposed
  rfd4_style motif, elements       conventional element colors
  rfd4_info motif                  source, condition counts, scaffold-length ranges
  rfd4_cache_info                   cache hits and last model-read timing

Each object has named selections under <object>_annotations. Coordinates and
native bonds are preserved. Only asymmetric-unit/model-1/first-altloc physical
atoms are shown; placeholders and zero-occupancy atoms are omitted. Labels show
source residue IDs, not final designed positions. Additional condition types are
retained in the source CIF; this viewer colors the five listed categories.

Valid JSON display caches load directly in PyMOL without a native subprocess.
On a cache miss, ~/git/RFD4-Proteina-dev/.venv/bin/python parses on CPU once.
Overrides: RFD4_REPO/RFD4_PYTHON/RFD4_CLUSTER environment variables, or repo=/python=
on rfd4_load. No remote PyMOL server or model/checkpoint loading is required.
""")


cmd.extend("rfd4_load", rfd4_load)
cmd.extend("rfd4_style", rfd4_style)
cmd.extend("rfd4_info", rfd4_info)
cmd.extend("help_rfd4", help_rfd4)
cmd.extend("rfd4_cache_info", rfd4_cache_info)


def _automatic_load(*args, **kwargs):
    original = _RFD4_STATE["original_load"]
    filename = args[0] if args else kwargs.get("filename", "")
    if not _is_rfd4_cif(filename):
        return original(*args, **kwargs)
    bound = inspect.signature(original).bind_partial(*args, **kwargs).arguments
    if int(bound.get("state", 0)) not in (0, 1):
        raise ValueError("Annotated inputs are single-state instructions; load each as a separate object.")
    rfd4_load(filename, object=bound.get("object", ""), mode="elements", labels=0,
              _self=bound.get("_self", cmd), _automatic=True, _zoom=bound.get("zoom", -1))
    return 1


def _install_autoload():
    current = cmd.load
    # Splice beneath the user's existing hook, preserving autosolo, ordinary
    # PDB/CIF loading and per-object styling. Also wire PyMOL's textual command
    # table: command-line file opens use that table, not cmd.load directly.
    namespace = getattr(current, "__globals__", {})
    if current.__name__ == "_pymolrc_load_hook" and "_pymolrc_orig_load" in namespace:
        underlying = namespace["_pymolrc_orig_load"]
        if not getattr(underlying, "_rfd4_dispatch", False):
            _RFD4_STATE["original_load"] = underlying
        namespace["_pymolrc_orig_load"] = _automatic_load
    else:
        if not getattr(current, "_rfd4_dispatch", False):
            _RFD4_STATE["original_load"] = current
        cmd.load = _automatic_load
    cmd.keyword["load"][0] = cmd.load


_automatic_load._rfd4_dispatch = True
_install_autoload()
if not _RFD4_STATE.get("atexit_registered"):
    atexit.register(_stop_native_worker)
    _RFD4_STATE["atexit_registered"] = True
