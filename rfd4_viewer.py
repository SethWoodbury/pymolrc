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
_detect_spec = importlib.util.spec_from_file_location("_rfd4_cif_detect", _RFD4_DIR / "scripts/rfd4_cif_detect.py")
_RFD4_DETECT = importlib.util.module_from_spec(_detect_spec)
_detect_spec.loader.exec_module(_RFD4_DETECT)
_fast_spec = importlib.util.spec_from_file_location("_rfd4_cif_fast", _RFD4_DIR / "scripts/rfd4_cif_fast.py")
_RFD4_FAST = importlib.util.module_from_spec(_fast_spec)
_fast_spec.loader.exec_module(_RFD4_FAST)
_RFD4_STATE.setdefault("cache_stats", {"memory_hits": 0, "disk_hits": 0, "native_parses": 0,
                                      "last_seconds": 0.0, "last_route": "none"})
_RFD4_STATE["cache_stats"].setdefault("lightweight_reads", 0)


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
    # Supported, fully self-contained CIFs require only the standard CIF
    # tokenizer and NumPy. No native imports or reader-tree filesystem walk.
    # Explicit runtime overrides request that runtime's native parser instead.
    if not repo and not python and not any(os.environ.get(key) for key in ("RFD4_REPO", "RFD4_PYTHON", "RFD4_CLUSTER")):
        before = source.stat()
        try:
            data = _RFD4_FAST.export_model(source)
        except (_RFD4_FAST.UnsupportedCIF, ImportError):
            pass  # Unsupported schemas/chemistry retain the authoritative native path.
        else:
            after = source.stat()
            if (before.st_ino, before.st_size, before.st_mtime_ns, before.st_ctime_ns) != (after.st_ino, after.st_size, after.st_mtime_ns, after.st_ctime_ns):
                raise RuntimeError("Source CIF changed during display parsing; reload its current contents.")
            _RFD4_CACHE.validate_model(data)
            stats = _RFD4_STATE["cache_stats"]
            stats["lightweight_reads"] += 1
            stats.update(last_seconds=round(time.monotonic() - started, 4), last_route="lightweight_reads")
            return data
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
    return _RFD4_DETECT.is_rfd4_cif(filename)


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
    # Display-only ORI markers share the parent object so browsing, transforms,
    # saved sessions and visibility stay together. They are never native atoms.
    physical_keys = {(a["chain"], a["resi"], a["resn"], a["name"]) for a in data["atoms"]}
    centroid_keys = []
    for centroid in data["centroids"]:
        marker = Atom()
        marker.name, marker.symbol, marker.resn, marker.chain = "ORI", "X", "ORI", ""
        marker.resi = str(len(centroid_keys) + 1)
        while (marker.chain, marker.resi, marker.resn, marker.name) in physical_keys | set(centroid_keys):
            marker.resi = str(int(marker.resi) + 1)
        marker.segi, marker.hetatm, marker.q, marker.vdw = "RFD4_CENTROID", 1, 1.0, 1.5
        marker.coord = list(centroid["xyz"])
        marker.id = len(model.atom) + 1
        model.atom.append(marker)
        centroid_keys.append((marker.chain, marker.resi, marker.resn, marker.name))
    _self.load_model(model, obj, zoom=0, discrete=0)
    # PyMOL may sort atoms and replace serial IDs during load. Match identities
    # before building masks so native annotation indices can never shift atoms.
    observed = _self.get_model(obj).atom
    identities = {(atom.chain, atom.resi, atom.resn, atom.name): atom for atom in observed}
    expected_keys = [(atom["chain"], atom["resi"], atom["resn"], atom["name"]) for atom in data["atoms"]]
    if len(identities) != len(observed) or len(set(expected_keys)) != len(expected_keys) or set(identities) != set(expected_keys + centroid_keys):
        _self.delete(obj)
        raise RuntimeError("PyMOL did not preserve unique native atom identities; no annotation masks were applied.")
    for atom, key in zip(data["atoms"], expected_keys):
        atom["pymol_index"] = identities[key].index
    data["centroid_indices"] = [identities[key].index for key in centroid_keys]
    _RFD4_OBJECTS[obj] = data
    if not _automatic:
        _ensure_selections(obj, _self)
    rfd4_style(obj, mode, labels=labels, _self=_self)
    if int(_zoom) == 1 or (int(_zoom) < 0 and first_object):
        _self.zoom(obj)
    if not _automatic:
        rfd4_info(obj, _self=_self)
    else:
        print("[RFD4] %s: %d physical atoms, %d bonds, %d centroid markers; %d scaffold/unoccupied atoms hidden." %
              (obj, len(data["atoms"]), len(data["bonds"]), len(data["centroids"]), data["excluded_atoms"]))
    return obj


def _ensure_selections(obj, _self=cmd, coordinates_only=False):
    data = _RFD4_OBJECTS[obj]
    if any("internal" not in atom for atom in data["atoms"]):
        raise ValueError("This object predates coordinate-group metadata; reload its source CIF into a new object.")
    ready = data.get("selections_created") or (coordinates_only and data.get("coordinate_selections_created"))
    if ready and _selection(obj, "fixed") in _self.get_names("all"):
        return
    selections = {
        "fixed": lambda atom: atom["fixed"], "absolute": lambda atom: atom["fixed"],
        "internal": lambda atom: atom["internal"],
        "reference": lambda atom: not atom["fixed"] and not atom["internal"],
        "sequence_fixed": lambda atom: atom["sequence_fixed"],
        "sequence_free": lambda atom: not atom["sequence_fixed"], "unindexed": lambda atom: atom["unindexed"],
        "hb_donor": lambda atom: atom["hbond"] in {"DONOR", "BOTH"},
        "hb_acceptor": lambda atom: atom["hbond"] in {"ACCEPTOR", "BOTH"},
        "hb_both": lambda atom: atom["hbond"] == "BOTH",
        "buried": lambda atom: atom["rasa"] == "BURIED", "partial": lambda atom: atom["rasa"] == "INTERMEDIATE",
        "exposed": lambda atom: atom["rasa"] == "EXPOSED",
    }
    if coordinates_only:
        selections = {key: predicate for key, predicate in selections.items()
                      if key in {"fixed", "absolute", "internal", "reference"}}
    for suffix, predicate in selections.items():
        _select_atoms(obj, suffix, predicate, data, _self)
    _self.group(obj + "_annotations", " ".join(_selection(obj, suffix) for suffix in selections), action="add")
    data["coordinate_selections_created"] = True
    if not coordinates_only:
        data["selections_created"] = True


def _element_palette(_self):
    """Use the user's named palette when installed, with standalone fallbacks."""
    choices = {
        "protein": ("sc05", "green"), "backbone": ("bb05", "forest"),
        "ligand": ("lig_amethyst3", "violet"), "metal": ("elemZn", "gray60"),
        "N": ("elemN", "blue"), "O": ("elemO", "red"),
        "S": ("elemS", "yellow"), "P": ("elemP", "orange"), "H": ("elemH", "white"),
    }
    return {key: named if _self.get_color_index(named) >= 0 else fallback
            for key, (named, fallback) in choices.items()}


def _shade_free_identity(obj, _self):
    """Lighten protein carbons only when the whole residue has free identity."""
    fixed_residues = {(a["chain"], a["resi"], a["resn"]) for a in _RFD4_OBJECTS[obj]["atoms"]
                      if a["sequence_fixed"]}
    indices = [str(a["pymol_index"]) for a in _RFD4_OBJECTS[obj]["atoms"]
               if a["element"].upper() == "C" and (a["chain"], a["resi"], a["resn"]) not in fixed_residues]
    if not indices:
        return
    colors = {}
    selection = "model %s and polymer.protein and index %s" % (obj, "+".join(indices))
    _self.iterate(selection, "colors.setdefault(color, []).append(str(index))", space={"colors": colors})
    for color, atom_indices in colors.items():
        rgb = tuple(value + 0.30 * (1.0 - value) for value in _self.get_color_tuple(color))
        name = "rfd4_free_" + "".join("%02x" % round(value * 255) for value in rgb)
        _self.set_color(name, rgb, quiet=1)
        _self.color(name, "model %s and index %s" % (obj, "+".join(atom_indices)))


def _annotated_object(object, _self):
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
    return obj


def _style_centroids(obj, _self):
    data = _RFD4_OBJECTS[obj]
    if "centroids" not in data:
        print("[RFD4] Reload the source CIF to display its centroid targets.")
        return
    indices = data.get("centroid_indices", [])
    if not indices:
        print("[RFD4] No centroids.")
        return
    selection = _selection(obj, "centroids")
    _self.select(selection, "model %s and index %s" % (obj, "+".join(map(str, indices))), quiet=1)
    _self.disable(selection)
    _self.group(obj + "_annotations", selection, action="add")
    _self.hide("everything", selection)
    _self.label(selection, "")
    if not data.get("show_centroids", True):
        print("[RFD4] Centroid markers hidden.")
        return
    _self.color("hotpink", selection)
    _self.set("sphere_scale", 1.0, selection)
    _self.set("sphere_transparency", 0.55, selection)
    _self.set("label_color", "hotpink", selection)
    _self.show("spheres", selection)
    _self.show("nonbonded", selection)
    for index, centroid in zip(indices, data["centroids"]):
        _self.label("model %s and index %d" % (obj, index), repr("Centroid " + centroid["chain"]))
    print("[RFD4] Centroids: " + ", ".join(c["chain"] for c in data["centroids"]) +
          " (magenta ORI markers; conditioning targets, not atoms).")


def rfd4_centroids(object="", show=1, _self=cmd):
    """Show/hide chain-centroid markers: rfd4_centroids [object], 1|0."""
    obj = _annotated_object(object, _self)
    _RFD4_OBJECTS[obj]["show_centroids"] = bool(int(show))
    _style_centroids(obj, _self)


def rfd4_style(object="", mode="coordinates", labels=1, context=1, element_colors=1, _self=cmd):
    """Switch annotation colors: coordinates | sequence | unindexed | hbond | burial | elements."""
    obj = _annotated_object(object, _self)
    _ensure_selections(obj, _self, coordinates_only=mode in {"coordinates", "elements"})
    palettes = {
        "coordinates": [("absolute", "cyan"), ("internal", "orange")],
        "sequence": [("sequence_fixed", "teal"), ("sequence_free", "gray70")],
        "unindexed": [("unindexed", "violet")],
        "hbond": [("hb_donor", "marine"), ("hb_acceptor", "salmon"), ("hb_both", "magenta")],
        "burial": [("buried", "forest"), ("partial", "orange"), ("exposed", "skyblue")],
        "elements": [],
    }
    if mode not in palettes:
        raise ValueError("Mode must be " + ", ".join(palettes))
    coordinate_view = mode in {"coordinates", "elements"}
    chemical_colors = mode == "elements" or (mode == "coordinates" and int(element_colors))
    show_context = coordinate_view and int(context)
    _self.hide("everything", obj)
    shown = (_selection(obj, "absolute") + " or " + _selection(obj, "internal")) if coordinate_view else obj
    _self.show("sticks", shown)
    _self.show("nonbonded", shown)
    if show_context:
        _self.show("lines", _selection(obj, "reference"))
        _self.show("nonbonded", _selection(obj, "reference"))
    _self.set("stick_radius", 0.13, obj)
    # Clear older viewer overrides so metal spheres inherit live global changes.
    metal_elements = "Zn+Fe+Mg+Mn+Cu+Co+Ni+Ca+Na+K+Mo+V+W+Cd+Hg+Pt+Pd+Au+Ag+Al+Ba+Sr+Li+Cs+Rb+Cr"
    _self.unset("sphere_scale", obj)
    _self.unset("sphere_scale", "model " + obj)
    _self.set("sphere_scale", 0.25, "model %s and not (elem %s)" % (obj, metal_elements))
    _self.color("gray60", obj)
    if chemical_colors:
        colors = _element_palette(_self)
        _self.color("atomic", obj)
        _self.color("cyan" if mode == "coordinates" else colors["protein"],
                    "model %s and polymer.protein and elem C" % obj)
        if mode == "elements":
            _self.color(colors["backbone"], "model %s and polymer.protein and elem C and backbone" % obj)
        _self.color(colors["ligand"], "model %s and not polymer.protein and elem C" % obj)
        for element in ("N", "O", "S", "P", "H"):
            _self.color(colors[element], "model %s and elem %s" % (obj, element))
        if mode == "coordinates":
            # Internal geometry takes precedence for carbon atoms in any component.
            _self.color("orange", _selection(obj, "internal") + " and elem C")
        _shade_free_identity(obj, _self)
    else:
        for suffix, color in palettes[mode]:
            _self.color(color, _selection(obj, suffix))
    if mode == "coordinates":
        _self.show("spheres", shown)
    elif mode == "hbond":
        _self.show("spheres", _selection(obj, "hb_donor") + " or " + _selection(obj, "hb_acceptor"))
    # Atom elements avoid the protein CA / calcium and CD / cadmium collisions.
    # Metal atoms keep the global sphere size; only nonmetals get small spheres.
    visible = "(" + shown + ")"
    if show_context:
        visible += " or " + _selection(obj, "reference")
    metals = "(%s) and elem %s" % (shown, metal_elements)
    if _self.count_atoms(metals):
        _self.show("spheres", metals)
        if chemical_colors:
            _self.color(colors["metal"], metals)
    if show_context:
        # The startup show hook adds metal spheres to lines/nonbonded. Undo that
        # for references: unconstrained atoms stay lines or unbonded crosses.
        _self.hide("spheres", _selection(obj, "reference"))
    _self.label(obj, "")
    if int(labels):
        _self.label("(%s) and name CA" % visible, '"%s/%s%s" % (chain,resn,resi)')
    _style_centroids(obj, _self)
    if chemical_colors:
        legend = "element colors; ligand/cofactor C=%s; metals=%s; free-identity protein C=lighter shade" % (colors["ligand"], colors["metal"])
        legend += "; protein absolute C=cyan; internal C=orange" if mode == "coordinates" else "; protein C=" + colors["protein"]
    else:
        legend = ", ".join(suffix + "=" + color for suffix, color in palettes[mode])
    print("[RFD4] " + mode + ": " + legend)
    if mode == "coordinates":
        atoms = _RFD4_OBJECTS[obj]["atoms"]
        absolute, internal = sum(a["fixed"] for a in atoms), sum(a["internal"] for a in atoms)
        print("[RFD4] %d absolute targets; %d internal-group atoms (group placement free); "
              "%d reference-only atoms %s. Sequence-fixed does not mean coordinate-fixed." %
              (absolute, internal, len(atoms) - absolute - internal, "in lines (unbonded atoms as crosses)" if show_context else "hidden"))


def rfd4_info(object, _self=cmd):
    """Print condition counts, retained source path, and expandable scaffold instructions."""
    data = _RFD4_OBJECTS.get(str(object))
    if data is None:
        raise ValueError("Load this object with rfd4_load first.")
    if any("internal" not in atom for atom in data["atoms"]):
        raise ValueError("This object predates coordinate-group metadata; reload its source CIF into a new object.")
    atoms = data["atoms"]
    print("[RFD4] Source: " + data["source"])
    print("[RFD4] %d physical atoms; %d native bonds; %d placeholders/unoccupied atoms omitted." %
          (len(atoms), len(data["bonds"]), data["excluded_atoms"]))
    for field in ("fixed", "sequence_fixed", "unindexed"):
        print("  %-15s %d atoms%s" % ("absolute" if field == "fixed" else field, sum(atom[field] for atom in atoms),
              "" if data["present"][field] else " (annotation absent; native default)"))
    print("  internal        %d atoms (relative geometry; group placement is not fixed)." % sum(a["internal"] for a in atoms))
    print("  reference-only  %d retained atoms; lines in the default coordinate view (context=0 hides them)." %
          sum(not a["fixed"] and not a["internal"] for a in atoms))
    if "centroids" not in data:
        print("  centroid metadata unavailable; reload the source CIF.")
    elif not data["centroids"]:
        print("  No centroids.")
    else:
        for centroid in data["centroids"]:
            print("  centroid %-8s xyz = %s Å (stored target)" % (centroid["chain"], centroid["xyz"]))
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
  pymol /path/input.cif             from a terminal: show coordinate-conditioned atoms
  load /path/input.cif              inside PyMOL: automatic native parsing
  rfd4_style , hbond                color the single visible annotated object
  rfd4_load /absolute/path/input.cif, motif
  rfd4_style motif, coordinates    element colors; cyan protein / purple ligand C; orange internal C
  rfd4_style motif, coordinates, element_colors=0   flat cyan absolute / orange internal atoms
  rfd4_style motif, coordinates, context=0   hide reference-only context lines
  rfd4_style motif, sequence       teal fixed identities; gray redesignable identities
  rfd4_style motif, unindexed      violet freely placed sequence groups
  rfd4_style motif, hbond          blue donor / salmon acceptor / magenta both
  rfd4_style motif, burial         green buried / orange partial / blue exposed
  rfd4_style motif, elements       ordinary protein/ligand colors; free-identity carbons lighter
  rfd4_centroids motif, 0          hide magenta centroid markers (1 = show)
  rfd4_info motif                  source, condition counts, scaffold-length ranges
  rfd4_cache_info                   cache hits and last model-read timing

Each object has named selections under <object>_annotations. Coordinates and
native bonds are preserved. Only asymmetric-unit/model-1/first-altloc physical
atoms are shown; placeholders and zero-occupancy atoms are omitted. Labels show
source residue IDs, not final designed positions. Additional condition types are
retained in the source CIF; this viewer colors the five listed categories.
Coordinate targets are sticks/small spheres; unconstrained reference atoms are
lines (unbonded atoms are crosses). Free-identity protein carbons use a lighter
shade of their primary color. Heteroatom colors do not encode sequence identity.
Centroid targets appear as labeled translucent magenta ORI spheres/crosses,
one per conditioned chain. They are display-only markers in the source frame,
not physical atoms or added inference tokens. No targets prints "No centroids".
Use the original annotated CIF for inference; PyMOL exports include display markers.

Supported self-contained CIFs use a lightweight CIF reader directly in PyMOL.
Other schemas use verified JSON display caches; a native cache miss runs
~/git/RFD4-Proteina-dev/.venv/bin/python on CPU once.
Overrides: RFD4_REPO/RFD4_PYTHON/RFD4_CLUSTER environment variables, or repo=/python=
on rfd4_load. No remote PyMOL server or model/checkpoint loading is required.
""")


cmd.extend("rfd4_load", rfd4_load)
cmd.extend("rfd4_style", rfd4_style)
cmd.extend("rfd4_centroids", rfd4_centroids)
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
    rfd4_load(filename, object=bound.get("object", ""), mode="coordinates", labels=0,
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
