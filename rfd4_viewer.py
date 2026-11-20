"""Opt-in, read-only PyMOL loading and styling of RFD4 annotated CIFs.

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
def _viewer_module(name):
    spec = importlib.util.spec_from_file_location(name, _RFD4_DIR / "scripts" / (name + ".py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


_RFD4_CACHE = _viewer_module("rfd4_model_cache")
_RFD4_DETECT = _viewer_module("rfd4_cif_detect")
_RFD4_FAST = _viewer_module("rfd4_cif_fast")
if "paint_component_carbons" not in globals():
    paint_component_carbons = _viewer_module("pymol_components").paint_component_carbons
_RFD4_STATE.setdefault("cache_stats", {"memory_hits": 0, "disk_hits": 0, "native_parses": 0,
                                      "last_seconds": 0.0, "last_route": "none"})
_RFD4_STATE["cache_stats"].setdefault("lightweight_reads", 0)


def _selection(obj, suffix):
    return obj + "_" + suffix


def _select_atoms(obj, suffix, predicate, data, _self):
    indices = [atom["pymol_index"] for atom in data["atoms"] if predicate(atom)]
    name = _selection(obj, suffix)
    _self.select_list(name, obj, indices, mode="index")
    _self.disable(name)


def _stop_native_worker():
    reader = _RFD4_STATE.get("reader")
    if reader is None and _RFD4_STATE.get("worker") is not None:
        # A viewer reloaded in a live session can inherit an older worker.
        reader = _RFD4_CACHE.NativeReader({})
        reader.worker = _RFD4_STATE["worker"]
    if reader is not None:
        reader.close()
    _RFD4_STATE.update(reader=None, worker=None)


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
    worker_key = reader["sha256"]
    source_hash = hashlib.sha256(source.read_bytes()).hexdigest()
    key = (source_hash, worker_key)
    stats = _RFD4_STATE["cache_stats"]

    def record(data, route):
        cache = _RFD4_STATE["cache"]
        cache[key] = data
        cache.move_to_end(key)
        if len(cache) > 32:
            cache.popitem(last=False)
        stats[route] += 1
        stats.update(last_seconds=round(time.monotonic() - started, 4), last_route=route)
        data = deepcopy(data)
        data["source"] = str(source)
        return data

    with _RFD4_STATE["lock"]:
        if key in _RFD4_STATE["cache"]:
            return record(_RFD4_STATE["cache"][key], "memory_hits")
        data = _RFD4_CACHE.read_cached(source, source_hash, worker_key)
        if data is not None:
            return record(data, "disk_hits")
        if _RFD4_STATE.get("reader") is None or _RFD4_STATE["worker_key"] != worker_key:
            _stop_native_worker()
            _RFD4_STATE.update(reader=_RFD4_CACHE.NativeReader(config), worker_key=worker_key)
        native = _RFD4_STATE["reader"]
        print("[RFD4] Reading uncached annotations with the native CPU parser.")
        try:
            data = native.model(source)
        finally:
            _RFD4_STATE["worker"] = native.worker
        if hashlib.sha256(source.read_bytes()).hexdigest() != source_hash:
            raise RuntimeError("Source CIF changed during native parsing; reload its current contents.")
        if _RFD4_CACHE.provenance(config)["sha256"] != worker_key:
            _stop_native_worker()
            raise RuntimeError("RFD4 reader installation changed during parsing; reload to use the new reader.")
        try:
            _RFD4_CACHE.write_cached(source_hash, reader, data)
        except OSError as error:
            print("[RFD4] Could not save the disposable display cache: " + str(error))
        return record(data, "native_parses")


def rfd4_cache_info():
    """Report disk/memory cache hits and reader timing for this PyMOL session."""
    print("[RFD4] Display cache: " + str(_RFD4_CACHE.cache_root()))
    print("[RFD4] " + json.dumps(_RFD4_STATE["cache_stats"], sort_keys=True))


def _is_rfd4_cif(filename):
    return _RFD4_DETECT.is_rfd4_cif(filename)


def rfd4_load(filename, object="", mode="coordinates", labels=0, repo="", python="", _self=cmd,
              _automatic=False, _zoom=-1):
    """Inspect one annotated CIF in either launch mode.

    rfd4_load /path/input.cif [, object=motif] [, mode=coordinates] [, labels=0]
    """
    source = Path(str(filename).strip('"\'')).expanduser().resolve(strict=True)
    labels, _zoom = int(labels), int(_zoom)
    if mode not in _RFD4_MODES:
        raise ValueError("Mode must be " + ", ".join(_RFD4_MODES))
    base = source.with_suffix("") if source.suffix.lower() == ".gz" else source
    stem = re.sub(r"[^A-Za-z0-9_]", "_", base.stem)
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
    reserved = set(physical_keys)
    for centroid in data["centroids"]:
        marker = Atom()
        marker.name, marker.symbol, marker.resn, marker.chain = "ORI", "X", "ORI", ""
        marker.resi = str(len(centroid_keys) + 1)
        while (marker.chain, marker.resi, marker.resn, marker.name) in reserved:
            marker.resi = str(int(marker.resi) + 1)
        marker.segi, marker.hetatm, marker.q, marker.vdw = "RFD4_CENTROID", 1, 1.0, 1.5
        marker.coord = list(centroid["xyz"])
        marker.id = len(model.atom) + 1
        model.atom.append(marker)
        key = (marker.chain, marker.resi, marker.resn, marker.name)
        centroid_keys.append(key)
        reserved.add(key)
    _self.load_model(model, obj, zoom=0, discrete=0)
    # PyMOL may sort atoms and replace serial IDs during load. Match identities
    # before building masks so native annotation indices can never shift atoms.
    identities = {}
    observed = _self.iterate(obj, "indices[(chain, resi, resn, name)] = index", space={"indices": identities})
    if len(physical_keys) != len(data["atoms"]) or len(identities) != observed or identities.keys() != reserved:
        _self.delete(obj)
        raise RuntimeError("PyMOL did not preserve unique native atom identities; no annotation masks were applied.")
    for atom in data["atoms"]:
        atom["pymol_index"] = identities[atom["chain"], atom["resi"], atom["resn"], atom["name"]]
    data["centroid_indices"] = [identities[key] for key in centroid_keys]
    _RFD4_OBJECTS[obj] = data
    try:
        if not _automatic:
            _ensure_selections(obj, _self)
        rfd4_style(obj, mode, labels=labels, _self=_self)
    except Exception:
        _self.delete(obj)
        _RFD4_OBJECTS.pop(obj, None)
        raise
    if _zoom == 1 or (_zoom < 0 and first_object):
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
    names = [_selection(obj, suffix) for suffix in selections]
    group = obj + "_annotations"
    new_group = group not in _self.get_names("group_objects")
    _self.group(group, " ".join(names), action="add")
    if new_group:
        # A visible selection-only group consumes a blank cell in grid mode.
        _self.disable(group)
        _self.group(group, action="close")
    data.setdefault("selection_names", set()).update(names)
    data["coordinate_selections_created"] = True
    if not coordinates_only:
        data["selections_created"] = True


def _element_palette(_self, legacy=False):
    """Look up installed theme colors, or use PyMOL colors when run standalone."""
    theme = "foliage" if legacy else getattr(_self, "_pymolrc_theme", "")
    palettes = globals().get("_PYMOLRC_THEMES", {})
    side, backbone, ligand = palettes.get(theme, ("green", "forest", "violet"))
    return {"protein": side, "backbone": backbone, "ligand": ligand, "soft": theme == "foliage"}


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
    _self.select_list(selection, obj, indices, mode="index")
    _self.disable(selection)
    _self.group(obj + "_annotations", selection, action="add")
    data.setdefault("selection_names", set()).add(selection)
    _self.hide("everything", selection)
    _self.label(selection, "")
    if not data.get("show_centroids", True):
        print("[RFD4] Centroid markers hidden.")
        return
    _self.color("hotpink", selection)
    _self.set("sphere_color", "hotpink", selection)
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


_RFD4_MODES = {
    "coordinates": (), "elements": (), "legacy": (),
    "sequence": (("sequence_fixed", "teal"), ("sequence_free", "gray70")),
    "unindexed": (("unindexed", "violet"),),
    "hbond": (("hb_donor", "marine"), ("hb_acceptor", "salmon"), ("hb_both", "magenta")),
    "burial": (("buried", "forest"), ("partial", "orange"), ("exposed", "skyblue")),
}
_RFD4_METAL_SEL = globals().get("_METAL_SEL", "metals")
for _name, _rgb in {"rfd4_gold": (246, 187, 69), "rfd4_ochre": (184, 146, 74),
                    "rfd4_lilac": (187, 145, 245),
                    "rfd4_slate": (105, 177, 216)}.items():
    cmd.set_color(_name, list(_rgb), quiet=1)


def _rfd4_recolor(obj, _self=cmd):
    """Apply a molecular theme to physical atoms without changing any representations."""
    data = _RFD4_OBJECTS[obj]
    options = data.get("style", {})
    mode = options.get("mode", "coordinates")
    if mode not in {"coordinates", "elements", "legacy"} or (mode != "elements" and not int(options.get("element_colors", 1))):
        return
    legacy = mode == "legacy"
    internal_color = "orange" if legacy else options.get("internal_color", "rfd4_lilac")
    internal = _selection(obj, "internal")
    scoped = "model " + obj
    if data.get("centroid_indices"):
        markers = "+".join(map(str, data["centroid_indices"]))
        scoped += " and not index " + markers
    colors = _element_palette(_self, legacy)
    _self.color("atomic", scoped)
    protein_color = "cyan" if legacy else colors["protein"]
    _self.color(protein_color, scoped + " and polymer.protein and elem C")
    if mode == "elements":
        _self.color(colors["backbone"], scoped + " and polymer.protein and elem C and backbone")
    if legacy:
        _self.color(colors["ligand"], scoped + " and not polymer.protein and elem C")
    elif "_paint_ligand_carbons" in globals():
        _paint_ligand_carbons(scoped, colors["ligand"], _self, protein_color=protein_color)
    else:
        paint_component_carbons(scoped, colors["ligand"], _self, protein_color=protein_color)
    if colors["soft"]:
        for symbol in ("N", "O", "S", "P", "H"):
            color = "elem" + symbol
            if _self.get_color_index(color) >= 0:
                _self.color(color, scoped + " and elem " + symbol)
        if _self.get_color_index("elemZn") >= 0:
            _self.color("elemZn", scoped + " and " + _RFD4_METAL_SEL)
    if mode in {"coordinates", "legacy"}:
        _self.color(internal_color, internal + " and elem C")
    _shade_free_identity(obj, _self)
    if mode == "elements":
        if getattr(_self, "_pymolrc_gradient", False) and "_cartoon_gradient" in globals():
            _cartoon_gradient(scoped, _self._pymolrc_theme, _self)
            data["cartoon_gradient"] = True
        elif data.pop("cartoon_gradient", False):
            _self.unset("cartoon_color", scoped + " and polymer.protein and name CA")


def rfd4_style(object="", mode="coordinates", labels=0, context="cartoon", element_colors=1, _self=cmd,
               fixed_scale=0.16, fixed_color="rfd4_gold", internal_color="rfd4_lilac",
               markers="auto", overview_min=40, fixed_cartoon="rfd4_ochre"):
    """Show conditioning on cartoons. context=cartoon|lines|0, or mode=legacy for the original view."""
    obj = _annotated_object(object, _self)
    if mode not in _RFD4_MODES:
        raise ValueError("Mode must be " + ", ".join(_RFD4_MODES))
    scale = float(fixed_scale)
    if not 0 < scale < float("inf"):
        raise ValueError("fixed_scale must be a finite positive sphere scale.")
    for color in (fixed_color, internal_color, fixed_cartoon):
        if _self.get_color_index(color) < 0:
            raise ValueError("Unknown PyMOL color " + str(color))
    context_name = {"0": "off", "1": "lines"}.get(str(context), str(context).lower())
    if context_name not in {"cartoon", "lines", "off"}:
        raise ValueError("context must be cartoon, lines (or 1), or off (or 0).")
    overview_min = int(overview_min)
    if markers not in {"auto", "all"} or overview_min < 1:
        raise ValueError("markers must be auto or all, and overview_min must be a positive C-alpha count.")
    legacy = mode == "legacy"
    coordinate_view = mode in {"coordinates", "elements", "legacy"}
    chemical_colors = mode == "elements" or (coordinate_view and int(element_colors))
    _ensure_selections(obj, _self, coordinates_only=coordinate_view)
    data = _RFD4_OBJECTS[obj]
    data["style"] = dict(mode=mode, labels=labels, context=context, element_colors=element_colors,
                         fixed_scale=fixed_scale, fixed_color=fixed_color, internal_color=internal_color,
                         markers=markers, overview_min=overview_min, fixed_cartoon=fixed_cartoon)
    if legacy:
        context_name = "lines" if context_name == "cartoon" else context_name
        fixed_color, internal_color, scale = "cyan", "orange", 0.25

    absolute, internal, reference = (_selection(obj, name) for name in ("absolute", "internal", "reference"))
    targets = "(%s or %s)" % (absolute, internal)
    scoped = "model " + obj
    with_ca = "(byres (%s and polymer.protein and name CA))" % scoped
    incoming = "(byres (%s and name N and neighbor (%s and name C)))" % (with_ca, with_ca)
    outgoing = "(byres (%s and name C and neighbor (%s and name N)))" % (with_ca, with_ca)
    shown = targets if coordinate_view else scoped
    overview = (mode == "coordinates" and context_name == "cartoon" and markers == "auto"
                and _self.count_atoms(absolute + " and polymer.protein and name CA") >= overview_min)
    if overview:
        # A continuous, wholly fixed peptide already has a gold trace. Retain
        # beads on partial residues, termini, disconnected motifs and ligands.
        interior = incoming + " and " + outgoing
        mixed = "(byres (%s or %s))" % (reference, internal)
        shown = "(%s and not ((%s) and not %s))" % (targets, interior, mixed)
    visible = scoped if not coordinate_view or context_name != "off" else targets
    metal_scope = shown if legacy else visible
    metal_atoms = "(%s) and %s" % (metal_scope, _RFD4_METAL_SEL)

    # This viewer handles metal visibility explicitly in each context.
    show = getattr(_self, "_pymolrc_metal_show_state", {}).get("show", _self.show)
    _self.hide("everything", obj)
    _self.unset("cartoon_color", scoped)
    if data.pop("cartoon_gradient", False):
        _self.unset("cartoon_color", scoped + " and polymer.protein and name CA")
    _self.unset("sphere_color", scoped)
    _self.unset("sphere_scale", obj)
    _self.unset("sphere_scale", scoped)
    _self.set("sphere_scale", scale, "%s and not %s" % (scoped, _RFD4_METAL_SEL))
    _self.set("stick_radius", 0.13, obj)
    show("sticks", shown)
    show("nonbonded", shown)
    if context_name == "cartoon":
        show("cartoon", scoped + " and polymer")
    if coordinate_view and context_name != "off":
        lines = reference
        if context_name == "cartoon":
            # Keep ligands and fragments without a cartoon trace visible.
            lines += " and not (%s or %s)" % (incoming, outgoing)
        show("lines", lines)
        show("nonbonded", lines)
        _self.hide("spheres", reference)

    if not coordinate_view and int(element_colors):
        # Store elemental sphere faces before applying condition colors to the
        # atoms. A metal donor can have a native sphere and a donor-colored stick.
        _self.color("atomic", metal_atoms)
        _self.alter(metal_atoms, "s.sphere_color = color")
    if not chemical_colors:
        _self.color("gray60", obj)
    if chemical_colors:
        _rfd4_recolor(obj, _self)
    elif coordinate_view:
        _self.color(fixed_color, absolute)
        _self.color(internal_color, internal)
    else:
        for suffix, color in _RFD4_MODES[mode]:
            _self.color(color, _selection(obj, suffix))

    if mode in {"coordinates", "legacy"}:
        show("spheres", shown)
        if not legacy:
            # Representation colors encode the coordinate mask even on O/N/P.
            # Atom/stick colors still identify the element.
            for selection, color in ((absolute, fixed_color), (internal, internal_color)):
                beads = selection + " and not " + _RFD4_METAL_SEL if chemical_colors else selection
                _self.set("sphere_color", color, beads)
            _self.set("cartoon_color", fixed_cartoon, absolute)
            _self.set("cartoon_color", internal_color, internal)
            _self.set("cartoon_color", "rfd4_slate", reference)
    elif mode == "hbond":
        show("spheres", _selection(obj, "hb_donor") + " or " + _selection(obj, "hb_acceptor"))
    if not legacy and (chemical_colors or (not coordinate_view and not int(element_colors))):
        _self.set("sphere_color", -1, metal_atoms)  # follow the atom's element color
    show("spheres", metal_atoms)
    _self.label(obj, "")
    if int(labels):
        _self.label("(%s) and name CA" % visible, '"%s/%s%s" % (chain,resn,resi)')
    _style_centroids(obj, _self)
    if mode in {"coordinates", "legacy"}:
        atoms = data["atoms"]
        n_abs = sum(a["fixed"] for a in atoms)
        n_int = sum(a["internal"] for a in atoms)
        print("[RFD4] %d absolute (%s), %d internal (%s), %d reference (%s context)." %
              (n_abs, fixed_color, n_int, internal_color, len(atoms) - n_abs - n_int, context_name))
        if overview:
            print("[RFD4] Fully fixed stretches are cartoon-only. Use rfd4_style , coordinates, markers=all for every bead.")
    else:
        print("[RFD4] " + mode + "  " + ", ".join(key + "=" + color for key, color in _RFD4_MODES[mode]))

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
    print("  reference-only  %d retained atoms; slate cartoon in coordinate view (context=0 hides them)." %
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
RFD4 annotation inspection
  pymol prediction.cif                       ordinary PyMOL load
  pymolrfd4 prediction.cif annotated.cif     inspect annotated files in a mixed session
  rfd4_load input.cif, motif                 inspect one file explicitly
  rfd4_autoload on | off | status             routing for subsequent loads

  rfd4_style motif, coordinates              ochre absolute / lilac internal / blue reference cartoon
  rfd4_style motif, coordinates, markers=all show every target as a bead
  rfd4_style motif, coordinates, fixed_scale=0.12, fixed_color=yellow
  rfd4_style motif, coordinates, context=0   hide reference atoms (lines = atomic context)
  rfd4_style motif, sequence                 teal fixed / gray free identity
  rfd4_style motif, unindexed                violet unindexed atoms
  rfd4_style motif, hbond                    donor / acceptor / both
  rfd4_style motif, burial                   buried / intermediate / exposed
  rfd4_style motif, elements                 molecular palette and element-colored sticks
  rfd4_style motif, legacy                   original cyan/orange atom view with reference lines
  rfd4_centroids motif, 0                    hide centroid targets (1 shows them)
  rfd4_info motif                            source, counts and segment ranges
  rfd4_cache_info                            reader route and timing

Omit the name when one annotated object is visible, e.g. rfd4_style , sequence.
Cartoons follow the alpha-carbon mask. Small spheres mark individual conditioned
atoms. With 40 or more absolute alpha carbons, continuous fully fixed stretches
use cartoons alone. Adjust overview_min, fixed_scale, fixed_cartoon, fixed_color
or internal_color in rfd4_style. labels=1 shows residue IDs. The original CIF
holds the annotations; ORI centroids are display markers. See RFD4_QUICKSTART.md.
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
    fmt = args[3] if len(args) > 3 else kwargs.get("format", "")
    if not _RFD4_STATE["autoload"] or str(fmt).lower() not in {"", "cif", "mmcif"} or not _is_rfd4_cif(filename):
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


def rfd4_autoload(mode="status"):
    """Toggle automatic annotation inspection for subsequent loads in this session."""
    value = str(mode).strip().lower()
    if value in {"on", "off"}:
        _RFD4_STATE["autoload"] = value == "on"
    elif value != "status":
        raise ValueError("Use rfd4_autoload on, off, or status.")
    print("[RFD4] Automatic inspection " + ("on" if _RFD4_STATE["autoload"] else "off"))


# Explicitly running the standalone kit still enables it on vanilla PyMOL.
# The full startup config reserves interception for the detail launcher.
_RFD4_STATE.setdefault("autoload", os.environ.get("PYMOL_RFD4", "") == "1" or
                       not hasattr(cmd, "_pymolrc_originals"))
cmd.extend("rfd4_autoload", rfd4_autoload)
_automatic_load._rfd4_dispatch = True
_install_autoload()
if not _RFD4_STATE.get("atexit_registered"):
    atexit.register(_stop_native_worker)
    _RFD4_STATE["atexit_registered"] = True
