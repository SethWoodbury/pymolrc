"""Named molecular palettes and the interactive PyMOL theme command."""
from collections import defaultdict
import os
from pymol import cmd

# Sidechain carbon, backbone carbon, ligand carbon. Both green themes refer to
# the user's original color names so manual set_color overrides still work.
_FOLIAGE_CARBONS = ("sc05", "bb05", "lig_amethyst3")
_PYMOLRC_THEMES = {
    "mineral": ("mineral_taupe", "mineral_stone", "mineral_lagoon"),
    "fern": ("fern_sage", "fern_leaf", "fern_orchid"),
    "classic": _FOLIAGE_CARBONS,
    "fjord": ("fjord_mint", "fjord_teal", "fjord_lilac"),
    "clay": ("clay_peach", "clay_terracotta", "clay_aqua"),
    "iris": ("iris_lavender", "iris_violet", "iris_jade"),
    "foliage": _FOLIAGE_CARBONS,
}
_CARTOON_RAMPS = {
    "mineral": ("0x786F62", "mineral_stone", "0xA09886"),
    "fern": ("0x456950", "fern_leaf", "0x77977A"),
}
_SOFT_ELEMENTS = {"N": "elemN", "O": "elemO", "S": "elemS", "P": "elemP", "H": "elemH"}
if not hasattr(cmd, "_pymolrc_theme"):
    _requested_theme = os.environ.get("PYMOLRC_THEME", "mineral").lower().strip()
    if _requested_theme not in _PYMOLRC_THEMES:
        print("[pymolrc] Unknown PYMOLRC_THEME %r. Using mineral." % _requested_theme)
        _requested_theme = "mineral"
    cmd._pymolrc_theme = _requested_theme
if not hasattr(cmd, "_pymolrc_gradient"):
    cmd._pymolrc_gradient = os.environ.get("PYMOLRC_GRADIENT", "off").lower().strip() in ("on", "1")
if not hasattr(cmd, "_pymolrc_ligand_types"):
    cmd._pymolrc_ligand_types = {}
    cmd._pymolrc_split_ligands = False


def _ligand_accents(primary):
    # Put the current theme's ligand first. Repeat colors only after the cycle.
    return tuple(dict.fromkeys((primary, "fern_orchid", "mineral_lagoon", "ligand_olive")))


def _annotation_groups(_self):
    """Recognize saved RFD4 views even after a PSE has lost its Python metadata."""
    selections = set(_self.get_names("selections"))
    parents = (g.removesuffix("_annotations") for g in _self.get_names("group_objects") if g.endswith("_annotations"))
    return {p for p in parents if all(p + suffix in selections for suffix in ("_absolute", "_internal", "_reference"))}


def _paint_ligand_carbons(scope, primary, _self, uniform=False, protein_color=None):
    split = _self._pymolrc_split_ligands and not uniform
    paint_component_carbons(scope, primary, _self, protein_color,
                            _self._pymolrc_ligand_types if split else None, _ligand_accents(primary))


def _cartoon_gradient(scope, palette, _self):
    """Apply a residue-rank gradient to the cartoon while preserving atom colors."""
    sidechain, backbone, _ = _PYMOLRC_THEMES[palette]
    stops = [_self.get_color_tuple(c) for c in _CARTOON_RAMPS.get(palette, (backbone, sidechain))]
    ca = "(%s) and polymer.protein and name CA and elem C" % scope
    chains = defaultdict(dict)
    _self.iterate(ca, "chains[model, segi, chain].setdefault((resv, resi), []).append(index)", space={"chains": chains})
    colors = {}
    for (model, _segi, _chain), residues in chains.items():
        for rank, residue in enumerate(sorted(residues)):
            fraction = rank / (len(residues) - 1) if len(residues) > 1 else 0.5
            position = fraction * (len(stops) - 1)
            i = min(int(position), len(stops) - 2)
            weight = position - i
            rgb = tuple(round(255 * (a + (b - a) * weight)) for a, b in zip(stops[i], stops[i + 1]))
            color = _self.get_color_index("0x%02x%02x%02x" % rgb)
            for index in residues[residue]:
                colors[model, index] = color
    if colors:
        _self.alter(ca, "s.cartoon_color = colors[model, index]", space={"colors": colors})
        _self.recolor(ca)


def ligand_colors(mode="list", _self=cmd):
    """ligand_colors uniform|by_type|list. Applies to loaded and future organic ligand carbons."""
    mode = str(mode).lower().strip()
    if mode not in ("list", "uniform", "by_type"):
        raise ValueError("Choose uniform, by_type or list.")
    if mode != "list":
        _self._pymolrc_split_ligands = mode == "by_type"
    primary = _PYMOLRC_THEMES[_self._pymolrc_theme][2]
    if mode != "list":
        annotations = getattr(_self, "_rfd4_viewer_state", {}).get("objects", {})
        restored = _annotation_groups(_self) - annotations.keys()
        for obj in _self.get_object_list():
            if obj in restored:
                continue
            scope = "model " + obj
            if obj in annotations:
                style = annotations[obj].get("style", {})
                view = style.get("mode", "coordinates")
                if view == "coordinates" and int(style.get("element_colors", 1)):
                    scope += " and not " + _selection(obj, "internal")
                elif view != "elements":
                    continue
            _paint_ligand_carbons(scope, primary, _self)
    if _self._pymolrc_split_ligands:
        accents = _ligand_accents(primary)
        print("Ligand colors by residue name. The %d accents cycle in assignment order." % len(accents))
        for name, i in sorted(_self._pymolrc_ligand_types.items()):
            print("  " + name + "  " + accents[i % len(accents)])
    else:
        print("Uniform ligand carbon  " + primary)


def _color_structure(selection, _self=cmd, carbons=None, palette=None):
    """Color atoms without changing the camera, representations or settings on other objects."""
    palette = palette or _self._pymolrc_theme
    sidechain, backbone, ligand = carbons or _PYMOLRC_THEMES[palette]
    scope = "(%s)" % selection
    _self.color("atomic", scope)
    _self.color(sidechain, scope + " and polymer and elem C")
    _self.color(backbone, scope + " and polymer.protein and backbone and elem C")
    _paint_ligand_carbons(scope, ligand, _self, uniform=carbons is not None, protein_color=sidechain)
    if palette == "foliage":
        for element, color in _SOFT_ELEMENTS.items():
            _self.color(color, scope + " and elem " + element)
        _self.color("elemZn", scope + " and (" + _METAL_SEL + ")")
    _self.unset("cartoon_color", selection)
    if _self._pymolrc_gradient:
        if carbons is None:
            _cartoon_gradient(scope, palette, _self)
        else:
            _self.unset("cartoon_color", scope + " and polymer.protein and name CA")


def theme(name="list", _self=cmd, gradient=""):
    """theme [list|next|prev|name] [, gradient=on|off]. Applies to open and future structures."""
    names = list(_PYMOLRC_THEMES)
    current = _self._pymolrc_theme
    name = str(name).lower().strip()
    if name in ("", "list") and gradient == "":
        for key, colors in _PYMOLRC_THEMES.items():
            elements = "soft" if key == "foliage" else "standard"
            print((" * " if key == current else "   ") + key + "  " + " / ".join(colors)
                  + "  (" + elements + " element colors)")
        print("Cartoon gradient " + ("on" if _self._pymolrc_gradient else "off"))
        return current
    if name in ("", "list"):
        name = current
    if name in ("next", "prev"):
        step = 1 if name == "next" else -1
        name = names[(names.index(current) + step) % len(names)]
    if name not in _PYMOLRC_THEMES:
        raise ValueError("Choose " + ", ".join(names) + ", next, prev or list.")
    previous_gradient = _self._pymolrc_gradient
    if gradient != "":
        value = str(gradient).lower().strip()
        if value not in ("on", "off", "1", "0"):
            raise ValueError("Use gradient=on or gradient=off.")
        _self._pymolrc_gradient = value in ("on", "1")
    _self._pymolrc_theme = name
    annotations = getattr(_self, "_rfd4_viewer_state", {}).get("objects", {})
    restored = _annotation_groups(_self) - annotations.keys()
    for obj in _self.get_object_list():
        if obj in annotations:
            _rfd4_recolor(obj, _self)
        elif obj not in restored:
            if previous_gradient and not _self._pymolrc_gradient:
                _self.unset("cartoon_color", "model %s and polymer.protein and name CA" % obj)
            _color_structure(obj, _self)
            _color_catres(obj, _self, show_sticks=False)
    print("Theme " + name + "  (" + " / ".join(_PYMOLRC_THEMES[name]) + ")")
    return name


cmd.extend("theme", theme)
cmd.extend("ligand_colors", ligand_colors)
