"""Named molecular palettes and the interactive PyMOL theme command."""
from pymol import cmd

# Sidechain carbon, backbone carbon, ligand carbon. Both green themes refer to
# the user's original color names so manual set_color overrides still work.
_FOLIAGE_CARBONS = ("sc05", "bb05", "lig_amethyst3")
_PYMOLRC_THEMES = {
    "classic": _FOLIAGE_CARBONS,
    "fjord": ("fjord_mint", "fjord_teal", "fjord_lilac"),
    "clay": ("clay_peach", "clay_terracotta", "clay_aqua"),
    "iris": ("iris_lavender", "iris_violet", "iris_jade"),
    "foliage": _FOLIAGE_CARBONS,
}
_SOFT_ELEMENTS = {"N": "elemN", "O": "elemO", "S": "elemS", "P": "elemP", "H": "elemH"}
if not hasattr(cmd, "_pymolrc_theme"):
    cmd._pymolrc_theme = "classic"


def _color_structure(selection, _self=cmd, carbons=None, palette=None):
    """Color atoms without changing the camera, representations or settings on other objects."""
    palette = palette or _self._pymolrc_theme
    sidechain, backbone, ligand = carbons or _PYMOLRC_THEMES[palette]
    scope = "(%s)" % selection
    _self.color("atomic", scope)
    _self.color(sidechain, scope + " and polymer and elem C")
    _self.color(backbone, scope + " and polymer.protein and backbone and elem C")
    _self.color(ligand, scope + " and organic and elem C")
    if palette == "foliage":
        for element, color in _SOFT_ELEMENTS.items():
            _self.color(color, scope + " and elem " + element)
        _self.color("elemZn", scope + " and (" + _METAL_SEL + ")")
    _self.unset("cartoon_color", selection)


def theme(name="list", _self=cmd):
    """theme [list|next|prev|classic|fjord|clay|iris|foliage]. Applies to open and future structures."""
    names = list(_PYMOLRC_THEMES)
    current = _self._pymolrc_theme
    name = str(name).lower().strip()
    if name in ("", "list"):
        for key, colors in _PYMOLRC_THEMES.items():
            elements = "soft" if key == "foliage" else "standard"
            print((" * " if key == current else "   ") + key + "  " + " / ".join(colors)
                  + "  (" + elements + " element colors)")
        return current
    if name in ("next", "prev"):
        step = 1 if name == "next" else -1
        name = names[(names.index(current) + step) % len(names)]
    if name not in _PYMOLRC_THEMES:
        raise ValueError("Choose " + ", ".join(names) + ", next, prev or list.")
    _self._pymolrc_theme = name
    annotations = getattr(_self, "_rfd4_viewer_state", {}).get("objects", {})
    for obj in _self.get_object_list():
        if obj in annotations:
            _rfd4_recolor(obj, _self)
        else:
            _color_structure(obj, _self)
            _color_catres(obj, _self, show_sticks=False)
    print("Theme " + name + "  (" + " / ".join(_PYMOLRC_THEMES[name]) + ")")
    return name


cmd.extend("theme", theme)
