"""Structure palettes, publication rendering and ball-and-stick views."""
from pymol import cmd


def publication_ray_trace(filename="publication_ray.png", width=3200, height=2400,
                          dpi=300, transparent=0, _self=cmd):
    # Apply temporary high-quality ray settings, save PNG, then restore.
    width = int(width)
    height = int(height)
    dpi = int(dpi)
    transparent = int(transparent)
    if not filename:
        filename = "publication_ray.png"
    if not str(filename).lower().endswith('.png'):
        filename = str(filename) + '.png'

    tracked = [
        'antialias', 'ray_trace_mode', 'ray_trace_disco_factor', 'ray_trace_gain',
        'ray_shadow', 'ray_texture', 'ray_opaque_background', 'stick_quality',
        'sphere_quality', 'cartoon_sampling'
    ]
    old = {}
    for key in tracked:
        old[key] = _self.get(key)

    try:
        _self.set('antialias', 4)
        _self.set('ray_trace_mode', 1)
        _self.set('ray_trace_disco_factor', 0)
        _self.set('ray_trace_gain', 0.12)
        _self.set('ray_shadow', 'off')
        _self.set('ray_texture', 0)
        _self.set('ray_opaque_background', 0 if transparent else 1)
        _self.set('stick_quality', 16)
        _self.set('sphere_quality', 3)
        _self.set('cartoon_sampling', 14)
        _self.ray(width, height)
        _self.png(filename, dpi=dpi, ray=0)
        print('publication_ray_trace: wrote %s (%dx%d, dpi=%d, transparent=%d)' % (
            filename, width, height, dpi, transparent
        ))
    finally:
        for key in tracked:
            _self.set(key, old[key])


def color_palette(selection="all", ligand_shell_radius=5.0, _self=cmd):
    """Align all objects to first, cycle foliage palettes by object,
    color ligands/elements, and show ligand-shell sidechain sticks."""
    sidechain_palette = ['sc10','sc09','sc08','sc07','sc06',
                         'sc05','sc04','sc03','sc02','sc01']
    backbone_palette = ['bb10','bb09','bb08','bb07','bb06',
                        'bb05','bb04','bb03','bb02','bb01']
    objects = _self.get_object_list(selection)
    if len(objects) == 0:
        return
    # Align all structures to the first one loaded
    target = objects[0]
    for obj in objects[1:]:
        try:
            _self.super(obj, target)
        except:
            pass
    # Carbon coloring:
    # - sidechain and backbone carbons get paired greens per object
    # - ligand carbons use amethyst for strong protein/ligand contrast
    for i, obj in enumerate(objects):
        sc_color = sidechain_palette[i % len(sidechain_palette)]
        bb_color = backbone_palette[i % len(backbone_palette)]
        _color_structure(obj, _self, (sc_color, bb_color, 'lig_amethyst3'), palette="foliage")
    # Metals as spheres (element/resn/atom-name) + catalytic (REMARK 666) sticks
    _self.show("spheres", "(%s) and (%s)" % (selection, _METAL_SEL))
    for _obj in objects:
        _color_catres(_obj, _self=_self, color="catres_hi")
    if _AF3_CLASH_FILTER_ENABLED:
        remove_clash_bonds(selection, _self=_self)
    if _AF3_LIGAND_CENTER_FILTER_ENABLED:
        remove_ligand_center_pseudobonds(selection, _self=_self)

    # Show ligand and nearby sidechains as sticks for fast active-site views.
    ligand_sel = '(%s) and organic and not (%s)' % (selection, _METAL_SEL)
    if _self.count_atoms(ligand_sel) > 0:
        _self.show('sticks', ligand_sel)
        _self.set('stick_radius', 0.25, ligand_sel)
        shell_sel = 'byres ((%s) and polymer.protein within %.2f of %s)' % (
            selection, float(ligand_shell_radius), ligand_sel)
        _self.select('_tmp_lig_shell', shell_sel)
        _self.hide('sticks', '_tmp_lig_shell and polymer.protein and resn GLY')
        _self.show('sticks', '_tmp_lig_shell and polymer.protein and not resn GLY and sidechain and not elem H')
        _self.set('stick_radius', 0.25, '_tmp_lig_shell and polymer.protein and not resn GLY and sidechain and not elem H')
        _self.show('sticks', '_tmp_lig_shell and polymer.protein and not resn GLY and name CA')
        _self.set('stick_radius', 0.25, '_tmp_lig_shell and polymer.protein and not resn GLY and name CA')
        _self.show('sticks', '_tmp_lig_shell and polymer.protein and resn PRO and name N')
        _self.set('stick_radius', 0.25, '_tmp_lig_shell and polymer.protein and resn PRO and name N')
        _self.delete('_tmp_lig_shell')
    _self.deselect()


def color_bb_rfdiffusion(selection="chain A", all_atom=0, backbone_only=0, _self=cmd):
    """Color `selection` N->C with the RFdiffusion gradient (dark blue -> navaho).

    By default only CARBON atoms are recolored, leaving non-carbon atoms untouched.
      selection       one or more chains / any selection (default 'chain A'; e.g. 'chain A+B')
      all_atom=1      recolor every element, not just carbons
      backbone_only=1 restrict to backbone atoms (sidechains keep their colors)
    Combos -> default: all carbons; backbone_only: backbone carbons;
    all_atom: every atom; all_atom+backbone_only: all backbone atoms.
    """
    palette = ("RFd_darkblue RFd_blue RFd_lightblue RFd_purple "
               "RFd_pink RFd_melon RFd_navaho")
    sel = "(%s)" % selection
    if not _rfd_bool(all_atom):
        sel += " and elem C"
    if _rfd_bool(backbone_only):
        sel += " and backbone"
    _self.spectrum("resi", palette, sel)


def color_bb_rfdiffusion4_proteina(selection="chain A", all_atom=0, backbone_only=0, _self=cmd):
    """Color `selection` N->C with the RFD4-Proteina Northern Nights gradient
    (teal -> blue -> violet -> purple), using the paper's molecular_spectrum stops.

    Same options as color_bb_rfdiffusion:
      selection       one or more chains / any selection (default 'chain A')
      all_atom=1      recolor every element, not just carbons
      backbone_only=1 restrict to backbone atoms (sidechains keep their colors)
    """
    palette = "rfd4p_teal rfd4p_blue rfd4p_violet rfd4p_purple"
    sel = "(%s)" % selection
    if not _rfd_bool(all_atom):
        sel += " and elem C"
    if _rfd_bool(backbone_only):
        sel += " and backbone"
    _self.spectrum("resi", palette, sel)


def gaussian_mode(selection="all", _self=cmd):
    """Quantum-chemistry / GaussView look: glossy ball-and-stick (every atom incl.
    H as a shiny sphere, metal larger), bond orders, and perspective. Colors are
    left as-is. Run `gaussian_off` to restore the normal look."""
    sel = "(%s)" % selection
    _self.show("sticks", sel)
    _self.show("spheres", sel)
    _self.set("stick_radius", 0.13, sel)
    _self.set("sphere_scale", 0.23, sel)
    _self.set("sphere_scale", 0.5, sel + " and metals")
    _self.set("specular", 1.0)
    _self.set("shininess", 55)
    _self.set("sphere_quality", 3)
    _self.set("stick_quality", 12)
    _self.set("valence", 1)
    _self.set("orthoscopic", 0)
    _self.rebuild()


def gaussian_off(selection="all", _self=cmd):
    """Undo gaussian_mode: restore matte / orthographic / normal sizes and re-apply
    the standard style."""
    sel = "(%s)" % selection
    _self.unset("shininess")
    _self.unset("sphere_quality")
    _self.unset("stick_quality")
    _self.set("specular", 0)
    _self.set("valence", 0)
    _self.set("orthoscopic", 1)
    _self.unset("sphere_scale", sel)
    _self.unset("stick_radius", sel)
    _self.set("sphere_scale", 0.7)
    _self.set("stick_radius", 0.25)
    _self.hide("spheres", sel + " and not metals")
    _self.rebuild()
    style_all(selection, _self=_self)


def gaussian_spin(toggle="on", _self=cmd):
    """Gentle continuous spin for the 'animating' effect (GUI only).
    `gaussian_spin on` starts a looping y-roll; `gaussian_spin off` stops it."""
    import pymol.util as _util
    if _rfd_bool(toggle):
        _self.mset("1 x120")
        _util.mroll(1, 120, 1)
        _self.mplay()
    else:
        _self.mstop()
        _self.mclear()


cmd.extend('color_palette', color_palette)
cmd.extend('color_bb_rfdiffusion', color_bb_rfdiffusion)
cmd.extend('color_bb_rfdiffusion4_proteina', color_bb_rfdiffusion4_proteina)
cmd.extend('gaussian_mode', gaussian_mode)
cmd.extend('gaussian_off', gaussian_off)
cmd.extend('gaussian_spin', gaussian_spin)
cmd.extend('publication_ray_trace', publication_ray_trace)
