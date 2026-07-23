# ============================================================================
# crystal_kit.pml  —  PyMOL crystallography helpers (standalone, run-able)
#
#   run ~/pymolrc/crystallography/crystal_kit.pml
#
# Commands (type `help_crystal` for full args):
#   xtal <path>          load crystals (+ sibling MTZ maps), split subunits,
#                        align all subunits to each other and to a design model
#   xtal_align           re-superpose subunits (super -> cealign fallback)
#   xtal_density <sel>   2Fo-Fc mesh (+ Fo-Fc green/red diff) around a selection
#   xtal_rms             RMSD summary: every subunit vs reference and vs design
#   xtal_catres          highlight the design's REMARK 666 catalytic residues
#   help_crystal         print this reference
#
# No dependency on ~/.pymolrc. Field conventions: 2Fo-Fc main map at ~1.0 sigma
# (blue mesh). Fo-Fc difference map at +/-3.0 sigma (green = positive/unmodeled,
# red = negative/overbuilt). Density levels are in sigma (normalize_ccp4_maps on).
# ============================================================================
python

from pymol import cmd, util
import os as _cx_os
import glob as _cx_glob
import re as _cx_re

# module-level registry populated by xtal()
_CX = {"crystals": [], "subunits": [], "maps": {}, "design": None, "ref": None}

# colors (carbon identity per group; non-carbon atoms keep CPK element colors).
# Palette: crystals = cool grey backbones + warm gold catalytic; design = foliage
# green backbone + aqua-cyan catalytic + amethyst ligand. Warm-vs-cool separates the
# crystal and design active sites when superposed; gold-vs-cyan is colorblind-safe.
cmd.set_color("cx_green",      [120, 172, 115])  # design protein carbons (foliage)
cmd.set_color("cx_amethyst",   [176, 126, 202])  # design ligand/cofactor carbons
cmd.set_color("cx_cat_design", [ 90, 202, 214])  # design catalytic carbons (aqua-cyan)
cmd.set_color("cx_cat_xtal",   [238, 210, 118])  # crystal catalytic carbons (goldenrod)
cmd.set_color("cx_catres",     [230,  40, 200])  # legacy magenta (still selectable)
cmd.set_color("elemZn",        [109, 124, 142])  # metal blue-grey (matches pymolrc)
cmd.set_color("cx_grey1", [203, 206, 211])       # crystal backbone -- light
cmd.set_color("cx_grey2", [157, 160, 166])       # crystal backbone -- mid
cmd.set_color("cx_grey3", [119, 122, 129])       # crystal backbone -- deep
_CX_XTAL_COLORS = ["cx_grey1", "cx_grey2", "cx_grey3"]

# Stick / metal sizing consistent with the user's pymolrc so kit objects match the
# rest of their look (normal-thickness sidechains, metal spheres at 0.7), rather
# than the thinner kit-specific sizes used before. Metals inherit the global
# sphere_scale exactly like `show_metals` in the pymolrc does.
cmd.set("stick_radius", 0.25)
cmd.set("sphere_scale", 0.7)
cmd.set("valence", 0)


def _cx_bool(v):
    return str(v).strip().lower() in ('1', 'true', 'on', 'yes', 't', 'y')


def _cx_short(base, idx):
    """A short, typeable object name from a crystal filename."""
    m = _cx_re.search(r'[Pp]in\d+', base)
    if m:
        return "xtal_" + m.group(0)
    return "xtal%d" % idx


def _cx_remark666(pdb_path):
    """Parse REMARK 666 MATCH MOTIF lines -> list of (chain, resi) strings."""
    out = []
    try:
        fh = open(_cx_os.path.expanduser(pdb_path))
    except Exception:
        return out
    for line in fh:
        if line.startswith("REMARK 666") and "MATCH MOTIF" in line:
            f = line.split()
            try:
                i = f.index("MOTIF")
                out.append((f[i + 1], f[i + 3]))   # chain, resi
            except (ValueError, IndexError):
                pass
    fh.close()
    return out


def _cx_super(mobile, target, _self=cmd):
    """super with a cealign fallback; returns (rmsd, n_atoms) or (None, 0)."""
    try:
        r = _self.super(mobile, target)
        if r and r[1] > 0:
            return (r[0], r[1])
    except Exception:
        pass
    try:
        r = _self.cealign(target, mobile)
        return (r.get("RMSD"), r.get("alignment_length", 0))
    except Exception:
        return (None, 0)


def _cx_style(_self=cmd):
    """Default look: crystals colored per-crystal, design in foliage green+amethyst.
    Stick/sphere sizes come from the global stick_radius (0.25) / sphere_scale (0.7)
    set at kit load, so it matches the rest of the pymolrc look."""
    _self.hide("everything", "all")
    for i, cname in enumerate(_CX["crystals"]):
        col = _CX_XTAL_COLORS[i % len(_CX_XTAL_COLORS)]
        subs = [s for s in _CX["subunits"] if s.startswith(cname + "_")] or [cname]
        for s in subs:
            _cx_style_object(s, col, _self)
    if _CX["design"]:
        d = _CX["design"]
        _self.show("cartoon", d)
        _self.show("sticks", "%s and not polymer and not solvent" % d)
        _self.color("cx_green", "%s and polymer and elem C" % d)
        _self.color("cx_amethyst", "%s and not polymer and elem C" % d)
        util.cnc(d)
        _self.color("elemZn", "%s and metals" % d)
    _self.set("cartoon_transparency", 0.0)


def _cx_style_object(s, col, _self=cmd):
    """Cartoon + ligand/hetero sticks + metal spheres on one object, in color `col`.
    Sizes inherit the global stick_radius/sphere_scale (pymolrc-consistent)."""
    _self.show("cartoon", s)
    _self.show("sticks", "%s and not polymer and not solvent and not metals" % s)
    _self.show("spheres", "%s and metals" % s)
    _self.color(col, "%s and elem C" % s)
    util.cnc("%s" % s)
    _self.color("elemZn", "%s and metals" % s)


def _cx_style_originals(_self=cmd):
    """Keep each untouched, un-aligned original crystal as a toggle-able reference:
    style it (per-crystal color, cartoon + ligand sticks + metal spheres) so it looks
    right when shown, then DISABLE it so it stays out of view until the user enables
    it, and move the originals to the top of the object panel."""
    tops = []
    for i, c in enumerate(_CX["crystals"]):
        col = _CX_XTAL_COLORS[i % len(_CX_XTAL_COLORS)]
        _cx_style_object(c, col, _self)
        _self.disable(c)
        tops.append(c)
    if tops:
        _self.order(" ".join(tops), location="top")


def _cx_work_sel():
    """The aligned working set (subunits + design), NOT the un-aligned originals — so
    `zoom` frames the superposed cluster instead of the spread-out reference copies."""
    base = _CX["subunits"] or _CX["crystals"]
    work = list(base) + ([_CX["design"]] if _CX["design"] else [])
    return ("(" + " or ".join(work) + ")") if work else "all"


def xtal(path, design="", split=1, align=1, to_design=1, catres=1,
         focus="all", surface_metals=1, _self=cmd):
    """Load crystal structures (+ their sibling MTZ maps), split subunits, and
    align everything. The everything-button.

      path      : a directory of crystals, a glob, or a single .pdb.
      design    : path to a design/reference model (its REMARK 666 catalytic
                  residues get highlighted); optional.
      split     : 1 = split each crystal into per-chain subunit objects (default).
      align     : 1 = superpose all subunits onto a reference (default).
      to_design : 1 = also superpose them onto the design model (default).
      catres    : 1 = highlight the design's REMARK 666 residues in the crystals.
      focus     : 'all' (default) or 'active' to zoom straight to the catalytic pocket.
      surface_metals : 1 = show all metals (default); 0 = hide metals outside the
                  active site (surface Zn are kept in by default -- they're useful).
    Loads <name> per crystal (+ <name>_map holding <name>_map.2fofc/.fofc), and
    <name>_<chain> per subunit. Type `help_crystal` for the full command set.
    """
    _self.set("normalize_ccp4_maps", 1)
    path = _cx_os.path.expanduser(path)
    if _cx_os.path.isdir(path):
        pdbs = sorted(_cx_glob.glob(_cx_os.path.join(path, "*.pdb")))
    elif ("*" in path or "?" in path):
        pdbs = sorted(_cx_glob.glob(path))
    else:
        pdbs = [path]
    # if the design model lives in the same folder, don't load it as a crystal too
    if design:
        _dpath = _cx_os.path.realpath(_cx_os.path.expanduser(design))
        pdbs = [p for p in pdbs if _cx_os.path.realpath(p) != _dpath]
    if not pdbs:
        print("xtal: no PDB files found at %s" % path); return
    _CX["crystals"] = []; _CX["subunits"] = []; _CX["maps"] = {}; _CX["ref"] = None
    for idx, p in enumerate(pdbs, 1):
        base = _cx_os.path.splitext(_cx_os.path.basename(p))[0]
        name = _cx_short(base, idx)
        while name in _CX["crystals"]:
            name += "b"
        _self.load(p, name)
        _CX["crystals"].append(name)
        mtz = _cx_os.path.splitext(p)[0] + ".mtz"
        if _cx_os.path.exists(mtz):
            mapobj = name + "_map"
            _self.load(mtz, mapobj)
            has = [n for n in _self.get_names("objects")
                   if n.startswith(mapobj + ".") and _self.get_type(n) == "object:map"]
            if has:
                _CX["maps"][name] = mapobj
            else:
                print("xtal: %s.mtz has no map coefficients (2fofc/fofc) — "
                      "looks like a data MTZ, not a refined map MTZ." % base)
        else:
            print("xtal: no sibling MTZ for %s (looked for %s)" % (name, _cx_os.path.basename(mtz)))
        if _cx_bool(split):
            for ch in _self.get_chains(name):
                if _self.count_atoms("%s and chain %s and polymer" % (name, ch)) < 20:
                    continue
                sub = "%s_%s" % (name, ch)
                # include this chain's active-site metals only -- NOT nearby protein
                # or modified residues from neighbouring chains (that double-counted
                # catalytic resi numbers, e.g. a neighbour's KCX 16).
                het = ("(%s and metals within 4 of (%s and chain %s and polymer))"
                       % (name, name, ch))
                _self.create(sub, "byres ((%s and chain %s) or %s)" % (name, ch, het))
                _CX["subunits"].append(sub)
            # keep the original, unsplit, un-aligned crystal as a reference object
            # (styled + disabled below via _cx_style_originals) instead of deleting it.
    if design:
        _self.load(_cx_os.path.expanduser(design), "design")
        _CX["design"] = "design"
        _CX["_cat666"] = _cx_remark666(design)
    _cx_style(_self)
    if _cx_bool(align):
        xtal_align(to_design=to_design, _self=_self, quiet=1)
    if _cx_bool(catres) and _CX["design"] and _CX.get("_cat666"):
        xtal_catres(_self=_self)
    if _CX["subunits"]:
        # originals survive as disabled, top-of-panel reference objects
        _cx_style_originals(_self)
    if not _cx_bool(surface_metals):
        _cx_hide_surface_metals(_self)
    if str(focus).lower() in ("active", "act", "site"):
        xtal_focus("active", _self=_self)
    else:
        _self.zoom(_cx_work_sel())
    n_sub = len(_CX["subunits"]) or len(_CX["crystals"])
    print("xtal: loaded %d crystal(s), %d subunit(s), maps for %d, design=%s. "
          "Try: xtal_rms | xtal_density <subunit> | xtal_catres | help_crystal"
          % (len(_CX["crystals"]), n_sub, len(_CX["maps"]),
             _CX["design"] or "none"))


def xtal_align(reference="", to_design=1, _self=cmd, quiet=0):
    """Superpose all crystal subunits onto a reference (default: first subunit),
    and optionally onto the design model. super -> cealign fallback.
      reference : object to align onto (default first subunit).
      to_design : 1 = also align the design onto the reference.
    """
    targets = _CX["subunits"] or _CX["crystals"]
    if not targets:
        print("xtal_align: nothing loaded — run `xtal <path>` first."); return
    ref = reference or targets[0]
    _CX["ref"] = ref
    rp = "(%s) and polymer" % ref
    for s in targets:
        if s == ref:
            continue
        _cx_super("(%s) and polymer" % s, rp, _self=_self)
    if _cx_bool(to_design) and _CX["design"]:
        # design shares the crystal sequence/numbering -> align by SEQUENCE so the
        # catalytic residues pair up correctly (super can pick a shifted register).
        d = _CX["design"]
        ok = False
        try:
            rr = _self.align("(%s) and polymer" % d, rp, cycles=5)
            ok = bool(rr) and rr[1] > 0
        except Exception:
            ok = False
        if not ok:
            _cx_super("(%s) and polymer" % d, rp, _self=_self)
    if not quiet:
        print("xtal_align: aligned %d subunit(s) onto %s%s"
              % (len(targets) - 1, ref, " (+ design)" if _CX["design"] else ""))


def xtal_density(sel, level=1.0, carve=1.8, diff=1, dlevel=3.0, surface=0, _self=cmd):
    """Show electron density around `sel` from the corresponding crystal's MTZ.

    Draws the 2Fo-Fc map (blue) at `level` sigma, and (if diff=1) the Fo-Fc
    difference map at +dlevel (green) / -dlevel (red). Density tracks the atoms
    even after alignment (the map is matrix-copied to sel's frame).

      sel     : a subunit object (e.g. xtal_Pin5_A) or any selection within one.
      level   : 2Fo-Fc contour, sigma. TIGHTER = raise (1.5); LARGER = lower (0.8).
      carve   : radius (A) of density kept around the atoms. TIGHTER = lower (1.2);
                LARGER context = raise (2.5).
      diff    : 1 = also draw the +/- dlevel Fo-Fc difference map.
      dlevel  : difference-map contour sigma (default 3.0).
      surface : 1 = solid transparent isosurface instead of mesh.
    """
    parent = None
    for c in _CX["crystals"]:
        if sel.startswith(c):
            parent = c; break
    if parent is None:
        # sel might be a subunit object name -> derive parent crystal
        for c in _CX["crystals"]:
            for s in _CX["subunits"]:
                if s.startswith(c + "_") and sel.startswith(s):
                    parent = c; break
    if parent is None or parent not in _CX["maps"]:
        print("xtal_density: no map found for '%s'. Loaded maps: %s"
              % (sel, list(_CX["maps"].keys()))); return
    mapobj = _CX["maps"][parent]
    selobj = sel.split()[0]
    tag = "dens_" + selobj
    # building map/mesh objects must NOT steal the camera: copying the full MTZ
    # map makes a big grid object and auto_zoom would jump the view out to it.
    _az = _self.get("auto_zoom")
    _self.set("auto_zoom", 0)
    # copy the crystal's maps into sel's aligned frame
    for suf in ("2fofc", "fofc"):
        src = "%s.%s" % (mapobj, suf)
        if src not in _self.get_names("all"):
            continue
        cp = "%s_%s" % (tag, suf)
        _self.delete(cp)
        _self.copy(cp, src)
        _self.matrix_copy(selobj, cp)
    _self.delete(tag + "_2fofc_msh"); _self.delete(tag + "_fofc_pos"); _self.delete(tag + "_fofc_neg")
    m2 = "%s_2fofc" % tag
    if m2 in _self.get_names("all"):
        nm = tag + "_2fofc_msh"
        if _cx_bool(surface):
            _self.isosurface(nm, m2, float(level), sel, carve=float(carve))
            _self.set("transparency", 0.5, nm)
        else:
            _self.isomesh(nm, m2, float(level), sel, carve=float(carve))
        _self.color("skyblue", nm); _self.set("mesh_width", 0.5)
    if _cx_bool(diff):
        mf = "%s_fofc" % tag
        if mf in _self.get_names("all"):
            _self.isomesh(tag + "_fofc_pos", mf,  float(dlevel), sel, carve=float(carve))
            _self.isomesh(tag + "_fofc_neg", mf, -float(dlevel), sel, carve=float(carve))
            _self.color("green", tag + "_fofc_pos")
            _self.color("red",   tag + "_fofc_neg")
    _self.set("auto_zoom", _az)
    print("xtal_density: %s  2Fo-Fc @ %.2f sigma (blue), Fo-Fc @ +/-%.1f (green/red), carve %.1f A"
          % (sel, float(level), float(dlevel), float(carve)))


def xtal_rms(reference="", design="", _self=cmd):
    """Print a RMSD summary: every subunit vs the reference subunit and vs the
    design model, plus per-crystal info (resolution, chains, ligand)."""
    targets = _CX["subunits"] or _CX["crystals"]
    if not targets:
        print("xtal_rms: nothing loaded — run `xtal <path>` first."); return
    ref = reference or _CX["ref"] or targets[0]
    dsn = design or _CX["design"]
    print("=== RMSD summary (Angstroms) ===")
    print("%-16s %12s %12s" % ("subunit", "vs %s" % ref, "vs design"))
    for s in targets:
        r1 = "-" if s == ref else _cx_rms_ca(s, ref, _self)
        r2 = _cx_rms_ca(s, dsn, _self) if dsn else "-"
        print("%-16s %12s %12s" % (s, r1, r2))
    print("=== per-crystal ===")
    for c in _CX["crystals"]:
        subs = [s for s in _CX["subunits"] if s.startswith(c + "_")]
        sel = "(" + " or ".join(subs) + ")" if subs else c
        zn = _self.count_atoms("%s and resn ZN" % sel)
        lig = "YYE" if _self.count_atoms("%s and resn YYE" % sel) else "none"
        print("  %-14s subunits=%d  Zn=%d  ligand=%s" % (c, len(subs) or 1, zn, lig))


def _cx_rms_ca(a, b, _self):
    # RMSD via a throwaway super copy: chain-agnostic, and never moves the real
    # objects. super does its own structural matching so different chain IDs
    # (NCS copies A/B/C/D) pair up correctly.
    tmp = "_cx_rmstmp"
    _self.delete(tmp)
    try:
        _self.create(tmp, "(%s) and polymer" % a)
        r = _cx_super(tmp, "(%s) and polymer" % b, _self=_self)
        return ("%.3f" % r[0]) if r[0] is not None else "n/a"
    except Exception:
        return "n/a"
    finally:
        _self.delete(tmp)


def xtal_catres(design="", color="cx_cat_xtal", design_color="cx_cat_design",
                rep="sticks", label=0, _self=cmd):
    """Highlight the REMARK 666 catalytic residues in every crystal subunit AND the
    design -- but the design in its OWN color, so its active site reads apart from the
    crystals' when they superpose. Carbons recolored; heteroatoms stay CPK.
      color        : crystal catalytic carbons (default goldenrod cx_cat_xtal).
      design_color : design catalytic carbons (default aqua-cyan cx_cat_design).
      rep          : representation for the catalytic sidechains (default sticks).
      label        : 1 = add one-letter+resi labels.
    """
    cats = _cx_remark666(design) if design else _CX.get("_cat666")
    if not cats:
        print("xtal_catres: no REMARK 666 residues (pass design=<path>, or run `xtal` with a design)."); return
    resi_or = "+".join(sorted(set(r for _, r in cats), key=lambda x: int(x)))
    objs = list(_CX["subunits"]) + ([_CX["design"]] if _CX["design"] else [])
    n = 0
    for o in objs:
        s = "%s and resi %s and not (metals or solvent)" % (o, resi_or)
        if _self.count_atoms(s) == 0:
            continue
        col = design_color if o == _CX["design"] else color
        _self.show(rep, s)
        _self.color(col, "%s and elem C" % s)
        util.cnc(s)
        n += 1
        if _cx_bool(label):
            _self.label("%s and name CA" % s, '"%s%s" % (resn, resi)')
    print("xtal_catres: highlighted residues %s in %d object(s) [crystals=%s, design=%s]"
          % (resi_or, n, color, design_color))


def _cx_active_sel():
    """Selection string for the catalytic pocket across all loaded objects."""
    cats = _CX.get("_cat666")
    objs = list(_CX["subunits"]) + ([_CX["design"]] if _CX["design"] else [])
    if not cats or not objs:
        return None
    resi_or = "+".join(sorted(set(r for _, r in cats), key=lambda x: int(x)))
    # exclude metals/waters that happen to share a catalytic resi number (e.g. a Zn
    # numbered 16 collides with the catalytic Lys/KCX 16).
    return "((" + " or ".join(objs) + ") and resi " + resi_or + " and not (metals or solvent))"


def _cx_hide_surface_metals(_self=cmd):
    """Hide metals that are not in the catalytic pocket (keep active-site metals)."""
    act = _cx_active_sel()
    if not act:
        return
    _self.hide("everything", "(metals) and not ((metals) within 6 of (%s))" % act)


def xtal_focus(target="active", _self=cmd):
    """Zoom the view. `target=active` -> the catalytic pocket (default); `all` ->
    everything; or pass any selection string."""
    t = str(target).lower()
    if t in ("active", "act", "site"):
        cats = _CX.get("_cat666")
        ref = _CX.get("ref") or (_CX["subunits"][0] if _CX["subunits"] else None)
        if not cats or not ref:
            print("xtal_focus: no catalytic residues known (run `xtal` with a design)."); return
        resi_or = "+".join(sorted(set(r for _, r in cats), key=lambda x: int(x)))
        _self.zoom("%s and resi %s and not (metals or solvent)" % (ref, resi_or), 4)
    elif t == "all":
        # frame the aligned working set, not the un-aligned reference originals
        _self.zoom(_cx_work_sel())
    else:
        _self.zoom(target, 5)
    print("xtal_focus: %s" % target)


def help_crystal(_self=cmd):
    print("""
=== crystallography kit (crystal_kit.pml) ===
Conventions: 2Fo-Fc = main map (blue mesh, ~1.0 sigma); Fo-Fc = difference map
(green = +3 sigma unmodeled / red = -3 sigma overbuilt). Levels are in sigma.

xtal <path> [, design=] [, split=1] [, align=1] [, to_design=1] [, catres=1]
            [, focus=all] [, surface_metals=1]
    Load crystals + sibling MTZ maps, split into per-chain subunits, align all
    to a reference and to the design, style, highlight catalytic residues.
    MANDATORY: path (dir / glob / .pdb).  Everything else optional.
    focus=active zooms to the catalytic pocket; surface_metals=0 hides non-active Zn.

xtal_focus [, target=active]
    Zoom the view: active = catalytic pocket (default), all = everything, or a selection.

xtal_align [, reference=] [, to_design=1]
    Re-superpose all subunits onto a reference (default first subunit).

xtal_density <sel> [, level=1.0] [, carve=1.8] [, diff=1] [, dlevel=3.0] [, surface=0]
    2Fo-Fc mesh + Fo-Fc green/red difference around a subunit/selection.
    MANDATORY: sel.  TIGHTER: raise level / lower carve.  LARGER: the reverse.

xtal_rms [, reference=] [, design=]
    RMSD of every subunit vs reference and vs design, + per-crystal info.

xtal_catres [, design=<path>] [, color=cx_cat_xtal] [, design_color=cx_cat_design] [, rep=sticks] [, label=0]
    Highlight REMARK 666 catalytic residues: crystals gold, design cyan (distinct).

help_crystal   — this text.
""")


cmd.extend("xtal", xtal)
cmd.extend("xtal_align", xtal_align)
cmd.extend("xtal_density", xtal_density)
cmd.extend("xtal_rms", xtal_rms)
cmd.extend("xtal_catres", xtal_catres)
cmd.extend("xtal_focus", xtal_focus)
cmd.extend("help_crystal", help_crystal)

python end
