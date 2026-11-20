"""Structure browsing, load-time styling and PyMOL command hooks."""
import gzip
import math
import os
from pymol import cmd

_AF3_STRUCTURE_NAV_INDEX = -1

def structure_step(step=1, selection="all", _self=cmd):
    # Cycle through loaded objects and show one structure at a time.
    global _AF3_STRUCTURE_NAV_INDEX
    objects = _self.get_object_list(selection)
    if not objects:
        print('structure_step: no objects loaded.')
        return
    enabled = set(_self.get_names('public_objects', enabled_only=1))
    if _AF3_STRUCTURE_NAV_INDEX < 0 or _AF3_STRUCTURE_NAV_INDEX >= len(objects):
        _AF3_STRUCTURE_NAV_INDEX = 0
        for i, obj in enumerate(objects):
            if obj in enabled:
                _AF3_STRUCTURE_NAV_INDEX = i
                break
    _AF3_STRUCTURE_NAV_INDEX = (_AF3_STRUCTURE_NAV_INDEX + int(step)) % len(objects)
    target = objects[_AF3_STRUCTURE_NAV_INDEX]
    if _AUTOSTYLE and target not in _STYLED and target in _KNOWN_OBJS:
        _style_one(target, _self=_self)   # lazy styling: style on first view (fast browsing)
    _self.disable('all')
    _self.enable(target)
    # camera intentionally preserved (no orient/center) so pgup/pgdn don't reset the view
    print('structure_step: %s (%d/%d)' % (target, _AF3_STRUCTURE_NAV_INDEX + 1, len(objects)))
    return target


def structure_next(selection="all", _self=cmd):
    return structure_step(1, selection, _self=_self)


def structure_prev(selection="all", _self=cmd):
    return structure_step(-1, selection, _self=_self)


def structure_show_all(selection="all", _self=cmd):
    if _AUTOSTYLE:
        for obj in _self.get_object_list(selection):
            if obj not in _STYLED and obj in _KNOWN_OBJS:
                _style_one(obj, _self=_self)
    _self.enable(selection)
    print('structure_show_all: enabled all objects.')
    return _self.count_atoms(selection)


def bind_structure_navigation_keys(_self=cmd):
    # PgUp/PgDn (and up/down) browse structures: walk a lazy-load list if active,
    # else cycle already-loaded objects. Camera preserved (no orient).
    _self.set_key('pgup', 'browse_prev')
    _self.set_key('pgdn', 'browse_next')
    _self.set_key('up', 'browse_prev')
    _self.set_key('down', 'browse_next')
    return 1


# ---------------------------------------------------------------------------
#  Metals, REMARK 666 catalytic residues, sequence, bg toggle (added 2026-06-12)
# ---------------------------------------------------------------------------

# Distinct highlight color for catalytic (REMARK 666) residue carbons.
# Change anytime, e.g.:  set_color catres_hi = [250,199,44]   (or: color yellow, catres)
cmd.set_color('catres_hi', [255, 150, 40])

# Metal detection by ELEMENT and/or RESIDUE NAME and/or ATOM NAME.
# The atom-name clause is guarded with "not polymer.protein" because atom name
# "CA" is the protein alpha-carbon, NOT calcium (so naive name matching would
# sphere every residue). Element/resn matching is unambiguous.
_METAL_ELEMS = ('Zn+Fe+Mg+Mn+Cu+Co+Ni+Ca+Na+K+Mo+V+W+Cd+Hg+Pt+Pd+Au+Ag+Al+'
                'Ba+Sr+Li+Cs+Rb+Cr')
_METAL_RESNS = ('ZN+FE+FE2+MG+MN+MN3+CU+CU1+CO+3CO+NI+CA+NA+K+MO+CD+HG+PT+PD+'
                'AU+AG+AL+BA+SR+LI+CS+RB+CR+V+W')
_METAL_NAMES = ('ZN+FE+MG+MN+CU+CO+NI+NA+K+MO+CD+HG+PT+PD+AU+AG+AL+BA+SR+LI+'
                'CS+RB+CR')
# Organic elements. An atom whose NAME collides with a metal symbol -- "CD"
# (Lys/Arg/Pro/Glu/Gln or carboxy-Lys KCX delta-carbon), "CO" (carbonyl C),
# "NI", "SR", ... -- but whose ELEMENT is one of these is NOT a metal. It is a
# side-chain/backbone atom of a non-standard residue (KCX, MSE, ligands) that
# PyMOL flags `not polymer`, so the `not polymer` gate below does not protect it
# the way it protects the alpha-carbon "CA". Excluding organic-element atoms from
# the name-based clause keeps e.g. the KCX CD carbon from being drawn as a grey
# metal sphere, while genuine ions (correct metal element) still match.
_METAL_ORGANIC = 'C+N+O+H+S+P+Se+B+F+Cl+Br+I'
# Gate on `not polymer` FIRST (a fast per-atom flag): metals/ions are never part of a
# protein/nucleic polymer chain. This shrinks the expensive elem/resn/name token match
# from every atom to just the handful of hetero atoms (~60x faster styling on load),
# and it makes the alpha-carbon "CA" safe automatically (CA is always polymer).
_METAL_SEL = ('(not polymer) and ((elem %s) or (resn %s) or '
              '((name %s) and not elem %s))'
              % (_METAL_ELEMS, _METAL_RESNS, _METAL_NAMES, _METAL_ORGANIC))

_CATRES_BY_OBJ = {}
_STYLED = set()
_AUTOSTYLE = True
_AUTOALIGN = False
_AUTOSOLO = False

_AA3TO1 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q',
    'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
    'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
    'TYR': 'Y', 'VAL': 'V', 'MSE': 'M', 'SEC': 'U', 'PYL': 'O', 'HID': 'H',
    'HIE': 'H', 'HIP': 'H', 'HSD': 'H', 'HSE': 'H', 'HSP': 'H', 'CYX': 'C',
    'CYM': 'C', 'LYN': 'K', 'ASH': 'D', 'GLH': 'E',
}


def show_metals(selection="all", _self=cmd):
    """Show every metal (by element / residue name / atom name) as a sphere.
    Sphere size comes from the global `sphere_scale` (0.7) — no per-selection `set`
    here, which keeps this O(1)-ish and cheap during big multi-file loads."""
    sel = '(%s) and %s' % (selection, _METAL_SEL)
    _self.show('spheres', sel)
    color = "elemZn" if _self._pymolrc_theme == "foliage" else "atomic"
    _self.color(color, sel)


# PMLs often hide everything after load and then show sticks/lines. Isolated
# ions have no bonds to draw, so restore their spheres when that representation
# is requested. Do not wrap hide/disable or change colors, coordinates or bonds.
# Keep original functions on cmd so reloading this section cannot nest hooks.
if not hasattr(cmd, '_pymolrc_metal_show_state'):
    cmd._pymolrc_metal_show_state = {'show': cmd.show, 'show_as': cmd.show_as}
_METAL_SHOW_STATE = cmd._pymolrc_metal_show_state


def _pymolrc_add_metal_spheres(representation, selection, _self=cmd):
    if not _AUTOSTYLE:
        return
    from pymol import selector, viewing
    # Match PyMOL's one-argument selection shorthand and native abbreviations.
    if not selection and (representation in ('', 'all') or
                          '(' in representation or '/' in representation):
        selection, representation = representation, 'wire'
    repmask = viewing._rep_to_repmask(representation)
    if not repmask & (viewing.repmasks['sticks'] | viewing.repmasks['lines'] |
                      viewing.repmasks['nonbonded'] | viewing.repmasks['nb_spheres']):
        return
    selected = selector.process(selection) or 'all'
    _METAL_SHOW_STATE['show']('spheres', '(%s) and (%s)' % (selected, _METAL_SEL), _self=_self)


def _pymolrc_show_hook(representation='wire', selection='', *, _self=cmd):
    result = _METAL_SHOW_STATE['show'](representation, selection, _self=_self)
    _pymolrc_add_metal_spheres(representation, selection, _self=_self)
    return result


def _pymolrc_show_as_hook(representation='wire', selection='', *, _self=cmd):
    result = _METAL_SHOW_STATE['show_as'](representation, selection, _self=_self)
    _pymolrc_add_metal_spheres(representation, selection, _self=_self)
    return result


cmd.show = _pymolrc_show_hook
cmd.show_as = _pymolrc_show_as_hook
cmd.keyword['show'][0] = cmd.show
cmd.keyword['as'][0] = cmd.show_as


_REMARK666_CACHE = {}   # path -> ((mtime, size), [(chain, resi), ...]); skips re-reading on reload


def _parse_remark666(path):
    """[(chain, resi), ...] of catalytic MATCH MOTIF residues from a PDB.
    REMARK 666 MATCH TEMPLATE <c> <n3> <r> MATCH MOTIF <chain> <n3> <resi> ...
    -> chain = field[9], resi = field[11]. Cheap: stops at first coordinate line.
    Cached by (mtime, size) so re-loading the same files (random_shuffle re-runs) is a stat,
    not a re-read of every header (matters on NFS)."""
    res = []
    # REMARK 666 is a PDB record, never CIF/PSE metadata. Avoid scanning entire
    # structures or binary sessions for it on every first load.
    if isinstance(path, str) and path.lower().removesuffix('.gz').endswith(('.cif', '.mmcif', '.pse', '.psw', '.pze')):
        return res
    if not path or not isinstance(path, str) or not os.path.isfile(path):
        return res
    key = None
    try:
        st = os.stat(path)
        key = (st.st_mtime, st.st_size)
        hit = _REMARK666_CACHE.get(path)
        if hit is not None and hit[0] == key:
            return hit[1]
    except OSError:
        return res
    try:
        opener = gzip.open if path.endswith('.gz') else open
        with opener(path, 'rt', errors='ignore') as fh:
            for line in fh:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    break
                if not line.startswith('REMARK 666 MATCH'):
                    continue
                sp = line.split()
                if len(sp) >= 12:
                    res.append((sp[9], sp[11]))
    except Exception:
        pass
    if key is not None:
        _REMARK666_CACHE[path] = (key, res)
    return res


def _catres_obj_sel(obj, cats):
    parts = ['(%s and chain %s and resi %s)' % (obj, ch, ri) for ch, ri in cats]
    return ' or '.join(parts) if parts else 'none'


def _catres_iter(_self=cmd):
    """Yield (obj, per-object catalytic selection) for live objects; prune dead ones.
    No persistent aggregate `catres` selection is built — one giant OR'd selection
    across many objects did not scale to big multi-file/glob loads. Per-object is
    robust and is all the styling needs."""
    live = set(_self.get_object_list())
    for obj in list(_CATRES_BY_OBJ.keys()):
        if obj not in live:
            _CATRES_BY_OBJ.pop(obj, None)
            continue
        yield obj, _catres_obj_sel(obj, _CATRES_BY_OBJ[obj])


def _color_catres(obj, _self=cmd, show_sticks=True):
    cats = _CATRES_BY_OBJ.get(obj)
    if not cats:
        return
    sel = _catres_obj_sel(obj, cats)
    # Highlight only SIDECHAIN carbons (backbone keeps the protein color), EXCEPT
    # PRO/GLY whose backbone carbons ARE highlighted. One merged color op; stick
    # radius comes from the global stick_radius (no per-selection `set`).
    _self.color('catres_hi',
                '((%s) and sidechain and elem C) or ((%s) and resn PRO+GLY and elem C)' % (sel, sel))
    if show_sticks:
        _self.show('sticks', '(%s) and (sidechain or name CA) and not hydro and not resn GLY' % sel)


def _style_one(obj, _self=cmd):
    """Apply the current theme and show ligands, metals and known matcher residues."""
    if obj in getattr(_self, '_rfd4_viewer_state', {}).get('objects', {}):
        _STYLED.add(obj)
        return
    _color_structure(obj, _self)
    _self.show("sticks", "(%s) and organic and not hydro and not (%s)" % (obj, _METAL_SEL))
    show_metals(obj, _self=_self)
    _color_catres(obj, _self=_self)
    _STYLED.add(obj)


def show_catres(selection="all", _self=cmd):
    """Show all known REMARK 666 catalytic residues as highlighted sticks (per object)."""
    n = 0
    for obj, sel in _catres_iter(_self=_self):
        _color_catres(obj, _self=_self)
        n += _self.count_atoms('(%s) and name CA' % sel)
    print('show_catres: %d catalytic residue(s).' % n)


def only_catres(_self=cmd):
    """Hide all protein sidechain sticks EXCEPT the catalytic (REMARK 666) ones.
    Operates per object (no aggregate selection), so it scales to big loads."""
    _self.hide('sticks', 'polymer.protein')
    n = 0
    for obj, sel in _catres_iter(_self=_self):
        stick = '(%s) and (sidechain or name CA) and not hydro and not resn GLY' % sel
        if _self.count_atoms(stick):
            _self.show('sticks', stick)
            _self.set('stick_radius', 0.25, stick)
        _self.color('catres_hi', '(%s) and sidechain and elem C' % sel)
        _self.color('catres_hi', '(%s) and resn PRO+GLY and elem C' % sel)
        n += _self.count_atoms('(%s) and name CA' % sel)
    if n:
        print('only_catres: showing %d catalytic residue(s) as sticks.' % n)
    else:
        print('only_catres: no catalytic (REMARK 666) residues found.')


def _seq_usage():
    print("seq: usage -- seq <selection | object | /path/to/file.pdb>\n"
          "  seq                       the mouse selection 'sele'\n"
          "  seq design                an object\n"
          "  seq chain A and resi 1-50 a sub-selection\n"
          "  seq ~/models/foo.pdb      load a structure file, print its sequence\n"
          "  aliases: sequence, get_seq, get_sequence")


def seq(selection="sele", _self=cmd):
    """Print the continuous single-line one-letter sequence (whole residues).
    `selection` may be the mouse selection `sele` (default), an object/selection
    string, OR a path to a structure file (.pdb/.cif/...), which is loaded
    temporarily just to read its sequence and then removed.
    USAGE: seq  |  seq design  |  seq chain A and resi 1-50  |  seq ~/x.pdb
    Aliases: sequence, get_seq, get_sequence."""
    import os as _o
    _STRUCT_EXT = (".pdb", ".ent", ".cif", ".mmcif", ".mcif", ".pdb.gz",
                   ".cif.gz", ".ent.gz", ".mmtf", ".mol", ".sdf")
    sel = str(selection)
    tmp = None
    path = _o.path.expanduser(sel)
    if _o.path.isfile(path):
        # a real file on disk -> load into a temp object (no view change), read, drop
        tmp = "_seq_file_tmp"
        _self.delete(tmp)
        try:
            _self.load(path, tmp, zoom=0)
        except Exception as e:
            print('seq: could not load file "%s": %s' % (selection, e)); _seq_usage(); return ''
        sel = tmp
    else:
        try:
            n = _self.count_atoms('(%s)' % sel)
        except Exception:
            n = -1
        if n <= 0:
            if sel.lower().endswith(_STRUCT_EXT):
                print('seq: file not found: "%s"' % selection)
            elif n == 0:
                print('seq: selection "%s" is empty (nothing picked?).' % selection)
            else:
                print('seq: bad selection "%s".' % selection)
            _seq_usage(); return ''
    try:
        collected = []
        _self.iterate('byres (%s) and name CA and polymer.protein' % sel,
                      'collected.append((model, segi, chain, resv, resn))',
                      space={'collected': collected})
        if not collected:
            print('seq: no protein residues (CA atoms) in "%s".' % selection); _seq_usage(); return ''
        seen = set(); ordered = []
        for rec in collected:
            k = rec[:4]
            if k in seen:
                continue
            seen.add(k); ordered.append(rec)
        s = ''.join(_AA3TO1.get(r[4], 'X') for r in ordered)
        print(s)
        print('seq: %d residues' % len(s))
        return s
    finally:
        if tmp is not None:
            _self.delete(tmp)


def bg_white(_self=cmd):
    """White background (keeps your render aesthetic; black ray-trace outline)."""
    _self.bg_color('white')
    _self.set('ray_trace_color', 'black')
    print('bg: white')


def bg_black(_self=cmd):
    """Black background (the default)."""
    _self.bg_color('black')
    _self.set('ray_trace_color', 'black')
    print('bg: black')


def autostyle(mode="status", _self=cmd):
    """Toggle auto display of metals (spheres) + catalytic residues (sticks) on load.
    USAGE: autostyle on|off|toggle|status"""
    global _AUTOSTYLE
    t = str(mode).strip().lower()
    if t in ('on', '1', 'true', 'enable'):
        _AUTOSTYLE = True
    elif t in ('off', '0', 'false', 'disable'):
        _AUTOSTYLE = False
    elif t == 'toggle':
        _AUTOSTYLE = not _AUTOSTYLE
    elif t not in ('status', ''):
        print('autostyle: use on|off|toggle|status')
    print('autostyle: %s' % ('on' if _AUTOSTYLE else 'off'))
    return int(_AUTOSTYLE)


def align_all(reference="", cycles=2, _self=cmd):
    """Structurally superpose ALL objects onto a reference (default: first loaded),
    so diverse structures share the same space. Sequence-independent `super` with a
    `cealign` fallback for very divergent pairs. Each object is aligned once to the
    reference (O(N)). Run once after loading a batch/glob.
    USAGE: align_all [reference [, cycles]]"""
    objs = _self.get_object_list()
    if len(objs) < 2:
        print('align_all: need >= 2 objects.')
        return
    ref = reference if (reference and reference in objs) else objs[0]
    try:
        cyc = int(cycles)
    except Exception:
        cyc = 2
    aligned, failed = 0, []
    for obj in objs:
        if obj == ref:
            continue
        ok = False
        try:
            r = _self.super(obj, ref, cycles=cyc)
            ok = bool(r) and r[1] > 0
        except Exception:
            ok = False
        if not ok:
            try:
                _self.cealign(ref, obj)
                ok = True
            except Exception:
                ok = False
        if ok:
            aligned += 1
        else:
            failed.append(obj)
    _self.zoom(ref, buffer=8)
    msg = 'align_all: superposed %d/%d onto %s' % (aligned, len(objs) - 1, ref)
    if failed:
        msg += '  (failed: %s%s)' % (', '.join(failed[:4]), '...' if len(failed) > 4 else '')
    print(msg)


def center_all(reference="", _self=cmd):
    """Instantly co-locate every object's center onto a reference's center
    (translation only, NO rotation). Much faster than align_all; good when you just
    want diverse structures roughly in the same spot. USAGE: center_all [reference]"""
    objs = _self.get_object_list()
    if len(objs) < 2:
        print('center_all: need >= 2 objects.')
        return
    ref = reference if (reference and reference in objs) else objs[0]
    try:
        cx, cy, cz = _self.centerofmass(ref)
    except Exception:
        print('center_all: could not compute reference center.')
        return
    n = 0
    for obj in objs:
        if obj == ref:
            continue
        try:
            ox, oy, oz = _self.centerofmass(obj)
            _self.translate([cx - ox, cy - oy, cz - oz], obj, camera=0)
            n += 1
        except Exception:
            pass
    _self.zoom(ref, buffer=8)
    print('center_all: co-centered %d object(s) on %s' % (n, ref))


def autoalign(mode="status", _self=cmd):
    """Toggle auto-superposition of each newly loaded object onto the first one,
    as it loads (hands-free 'same space' for globs; adds a little load time).
    USAGE: autoalign on|off|toggle"""
    global _AUTOALIGN
    t = str(mode).strip().lower()
    if t in ('on', '1', 'true', 'enable'):
        _AUTOALIGN = True
    elif t in ('off', '0', 'false', 'disable'):
        _AUTOALIGN = False
    elif t == 'toggle':
        _AUTOALIGN = not _AUTOALIGN
    elif t not in ('status', ''):
        print('autoalign: use on|off|toggle|status')
    print('autoalign: %s' % ('on' if _AUTOALIGN else 'off'))
    return int(_AUTOALIGN)


def autosolo(mode="status", _self=cmd):
    """Toggle 'solo load': newly loaded objects beyond the first load HIDDEN, so a big
    glob loads/renders fast; browse with pgup/pgdn (or `structure_show_all`).
    USAGE: autosolo on|off|toggle"""
    global _AUTOSOLO
    t = str(mode).strip().lower()
    if t in ('on', '1', 'true', 'enable'):
        _AUTOSOLO = True
    elif t in ('off', '0', 'false', 'disable'):
        _AUTOSOLO = False
    elif t == 'toggle':
        _AUTOSOLO = not _AUTOSOLO
    elif t not in ('status', ''):
        print('autosolo: use on|off|toggle|status')
    print('autosolo: %s' % ('on' if _AUTOSOLO else 'off'))
    return int(_AUTOSOLO)


def _style_catres_batched(objs, _self=cmd, chunk=40):
    """Style catalytic residues for many objects in O(N/chunk) ops (not per-object): build a
    TRANSIENT combined selection per chunk (object-scoped clauses; residues grouped per chain
    so each clause is short), show+color, then delete it. Transient + object-scoped avoids
    the stale-aggregate-selection ('only the first') and cross-object-chain bugs."""
    items = [(o, _CATRES_BY_OBJ[o]) for o in objs if _CATRES_BY_OBJ.get(o)]
    if not items:
        return
    for i in range(0, len(items), chunk):
        parts = []
        for o, cats in items[i:i + chunk]:
            bychain = {}
            for ch, ri in cats:
                bychain.setdefault(ch, []).append(str(ri))
            for ch, ris in bychain.items():
                parts.append('(%s and chain %s and resi %s)' % (o, ch, '+'.join(ris)))
        if not parts:
            continue
        _self.select('_cat_tmp', ' or '.join(parts))
        _self.show('sticks', '_cat_tmp and (sidechain or name CA) and not hydro and not resn GLY')
        _self.color('catres_hi', '(_cat_tmp and sidechain and elem C) or (_cat_tmp and resn PRO+GLY and elem C)')
    _self.delete('_cat_tmp')
    _self.deselect()


def style_all(selection="all", _self=cmd):
    """Apply the current theme to ordinary structures in one batched operation."""
    objs = _self.get_object_list() if selection in ('all', '(all)', '*') else _self.get_object_list(selection)
    annotated = getattr(_self, '_rfd4_viewer_state', {}).get('objects', {})
    protected = [obj for obj in objs if obj in annotated]
    objs = [obj for obj in objs if obj not in annotated]
    if not objs:
        print('style_all: annotated condition styles preserved.' if protected else 'style_all: no objects.')
        return
    s = selection
    if protected:
        s = '(%s) and not (%s)' % (selection, ' or '.join('model ' + obj for obj in protected))
    previous_updates = _self.get_setting_int('suspend_updates')
    _self.set('suspend_updates', 1)
    try:
        _color_structure(s, _self)
        _self.show('cartoon', '(%s) and polymer' % s)                              # representations
        _self.show('sticks', '(%s) and organic and not (%s) and not hydro' % (s, _METAL_SEL))
        _self.show('spheres', '(%s) and (%s)' % (s, _METAL_SEL))
        _self.hide('everything', '(%s) and hydro' % s)
        _style_catres_batched(objs, _self=_self)                                   # catalytic (chunked, transient)
        _STYLED.update(objs)
    finally:
        _self.set('suspend_updates', previous_updates)
    print('style_all: styled %d object(s).' % len(objs))


_BROWSE_PRESTYLE = [True]    # default ON: browse_list pre-styles ALL up front (one-time load pause, then instant scroll); `browse_prestyle off` for instant start + tiny per-view hitch


def browse_list(listfile, _self=cmd):
    """Load ALL structures in a list (so they all appear in the object list) and show ONLY
    the first. The remaining objects stay hidden. Use `browse_prestyle off` to defer styling
    each ordinary structure until it is first shown.
    USAGE: browse_list <listfile>"""
    global _AUTOSOLO
    try:
        with open(listfile) as fh:
            files = [l.strip() for l in fh if l.strip()]
    except Exception as e:
        print('browse_list: cannot read %s: %s' % (listfile, e))
        return
    if not files:
        print('browse_list: empty list.')
        return
    solo_prev = _AUTOSOLO
    _AUTOSOLO = True             # first -> styled+shown; the rest -> loaded but hidden + styling deferred
    loaded = 0
    try:
        for f in files:
            try:
                _self.load(f)
                loaded += 1
            except Exception as error:
                print("browse_list: could not open %s (%s)" % (f, error))
    finally:
        _AUTOSOLO = solo_prev
    objs = _self.get_object_list()
    first = objs[0] if objs else None
    if _BROWSE_PRESTYLE[0]:
        style_all('all', _self=_self)          # pre-style every object now (one-time pause)
    if first:
        _self.disable('all')
        _self.enable(first)                    # ONLY the first is shown; the rest sit in the sidebar
        _self.zoom(first, buffer=5)
    extra = '  (pre-styled all)' if _BROWSE_PRESTYLE[0] else ''
    print('browse_list: opened %d of %d files%s. PgUp/PgDn to cycle.' % (loaded, len(files), extra))
    return first


def browse_prestyle(mode="status", _self=cmd):
    """Toggle whether browse_list/random_shuffle pre-styles ALL structures up front
    (smooth scrolling, one-time load pause) vs. styling each on first view (instant start,
    tiny per-view pause). USAGE: browse_prestyle on|off|toggle"""
    t = str(mode).strip().lower()
    if t in ('on', '1', 'true', 'enable'):
        _BROWSE_PRESTYLE[0] = True
    elif t in ('off', '0', 'false', 'disable'):
        _BROWSE_PRESTYLE[0] = False
    elif t == 'toggle':
        _BROWSE_PRESTYLE[0] = not _BROWSE_PRESTYLE[0]
    elif t not in ('status', ''):
        print('browse_prestyle: use on|off|toggle|status')
    print('browse_prestyle: %s' % ('on' if _BROWSE_PRESTYLE[0] else 'off'))
    return int(_BROWSE_PRESTYLE[0])


def browse_step(step=1, _self=cmd):
    return structure_step(step, _self=_self)


def browse_next(_self=cmd):
    return structure_step(1, _self=_self)


def browse_prev(_self=cmd):
    return structure_step(-1, _self=_self)


# --- Lean load/fetch hooks: per-object metals + catres only; NO global recenter
#     (this replaces the old center('all')-on-every-load that caused the lag). ---
if not hasattr(cmd, "_pymolrc_originals"):
    def _original_command(name):
        current = getattr(cmd, name)
        if current.__name__ == "_pymolrc_" + name + "_hook":
            current = current.__globals__["_pymolrc_orig_" + name]
        if getattr(current, "_rfd4_dispatch", False):
            current = cmd._rfd4_viewer_state["original_load"]
        return current

    cmd._pymolrc_originals = tuple(_original_command(name) for name in ("load", "fetch", "delete"))
_pymolrc_orig_load, _pymolrc_orig_fetch, _pymolrc_orig_delete = cmd._pymolrc_originals


def _pymolrc_autostyle_new(filename, new_objs, ref=None, _self=cmd):
    cats = _parse_remark666(filename)
    if ref is None:
        all_objs = _self.get_object_list()
        ref = all_objs[0] if all_objs else None
    for obj in new_objs:
        _STYLED.discard(obj)                  # fresh load: clear any stale styled/catres state
        _CATRES_BY_OBJ.pop(obj, None)         # (handles delete + reload of the same name)
        if cats:                              # record catalytic residues (cheap)
            _CATRES_BY_OBJ[obj] = list(cats)
        if _AUTOALIGN and ref and obj != ref:
            try:
                _self.super(obj, ref, cycles=1)
            except Exception:
                pass
        if not _AUTOSTYLE:
            continue
        if obj in getattr(_self, '_rfd4_viewer_state', {}).get('objects', {}):
            # The annotated loader owns its condition colors/representations.
            # Treat it as styled so navigation cannot repaint the conditions.
            _STYLED.add(obj)
            if _AUTOSOLO and ref and obj != ref:
                _self.disable(obj)
            continue
        if _AUTOSOLO and ref and obj != ref:
            _self.disable(obj)                # hidden; styled lazily on first view (pgup/pgdn)
        else:
            _style_one(obj, _self=_self)      # style now (default: every object on load)


# Track only objects managed by these hooks. Before/after names distinguish a
# new load from existing review objects created by load_raw or other plugins.
_KNOWN_OBJS = set()


def _restore_finite_view(previous, new, had_objects):
    """PyMOL can auto-zoom to NaN when an instruction CIF has placeholder coordinates."""
    if all(map(math.isfinite, cmd.get_view())) or not all(map(math.isfinite, previous)):
        return
    cmd.set_view(previous, animate=0)
    if had_objects:
        return
    selections = []
    try:
        for obj in new:
            indices = []
            cmd.iterate_state(1, obj, "if finite(x) and finite(y) and finite(z): keep.append(index)",
                              space={"finite": math.isfinite, "keep": indices})
            if indices:
                name = cmd.get_unused_name("_pymolrc_finite")
                cmd.select_list(name, obj, indices, mode="index")
                selections.append(name)
        if selections:
            cmd.zoom(" or ".join(selections), animate=0)
    finally:
        if selections:
            _pymolrc_orig_delete(" ".join(selections))


def _pymolrc_load_hook(*args, **kwargs):
    before = set(cmd.get_object_list())
    fn = args[0] if args else kwargs.get('filename')
    fmt = args[3] if len(args) > 3 else kwargs.get('format', '')
    suffix = str(fn).lower().removesuffix('.gz')
    cif = str(fmt).lower() in ('cif', 'mmcif') or suffix.endswith(('.cif', '.mmcif'))
    previous_view = cmd.get_view() if cif else None
    result = _pymolrc_orig_load(*args, **kwargs)
    try:
        cur = cmd.get_object_list()
        session = str(fmt).lower() in ('pse', 'psw', 'pze') or suffix.endswith(('.pse', '.psw', '.pze'))
        if session:
            # A session already contains deliberate colors, representations,
            # enabled objects and camera. Do not restyle now or on navigation.
            _KNOWN_OBJS.clear()
            _KNOWN_OBJS.update(cur)
            _STYLED.clear()
            _STYLED.update(cur)
            getattr(cmd, '_rfd4_viewer_state', {}).get('objects', {}).clear()
            return result
        new = [o for o in cur if o not in before]
        if previous_view is not None:
            _restore_finite_view(previous_view, new, before)
        if new:
            _pymolrc_autostyle_new(fn, new, cur[0] if cur else None)   # ref passed (no re-fetch)
            _KNOWN_OBJS.update(new)                        # add-only (O(new)); delete hook prunes
    except Exception as e:
        print('[pymolrc] autostyle (load) skipped: %s' % e)
    return result


def _pymolrc_fetch_hook(*args, **kwargs):
    before = set(cmd.get_object_list())
    result = _pymolrc_orig_fetch(*args, **kwargs)
    try:
        cur = cmd.get_object_list()
        new = [o for o in cur if o not in before]
        if new:
            _pymolrc_autostyle_new(None, new, cur[0] if cur else None)   # fetched: no REMARK 666 MATCH
            _KNOWN_OBJS.update(new)
    except Exception as e:
        print('[pymolrc] autostyle (fetch) skipped: %s' % e)
    return result


def _pymolrc_delete_hook(name='all', *a, **k):
    # Keep _KNOWN_OBJS / per-object state in sync so delete + reload of the same name
    # re-styles correctly. Deletes are rare, so the O(M) reconcile here is fine.
    annotated = getattr(cmd, '_rfd4_viewer_state', {}).get('objects', {})
    previously_live = set(cmd.get_object_list()) if annotated else set()
    result = _pymolrc_orig_delete(name, *a, **k)
    try:
        live = set(cmd.get_object_list())
        _KNOWN_OBJS.intersection_update(live)
        _STYLED.intersection_update(live)
        for stale in annotated.keys() - live:
            record = annotated.pop(stale)
            if stale in previously_live:
                masks = set(cmd.get_names("selections")).intersection(record.get("selection_names", ()))
                if masks:
                    _pymolrc_orig_delete(" ".join(masks))
                group = stale + "_annotations"
                if group in cmd.get_names("group_objects"):
                    # Excise preserves any extra objects the user put in this group.
                    cmd.group(group, action="excise")
        for stale in _CATRES_BY_OBJ.keys() - live:
            del _CATRES_BY_OBJ[stale]
    except Exception:
        pass
    return result


cmd.delete = _pymolrc_delete_hook


cmd.load = _pymolrc_load_hook
cmd.fetch = _pymolrc_fetch_hook
cmd.keyword["load"][0] = cmd.load
cmd.keyword["delete"][0] = cmd.delete
cmd.keyword["fetch"][0] = cmd.fetch


cmd.extend('structure_step', structure_step)
cmd.extend('structure_next', structure_next)
cmd.extend('structure_prev', structure_prev)
cmd.extend('structure_show_all', structure_show_all)
cmd.extend('bind_structure_navigation_keys', bind_structure_navigation_keys)
cmd.extend('show_metals', show_metals)
cmd.extend('show_catres', show_catres)
cmd.extend('only_catres', only_catres)
cmd.extend('seq', seq)
cmd.extend('sequence', seq)
cmd.extend('get_seq', seq)
cmd.extend('get_sequence', seq)
cmd.extend('bg_white', bg_white)
cmd.extend('bg_black', bg_black)
cmd.extend('autostyle', autostyle)
cmd.extend('align_all', align_all)
cmd.extend('center_all', center_all)
cmd.extend('autoalign', autoalign)
cmd.extend('autosolo', autosolo)
cmd.extend('style_all', style_all)
cmd.extend('browse_list', browse_list)
cmd.extend('browse_next', browse_next)
cmd.extend('browse_prev', browse_prev)
cmd.extend('browse_step', browse_step)
cmd.extend('browse_prestyle', browse_prestyle)

bind_structure_navigation_keys(_self=cmd)
