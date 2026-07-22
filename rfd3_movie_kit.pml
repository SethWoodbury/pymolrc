# rfd3_movie_kit.pml -- standalone RFdiffusion3 trajectory movie/figure toolkit.
#
# A single, self-contained file with JUST the RFd3 commands from the pymolrc
# repo (https://github.com/SethWoodbury/pymolrc), decoupled from that repo's
# foliage color scheme, style_all, and catalytic-REMARK666 machinery. Doesn't
# touch or require your existing .pymolrc at all -- just adds these commands
# on top of it.
#
# INSTALL (pick one):
#   1) Download just this file and load it from your OWN .pymolrc:
#        curl -O https://raw.githubusercontent.com/SethWoodbury/pymolrc/main/rfd3_movie_kit.pml
#      then add one line to the end of your ~/.pymolrc:
#        run /path/to/rfd3_movie_kit.pml
#   2) Or just `run rfd3_movie_kit.pml` once per PyMOL session.
#
# See docs/RFD3_QUICKSTART.md in the repo for the short version of this, and
# docs/RFD3_FIGURES.md for the full "why" behind every non-obvious line here.
#
# Commands this adds: rfd3_movie, rfd3_check_frames, style_fixed, apply_camera,
# draw_connectors, catalytic_sel_from_motif, add_custom_bond, remove_bond,
# color_bb_rfdiffusion3.

set_color RFd_darkblue, [75,95,170]
set_color RFd_purple, [213,154,181]
set_color RFd_pink, [255,172,183]
set_color paper_teal, [79,185,175]

python


def _rfd_bool(v):
    return str(v).strip().lower() in ('1', 'true', 'on', 'yes', 't', 'y')


def color_bb_rfdiffusion3(selection="chain A", all_atom=0, backbone_only=0, _self=cmd):
    """Color `selection` N->C with the RFdiffusion3 gradient
    (pink -> purple -> teal -> dark blue; teal auto-centered at the chain midpoint).

      selection       one or more chains / any selection (default 'chain A')
      all_atom=1      recolor every element, not just carbons
      backbone_only=1 restrict to backbone atoms (sidechains keep their colors)
    USAGE: color_bb_rfdiffusion3 [sel] [, all_atom] [, backbone_only]"""
    _rv = []
    _self.iterate("(%s) and name CA and polymer" % selection,
                  "_rv.append(resv)", space={'_rv': _rv})
    if not _rv:
        return
    lo, hi = min(_rv), max(_rv)
    mid = (lo + hi) // 2
    sel = "(%s)" % selection
    if not _rfd_bool(all_atom):
        sel += " and elem C"
    if _rfd_bool(backbone_only):
        sel += " and backbone"
    _self.spectrum("resi", "RFd_pink RFd_purple paper_teal",
                   "%s and resi %d-%d" % (sel, lo, mid))
    _self.spectrum("resi", "paper_teal RFd_darkblue",
                   "%s and resi %d-%d" % (sel, mid + 1, hi))


def style_fixed(obj, trim_atoms=None, glu_drop="N+CA+C+O+CB", his_drop="N+CA+C+O",
                 fixed_color="orange", color_map=None, lig_rep="sticks",
                 stick_r=0.16, bb_stick_r=0.30, real_sphere=0.28, vplace_sphere=0.18,
                 _self=cmd):
    """Solid ball-and-stick on a trajectory's fixed theozyme (motif residues + ligand
    + metals; everything not resn UNK). Fixes the "fixed atoms never look opaque"
    problem AT ITS ROOT: small disconnected motif fragments (no backbone/no bonds)
    fall back to the `nonbonded`/`lines` reps, which ignore sphere_transparency and
    stick_transparency entirely -- no amount of `set ..._transparency, 0` fixes that
    on its own. Fix: hide nonbonded+lines, rebuild sticks/spheres fresh, THEN force
    transparency 0 last. Idempotent -- safe to call again any time you want to
    restyle (different color_map/lig_rep/etc) without reloading the trajectory.
    Works on the live trajectory object AND on a single-state cmd.create() copy.

      trim_atoms   dict of {resn: atom_names} to trim any FULLY fixed residue (every
                   sidechain atom fixed, not just CB) down to its reactive atoms only
                   -- e.g. a fully-fixed Glu/His where only the carboxylate/imidazole
                   should show, not the placeholder backbone too. Default:
                   {"GLU": glu_drop, "HIS": his_drop} (the two seen so far); add your
                   own {resn: "N+CA+C+O..."} entries for any other residue that needs
                   it, or pass {} to disable trimming entirely. RFdiffusion3
                   sidechains use generic V0-V8 placeholder names, so real backbone
                   atom names are the only thing safe to drop by name.
      glu_drop / his_drop  shorthand for the default GLU/HIS entries above; ignored
                   if `trim_atoms` is given explicitly.
      fixed_color  color for fixed-protein carbons not overridden by `color_map`.
      color_map    optional dict of {key: color} applied AFTER fixed_color, for
                   per-residue/selection overrides, e.g. {"HIS": "skyblue",
                   "chain B and resi 45": "magenta"}. A key is used as a selection
                   as-is if it looks like one (contains "and"/"or" or a keyword like
                   chain/resi/resn/name/index), else treated as a bare `resn <key>`.
      lig_rep      representation for the ligand -- 'sticks' (default), or 'lines'/
                   'nb_spheres'/etc for a lighter look.
    USAGE: style_fixed <obj>"""
    if trim_atoms is None:
        trim_atoms = {"GLU": glu_drop, "HIS": his_drop}
    if color_map is None:
        color_map = {}

    fx = "(%s) and not resn UNK" % obj
    lig = "(%s) and hetatm and not metals" % fx
    met = "(%s) and metals" % fx
    prot = "(%s) and not hetatm" % fx

    _self.hide("nonbonded", fx)
    _self.hide("lines", fx)
    _self.hide("spheres", fx)
    _self.show("sticks", fx)
    _self.show("spheres", prot)
    _self.show("spheres", met)

    for resn, drop in trim_atoms.items():
        extra = "(%s) and resn %s and name %s" % (fx, resn, drop)
        _self.hide("sticks", extra)
        _self.hide("spheres", extra)

    if lig_rep != "sticks":
        _self.hide("sticks", lig)
    _self.hide("lines", lig)
    _self.hide("nonbonded", lig)
    _self.show(lig_rep, lig)
    if lig_rep == "sticks":
        _self.set("stick_radius", 0.2, lig)

    _self.set("stick_radius", stick_r, prot)
    _self.set("stick_radius", bb_stick_r, "(%s) and name N+CA+C+O" % prot)
    _self.set("sphere_scale", real_sphere, "(%s) and name CA+CB" % prot)
    _self.set("sphere_scale", vplace_sphere, "(%s) and not name N+CA+C+O+CB" % prot)

    _self.color(fixed_color, "(%s) and elem C" % prot)
    _kw = ("chain ", "resi ", "resn ", "name ", "index ", " and ", " or ")
    for key, col in color_map.items():
        sel = key if any(tok in key for tok in _kw) else ("resn %s" % key)
        _self.color(col, "(%s) and (%s) and elem C" % (prot, sel))

    _self.set("stick_transparency", 0, fx)
    _self.set("sphere_transparency", 0, fx)
    print("style_fixed: %d fixed atoms styled opaque on %s." % (_self.count_atoms(fx), obj))


def rfd3_check_frames(obj, reverse=None, n_states=None, _self=cmd):
    """Sanity-check which end of a loaded RFd3 trajectory is noise vs folded, and
    (if `reverse` is given) whether that setting actually plays noise -> folded.

    rfd3_movie() does NOT distinguish a `_noisy_model_` trajectory file from a
    `_denoised_model_` one -- both load and render identically, and the filename
    is only ever used to auto-derive the design JSON's path. This checks the
    actual GEOMETRY instead of trusting the filename or a hardcoded assumption:
    RFd3 noise is a collapsed blob near the origin (small radius), the folded
    design is expanded (large radius) -- see docs/RFD3_FIGURES.md item 3. Measures
    the CA radius-from-centroid of the diffusing chain (resn UNK) at state 1 vs
    the last state, using literal state numbers (unaffected by any `reverse`
    frame remapping already applied via mset).

      reverse   pass the same value you used (or plan to use) for rfd3_movie's
                `reverse=` -- prints whether it matches what this trajectory
                actually needs, or is backwards.
    Returns the reverse value (0 or 1) this trajectory needs for noise -> folded
    playback, or None if inconclusive. USAGE: rfd3_check_frames <obj> [, reverse]
    """
    import math
    ns = n_states or _self.count_states(obj)

    def _radius(state):
        m = _self.get_model("(%s) and resn UNK and name CA" % obj, state=state)
        coords = [a.coord for a in m.atom]
        if len(coords) < 3:
            return None
        n = len(coords)
        cx = sum(c[0] for c in coords) / n
        cy = sum(c[1] for c in coords) / n
        cz = sum(c[2] for c in coords) / n
        return max(math.sqrt((c[0] - cx) ** 2 + (c[1] - cy) ** 2 + (c[2] - cz) ** 2) for c in coords)

    r1, rN = _radius(1), _radius(ns)
    if r1 is None or rN is None:
        print("rfd3_check_frames: no resn UNK CA atoms found on %s (wrong object, or not an RFd3 trajectory?)." % obj)
        return None
    print("rfd3_check_frames: state 1 radius %.1f A, state %d radius %.1f A." % (r1, ns, rN))
    if abs(r1 - rN) < 2.0:
        print("  Too close to call from the endpoints alone -- try comparing against a middle state, "
              "or trust the design-JSON metadata instead.")
        return None
    if r1 < rN:
        need_reverse = 0
        print("  state 1 = NOISE (smaller/collapsed), state %d = FOLDED (larger/expanded)." % ns)
    else:
        need_reverse = 1
        print("  state 1 = FOLDED (larger/expanded), state %d = NOISE (smaller/collapsed)." % ns)
    print("  For noise -> folded playback: reverse=%d." % need_reverse)
    if reverse is not None:
        got = 1 if _rfd_bool(reverse) else 0
        if got == need_reverse:
            print("  reverse=%d matches -- frame 1 is noise, frame %d is folded. Good." % (got, ns))
        else:
            print("  MISMATCH: reverse=%d was used, but this trajectory needs reverse=%d for noise -> folded. "
                  "Frames are currently playing backwards (folded -> noise)." % (got, need_reverse))
    return need_reverse


def rfd3_movie(traj_file, color_scheme=1, cloud=1, fixed_color="orange",
               fixed_sidechain=0, reverse=0, fixed_json="", name="rfd3_traj",
               ca_scale=0.4, cloud_scale=0.13, cloud_transparency=0.55,
               bonds=None, unbonds=None, _self=cmd, **style_kwargs):
    """Set up an RFdiffusion3 diffusion movie from ONE trajectory .cif.gz.

    A `*_noisy_model_*.cif.gz` / `*_denoised_model_*.cif.gz` file is itself a
    100-state diffusion trajectory. Loads that ONE file. DIFFUSING residues (resn
    UNK) show as CA spheres + a faint V0-V8 sidechain cloud, colored by the
    RFdiffusion3 gradient by default. FIXED atoms (the catalytic motif + cofactor +
    metal) are auto-detected, held static, and styled via style_fixed() (see its
    docstring for the full customization surface -- trim_atoms/color_map/lig_rep/
    etc, all forwardable here as extra keyword args). Then type `mplay`.

    On real trajectories seen so far, raw state 1 is the folded design and state N
    is noise -- backwards from what's intuitive -- so `reverse=1` is usually what
    you want for noise -> folded playback. This isn't hardcoded: rfd3_movie doesn't
    treat a `_noisy_model_` file any differently from a `_denoised_model_` one, so
    every call auto-runs rfd3_check_frames() below and prints whether your
    `reverse=` setting actually plays noise -> folded for THIS file. See
    docs/RFD3_FIGURES.md.

    If the design JSON is found -- auto-derived from the trajectory name
    (`..._model_N.json`) or given as `fixed_json=<path>` -- its `select_fixed_atoms`
    residue count is cross-checked against detection and any mismatch is warned.

      color_scheme    : 1/'rfd3' (default) = RFd3 gradient on the diffusing protein;
                        0/'cyan' = flat cyan; or supply ANY PyMOL color name for a
                        flat custom color (e.g. 'hotpink'), or a space-separated list
                        of color names (e.g. 'marine yellow red') for a custom N-stop
                        gradient N->C via cmd.spectrum -- same mechanism
                        color_bb_rfdiffusion3 uses internally, just with your colors.
      cloud           : 1 = faint V-atom sidechain cloud on; 0 = CA spheres only.
      ca_scale        : sphere_scale of the diffusing chain's CA trace (default 0.4).
      cloud_scale / cloud_transparency : size/transparency of the V-atom cloud
                        (defaults 0.13 / 0.55).
      fixed_color     : color for the fixed-protein carbons (default 'orange'); for
                        finer per-residue control use `color_map` (see style_fixed).
      fixed_sidechain : 1 = show only the sidechains of fixed residues (hide backbone).
      reverse         : 1 = play the trajectory in reverse (states N..1). Verified
                        automatically every call -- see above.
      fixed_json      : optional path to the design JSON (else auto-derived).
      bonds           : optional list of (sel_a, sel_b[, label]) tuples, forced
                        immediately via add_custom_bond() -- e.g. a covalent
                        TS-adduct PyMOL won't perceive on its own.
      unbonds         : optional list of (sel_a, sel_b[, label]) tuples, removed
                        immediately via remove_bond() -- e.g. a spurious bond PyMOL
                        wrongly perceived.
      **style_kwargs  : any other keyword (trim_atoms=, color_map=, lig_rep=,
                        stick_r=, ...) is passed straight through to style_fixed().
                        Call style_fixed() again by hand any time afterward to
                        restyle the fixed theozyme without reloading the trajectory.
    Fixed-atom detection uses `not resn UNK`. Uses `defer_builds_mode 3`; never
    toggle representations per-frame on this object once loaded -- pull a single
    state out with cmd.create() first (see docs/RFD3_FIGURES.md #4). `set
    defer_builds_mode, 0` to return to fast single-structure browsing.
    """
    import os as _o
    import json as _json
    _self.set("defer_builds_mode", 3)
    _self.delete(name)
    _self.load(_o.path.expanduser(traj_file), name)
    ns = _self.count_states(name)
    move = "(%s) and resn UNK" % name
    fixed = "(%s) and not resn UNK" % name
    Vsel = "%s and name V0+V1+V2+V3+V4+V5+V6+V7+V8" % move
    _self.hide("everything", name)
    _self.show("spheres", "%s and name CA" % move)
    _self.set("sphere_scale", ca_scale, "%s and name CA" % move)
    if _rfd_bool(cloud):
        _self.show("spheres", Vsel)
        _self.set("sphere_scale", cloud_scale, Vsel)
        _self.set("sphere_transparency", cloud_transparency, Vsel)
    _scheme_tokens = ("1", "true", "on", "yes", "t", "y",
                      "0", "false", "off", "no", "f", "n", "", "rfd3", "cyan")
    _custom_scheme = (color_scheme.strip()
                       if isinstance(color_scheme, str) and color_scheme.strip().lower() not in _scheme_tokens
                       else None)
    if _custom_scheme:
        if " " in _custom_scheme:
            _self.spectrum("resi", _custom_scheme, "%s and elem C" % move)
            if _rfd_bool(cloud):
                _self.spectrum("resi", _custom_scheme, Vsel)
        else:
            _self.color(_custom_scheme, "%s and elem C" % move)
            if _rfd_bool(cloud):
                _self.color(_custom_scheme, Vsel)
    elif _rfd_bool(color_scheme) or (isinstance(color_scheme, str) and color_scheme.strip().lower() == "rfd3"):
        color_bb_rfdiffusion3(move, all_atom=1, _self=_self)
    else:
        _self.color("cyan", "%s and name CA" % move)
        if _rfd_bool(cloud):
            _self.color("grey70", Vsel)
    style_fixed(name, fixed_color=fixed_color, _self=_self, **style_kwargs)
    if _rfd_bool(fixed_sidechain):
        _prot = "(%s) and not hetatm" % fixed
        _self.hide("everything", "%s and backbone" % _prot)
    for _item in (bonds or []):
        add_custom_bond(_item[0], _item[1], label=(_item[2] if len(_item) > 2 else ""), _self=_self)
    for _item in (unbonds or []):
        remove_bond(_item[0], _item[1], label=(_item[2] if len(_item) > 2 else ""), _self=_self)
    _jp = _o.path.expanduser(fixed_json) if fixed_json else ""
    if not _jp:
        _b = _o.path.expanduser(str(traj_file))
        for _suf in (".cif.gz", ".cif", ".pdb.gz", ".pdb"):
            if _b.endswith(_suf):
                _b = _b[:-len(_suf)]
                break
        _jp = _b.replace("_noisy_model_", "_model_").replace("_denoised_model_", "_model_") + ".json"
    _note = ""
    if _o.path.exists(_jp):
        try:
            _spec = _json.load(open(_jp)).get("specification", {})
            _nj = len(_spec.get("select_fixed_atoms", {}))
            _rset = set()
            _self.iterate("(%s) and not hetatm" % fixed, "_rset.add((chain, resi))",
                          space={'_rset': _rset})
            _nd = len(_rset)
            if _nj == _nd:
                _note = " JSON select_fixed_atoms: %d residues (matches)." % _nj
            else:
                _note = " WARNING: JSON lists %d fixed residues, detected %d." % (_nj, _nd)
        except Exception as _e:
            _note = " (could not read JSON: %s)" % _e
    _self.mset(("%d -1" % ns) if _rfd_bool(reverse) else ("1 -%d" % ns))
    _self.set("all_states", 0)
    _self.orient(name)
    _self.frame(1)
    print("rfd3_movie: %d states; %d fixed atoms (%s).%s Type `mplay` (noise -> folded)."
          % (ns, _self.count_atoms(fixed), fixed_color, _note))
    # Auto direction-check: doesn't matter whether traj_file is a _noisy_model_
    # or _denoised_model_ trajectory -- neither is treated specially above, so
    # verify the actual geometry instead of trusting the filename or a hardcoded
    # assumption. See rfd3_check_frames() above and docs/RFD3_FIGURES.md item 3.
    rfd3_check_frames(name, reverse=reverse, n_states=ns, _self=_self)


def add_custom_bond(sel_a, sel_b, label="", _self=cmd):
    """Force a bond PyMOL won't perceive on its own (e.g. across a protein/ligand
    boundary, or between two disconnected motif fragments) -- e.g. a covalent
    TS-adduct. cmd.bond writes topology once, so it holds across every state of a
    multi-state object; don't call this inside a per-frame loop. Requires each
    selection to resolve to exactly one atom. USAGE: add_custom_bond <sel_a>, <sel_b>
    [, label]"""
    na, nb = _self.count_atoms(sel_a), _self.count_atoms(sel_b)
    if na != 1 or nb != 1:
        print("add_custom_bond [%s]: SKIPPED, matched %d,%d atoms (need exactly 1,1)."
              % (label, na, nb))
        return False
    _self.bond(sel_a, sel_b)
    _self.show("sticks", "(%s) or (%s)" % (sel_a, sel_b))
    print("add_custom_bond [%s]: %.2f A." % (label, _self.get_distance(sel_a, sel_b)))
    return True


def remove_bond(sel_a, sel_b, label="", _self=cmd):
    """Remove a bond PyMOL perceived that you don't want -- the mirror of
    add_custom_bond() (e.g. clearing a spurious bond PyMOL guessed across a chain
    break, clash, or metal coordination it shouldn't have drawn). Requires each
    selection to resolve to exactly one atom. USAGE: remove_bond <sel_a>, <sel_b>
    [, label]"""
    na, nb = _self.count_atoms(sel_a), _self.count_atoms(sel_b)
    if na != 1 or nb != 1:
        print("remove_bond [%s]: SKIPPED, matched %d,%d atoms (need exactly 1,1)."
              % (label, na, nb))
        return False
    _self.unbond(sel_a, sel_b)
    print("remove_bond [%s]: bond removed." % label)
    return True


def catalytic_sel_from_motif(clean_obj, motif_sel, max_match=4.0, _self=cmd):
    """Identify a clean (real-sequence) output model's catalytic residues by
    coordinate-matching them to a trajectory's FIXED motif selection -- the theozyme
    is fixed, so it occupies near-identical coordinates in both the raw noisy
    trajectory and the final clean model. Far more precise than a distance cutoff
    off the ligand (e.g. "within 4.5 of resn LIG"), which grabs the whole first shell
    and over-colors. Returns a selection string, or '' if nothing matched within
    max_match Angstroms. USAGE: catalytic_sel_from_motif <clean_obj>, <motif_sel>
    [, max_match]"""
    import collections as _collections

    def _anchor(atoms):
        for want in ("CB", "CA"):
            for a in atoms:
                if a.name == want:
                    return a.coord
        n = len(atoms)
        return [sum(a.coord[i] for a in atoms) / n for i in range(3)]

    def _dist(a, b):
        return sum((a[i] - b[i]) ** 2 for i in range(3)) ** 0.5

    mm = _self.get_model(motif_sel)
    if not mm.atom:
        print("catalytic_sel_from_motif: motif_sel matched no atoms.")
        return ""
    motif = _collections.OrderedDict()
    for a in mm.atom:
        motif.setdefault((a.chain, a.resi), []).append(a)

    cm = _self.get_model("(%s) and polymer" % clean_obj)
    cres = _collections.OrderedDict()
    for a in cm.atom:
        cres.setdefault((a.chain, a.resi), []).append(a)
    cres_anchors = {k: _anchor(v) for k, v in cres.items()}

    matched = set()
    for _, atoms in motif.items():
        r = _anchor(atoms)
        best = min(((_dist(r, anc), k) for k, anc in cres_anchors.items()), default=None)
        if best and best[0] <= max_match:
            matched.add(best[1])
    if not matched:
        print("catalytic_sel_from_motif: no residue matched within %.1f A." % max_match)
        return ""

    by_chain = _collections.defaultdict(list)
    for ch, resi in matched:
        by_chain[ch].append(resi)

    def _resi_num(r):
        digits = "".join(c for c in r if c.isdigit())
        return int(digits) if digits else 0

    parts = ["(chain %s and resi %s)" % (ch, "+".join(sorted(resis, key=_resi_num)))
             for ch, resis in by_chain.items()]
    return "(%s) and (%s)" % (clean_obj, " or ".join(parts))


def apply_camera(obj, mode="orient", sel=None, custom_view=None, zoom_buffer=5.0,
                  zoom_to=None, span_states=False, n_states=None, turns=(), _self=cmd):
    """One entry point for camera control instead of ad hoc orient/zoom calls.
      mode="orient" -- auto-frame (whole-protein fold), the default.
      mode="zoom"   -- tight crop on `sel` (e.g. the theozyme) for active-site figures.
      mode="view"   -- exact `get_view` matrix passed in `custom_view` (18 floats;
                       run `print(cmd.get_view())` after manually posing a view once).
      zoom_to       -- after framing/rotating, re-center/crop on this selection.
      span_states + n_states -- for ONE fixed camera meant to hold across a whole
                       multi-state trajectory: takes rotation+zoom from the FOLDED
                       frame. RFdiffusion(3) noise is a collapsed blob near the
                       origin (small radius); the chain EXPANDS as it folds, so the
                       folded frame is the largest state and the one that can clip --
                       frame the camera to it, not to the noise.
      turns         -- e.g. [("y", 180)], applied in order after framing/rotation.
    USAGE: apply_camera <obj> [, mode] [, sel] [, zoom_buffer] [, turns]"""
    if mode == "view" and custom_view is not None:
        _self.set_view(custom_view)
        if zoom_to:
            _self.zoom(zoom_to, zoom_buffer)
    elif mode == "zoom":
        _self.zoom(sel or obj, zoom_buffer)
    else:
        target = sel or obj
        if span_states and n_states:
            _self.frame(n_states)
            _self.orient(target)
            _self.zoom(target, zoom_buffer)
        else:
            _self.orient(target)
        if zoom_to:
            _self.zoom(zoom_to, zoom_buffer)
    for axis, deg in turns:
        if deg:
            _self.turn(axis, deg)


def draw_connectors(obj, color="grey50", max_dist=14.0, _self=cmd):
    """Dashed lines from each fixed motif residue (in a trajectory's `not resn UNK`
    selection) to its nearest diffusing-chain CA -- a visual stand-in for
    "this theozyme residue sits here in the scaffold" when there's no covalent bond
    to draw (most RFdiffusion3 motif residues are disconnected fragments). Only
    meaningful on a FOLDED frame -- nearest-CA is close to random during noise.
    USAGE: draw_connectors <obj> [, color] [, max_dist]"""
    ca_atoms = _self.get_model("(%s) and resn UNK and name CA" % obj).atom
    ca = [(a.index, a.coord) for a in ca_atoms]
    if not ca:
        print("draw_connectors: no diffusing-chain CA atoms found on %s." % obj)
        return 0

    resis = []
    _self.iterate("(%s) and not resn UNK and not hetatm" % obj,
                  "resis.append((chain, resi))", space={"resis": resis})

    seen, n = set(), 0
    for chain, resi in resis:
        if (chain, resi) in seen:
            continue
        seen.add((chain, resi))
        anchor_atom = None
        for pick in ("name CB", "name CA", "name V0", "all"):
            asel = "(%s) and chain %s and resi %s and %s" % (obj, chain, resi, pick)
            if _self.count_atoms(asel):
                anchor_atom = _self.get_model(asel).atom[0]
                break
        if anchor_atom is None:
            continue
        ax, ay, az = anchor_atom.coord
        best = min(ca, key=lambda c: (ax - c[1][0]) ** 2 + (ay - c[1][1]) ** 2 + (az - c[1][2]) ** 2)
        if ((ax - best[1][0]) ** 2 + (ay - best[1][1]) ** 2 + (az - best[1][2]) ** 2) ** 0.5 > max_dist:
            continue
        dname = "conn_%d" % n
        n += 1
        _self.distance(dname, "(%s) and index %d" % (obj, anchor_atom.index),
                        "(%s) and index %d" % (obj, best[0]))
        _self.hide("labels", dname)
        for k, v in (("dash_color", color), ("dash_gap", 0.35), ("dash_length", 0.4), ("dash_radius", 0.05)):
            _self.set(k, v, dname)
    print("draw_connectors: drew %d connector(s) on %s." % (n, obj))
    return n


def help_rfd3(_self=cmd):
    """Command reference for this RFdiffusion3 movie/figure kit."""
    print("""
=== RFdiffusion3 movie / figure kit ===  (full write-up: docs/RFD3_FIGURES.md)

rfd3_movie <traj.cif.gz> [, reverse=0] [, color_scheme=1] [, cloud=1]
           [, fixed_sidechain=0] [, fixed_json=] [, ca_scale=0.4]
           [, cloud_scale=0.13] [, cloud_transparency=0.55] [, bonds=] [, unbonds=]
    MANDATORY: traj_file. Load ONE diffusion trajectory. Diffusing chain = CA
    spheres + faint sidechain cloud; fixed motif/cofactor/metal auto-styled.
    reverse=1 RECOMMENDED (frame 1 = noise, last = folded). color_scheme:
    1/'rfd3' gradient, 0/'cyan' flat, any color name = flat, "c1 c2 c3" = custom
    gradient. Extra kwargs (trim_atoms, color_map, lig_rep, ...) pass to style_fixed.
rfd3_check_frames <obj> [, reverse] [, n_states]
    Measure which end is noise vs folded from geometry (runs auto after rfd3_movie).
style_fixed <obj> [, trim_atoms=] [, fixed_color=orange] [, color_map=] [, lig_rep=sticks]
    MANDATORY: obj. Solid opaque ball-and-stick on the fixed theozyme. Idempotent.
apply_camera <obj> [, mode=orient] [, sel=] [, custom_view=] [, span_states=] [, n_states=] [, turns=]
    Camera control. span_states=1,n_states=N holds one folded-frame camera across a traj.
draw_connectors <obj> [, color=grey50] [, max_dist=14.0]
catalytic_sel_from_motif <clean_obj>, <motif_sel> [, max_match=4.0]
add_custom_bond <sel_a>, <sel_b> [, label]   |   remove_bond <sel_a>, <sel_b> [, label]
color_bb_rfdiffusion3 [sel]     RFd3 N->C gradient (pink->purple->teal->darkblue).

Quickstart: RFD3_QUICKSTART.md   Examples: examples/
""")


cmd.extend('color_bb_rfdiffusion3', color_bb_rfdiffusion3)
cmd.extend('style_fixed', style_fixed)
cmd.extend('rfd3_movie', rfd3_movie)
cmd.extend('rfd3_check_frames', rfd3_check_frames)
cmd.extend('add_custom_bond', add_custom_bond)
cmd.extend('remove_bond', remove_bond)
cmd.extend('catalytic_sel_from_motif', catalytic_sel_from_motif)
cmd.extend('apply_camera', apply_camera)
cmd.extend('draw_connectors', draw_connectors)
cmd.extend('help_rfd3', help_rfd3)
cmd.extend('help_rfd3_movie', help_rfd3)

print("rfd3_movie_kit: loaded 10 commands (rfd3_movie, rfd3_check_frames, style_fixed, "
      "apply_camera, draw_connectors, catalytic_sel_from_motif, add_custom_bond, "
      "remove_bond, color_bb_rfdiffusion3, help_rfd3).")

python end
