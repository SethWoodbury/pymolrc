#!/usr/bin/env python3
"""
Worked example: a Figure-1-style RFdiffusion3 trajectory figure, using the
consolidated .pymolrc commands (rfd3_movie, style_fixed, apply_camera,
draw_connectors, catalytic_sel_from_motif, add_custom_bond) instead of ad hoc
per-project styling code.

Run from Terminal (uses your own PyMOL install, not a sandbox):
    pymol -cq examples/make_rfd3_figure.py

Or from inside an interactive PyMOL session:
    run examples/make_rfd3_figure.py

Uses the trajectory bundled in example_pdbs/. See docs/RFD3_FIGURES.md for
the reasoning behind every non-obvious line here (reverse=1, the opacity
fix, the noise-is-collapsed camera trick, why cartoon toggles per-frame hang
PyMOL, etc.) -- this script deliberately doesn't repeat those explanations
inline.
"""
import os
from pymol import cmd

HERE = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else os.getcwd()
REPO = os.path.dirname(HERE)

TRAJ = os.path.join(REPO, "example_pdbs", "ZAPP_p1D1_i14_rfd3_noisy_trajectory.cif.gz")
OUTDIR = os.path.join(HERE, "output")
os.makedirs(OUTDIR, exist_ok=True)

NAME = "rfd3_traj"

# rfd3_movie's own docstring says reverse=0 gives noise->folded, but that's
# backwards on real trajectories -- state 1 is the true folded structure,
# state N is true noise. reverse=1 flips PyMOL's frame->state mapping so
# frame 1 = noise, frame n_states = folded, matching everything below.
# (docs/RFD3_FIGURES.md #1)
rfd3_movie(TRAJ, color_scheme=1, cloud=1, fixed_color="orange", name=NAME, reverse=1)

n_states = cmd.count_states(NAME)
fixed = "(%s) and not resn UNK" % NAME
move = "(%s) and resn UNK" % NAME

# Replaces this project's old hand-rolled opacity/trim/ligand-styling block --
# style_fixed does all of it in one call, and fixes the opacity bug at its
# root (docs/RFD3_FIGURES.md #2) instead of papering over it.
style_fixed(NAME)

# A little transparency on the diffusing CA trace so the theozyme reads through.
cmd.set("sphere_transparency", 0.35, "%s and name CA" % move)

# Dashed lines from each fixed motif residue to its nearest scaffold CA --
# only meaningful once the chain has actually folded (docs/RFD3_FIGURES.md #10
# note on draw_connectors), so this gets called again at the folded frame below.

cmd.bg_color("white")
cmd.set("ray_opaque_background", 1)
cmd.set("antialias", 2)
cmd.set("ray_trace_mode", 1)
cmd.set("ray_shadows", 0)

# ONE fixed camera held across the whole strip, taken from the FOLDED frame
# (the largest state -- RFd3 noise is a collapsed blob, not an expansion;
# docs/RFD3_FIGURES.md #3). Without this, panels drawn while the chain is
# still noise/small would zoom in tight, then "pull back" as it unfolds.
apply_camera(NAME, mode="orient", span_states=True, n_states=n_states, zoom_buffer=6.0)

# Representative frames spanning noise -> folded (Fig 1-style panel strip).
fractions = [0.0, 0.25, 0.5, 0.75, 1.0]
frames = sorted({max(1, round(f * (n_states - 1)) + 1) for f in fractions})

for i, fr in enumerate(frames):
    cmd.frame(fr)
    cmd.ray(2400, 1800)
    outpath = os.path.join(OUTDIR, "rfd3_snapshot_%02d_state%03d.png" % (i + 1, fr))
    cmd.png(outpath, dpi=600)
    print("saved %s" % outpath)

# --- Bonus panel: final (folded) frame, cartoon backbone + connector lines ---
# Built on a separate single-state COPY (cmd.create), never on the live
# multi-state object -- toggling cartoon on the deferred trajectory object
# hangs PyMOL (docs/RFD3_FIGURES.md #4). Also note the state index here is
# n_states, the REAL final state -- NOT reinterpreted through the reverse
# mapping, since cmd.create takes a literal state number (docs/RFD3_FIGURES.md #1).
cmd.frame(n_states)
FINAL = "rfd3_final"
cmd.create(FINAL, NAME, n_states, 1)
final_move = "(%s) and resn UNK" % FINAL
cmd.hide("everything", final_move)
cmd.set("cartoon_trace_atoms", 1, final_move)
cmd.set("ribbon_trace_atoms", 1, final_move)
cmd.show("cartoon", final_move)
color_bb_rfdiffusion3(final_move, all_atom=1, _self=cmd)
apply_camera(FINAL, mode="orient", zoom_buffer=6.0)
draw_connectors(FINAL, color="grey50", max_dist=14.0)
cmd.ray(2400, 1800)
final_path = os.path.join(OUTDIR, "rfd3_final_frame_cartoon.png")
cmd.png(final_path, dpi=600)
cmd.delete(FINAL)
print("saved %s" % final_path)

# --- Bonus panel: the theozyme alone (fixed catalytic atoms only) ---
cmd.hide("everything", move)
apply_camera(NAME, mode="zoom", sel=fixed, zoom_buffer=3.0)
cmd.ray(2400, 1800)
theozyme_path = os.path.join(OUTDIR, "rfd3_theozyme_only.png")
cmd.png(theozyme_path, dpi=600)
print("saved %s" % theozyme_path)

print("Done: %d trajectory snapshots + final-cartoon + theozyme-only in %s" % (len(frames), OUTDIR))
print("Object was left on the theozyme-only view -- re-run rfd3_movie() for the "
      "normal spheres+cloud look back if you want to `mplay` it.")

# --- Optional: if you also have a real, sequence-assigned "clean" output ---
# model of this design (not just the placeholder-sequence trajectory), you
# can identify its catalytic residues by coordinate-matching them to this
# trajectory's fixed motif (docs/RFD3_FIGURES.md #6), which is far more
# precise than a raw distance cutoff off the ligand:
#
#   cmd.load("/path/to/clean_output_model.pdb", "clean")
#   catres = catalytic_sel_from_motif("clean", fixed, max_match=4.0)
#   if catres:
#       cmd.show("sticks", catres)
#       cmd.color("orange", "(%s) and elem C" % catres)
#
# And if the design has an explicit covalent bond PyMOL won't perceive on its
# own (e.g. a TS-adduct linking a catalytic Lys to the ligand), force it once
# up front -- never inside a per-frame loop (docs/RFD3_FIGURES.md #4):
#
#   add_custom_bond("clean and resi 75 and name NZ", "clean and resn LIG and name C1",
#                    label="Lys75-NZ to ligand C1 (TS adduct)")
