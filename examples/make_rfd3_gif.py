#!/usr/bin/env python3
"""
Worked example: render every frame of an RFdiffusion3 trajectory and assemble
a looping GIF + truecolor MP4, using the consolidated .pymolrc commands plus
scripts/frames_to_movie.py for assembly.

Run from Terminal:
    pymol -cq examples/make_rfd3_gif.py

Or from inside an interactive PyMOL session:
    run examples/make_rfd3_gif.py

See docs/RFD3_FIGURES.md for why representations are never toggled per-frame
here (#4), and why GIF/MP4 assembly happens the way it does in
frames_to_movie.py (#8, #9).
"""
import os
import sys
from pymol import cmd

HERE = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else os.getcwd()
REPO = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(REPO, "scripts"))
from frames_to_movie import frames_to_gif, frames_to_mp4  # noqa: E402

TRAJ = os.path.join(REPO, "example_pdbs", "ZAPP_p1D1_i14_rfd3_noisy_trajectory.cif.gz")
OUTDIR = os.path.join(HERE, "output")
FRAMES_DIR = os.path.join(OUTDIR, "rfd3_gif_frames")
os.makedirs(FRAMES_DIR, exist_ok=True)

NAME = "rfd3_traj"
FPS = 15
RAY_TRACE = False  # True = much nicer, much slower for ~100 frames
WIDTH, HEIGHT = 900, 675

# reverse=1: see docs/RFD3_FIGURES.md #1 -- raw state 1 is the true folded
# structure, state N is true noise, backwards from what's intuitive.
rfd3_movie(TRAJ, color_scheme=1, cloud=1, fixed_color="orange", name=NAME, reverse=1)

n_states = cmd.count_states(NAME)
move = "(%s) and resn UNK" % NAME

style_fixed(NAME)
cmd.set("sphere_transparency", 0.35, "%s and name CA" % move)

cmd.bg_color("white")
cmd.set("ray_opaque_background", 1)
cmd.set("antialias", 2)
cmd.viewport(WIDTH, HEIGHT)

# One fixed camera for the whole GIF, taken from the folded (largest) frame --
# docs/RFD3_FIGURES.md #3. NOTE: representations are never touched per-frame
# below (no cartoon toggling etc.) -- see #4 for why that hangs PyMOL.
apply_camera(NAME, mode="orient", span_states=True, n_states=n_states, zoom_buffer=6.0)

frame_paths = []
for i in range(1, n_states + 1):
    cmd.frame(i)
    if RAY_TRACE:
        cmd.ray(WIDTH, HEIGHT)
    path = os.path.join(FRAMES_DIR, "frame_%03d.png" % i)
    cmd.png(path, width=WIDTH, height=HEIGHT, dpi=150, ray=int(RAY_TRACE))
    frame_paths.append(path)
    if i % 10 == 0 or i == n_states:
        print("rendered %d/%d" % (i, n_states))

print("All %d frames saved in %s" % (n_states, FRAMES_DIR))

# Hold the last (folded) frame a bit longer than the rest so it reads as the
# "destination," then hand off to the pure-Python assembly helpers -- no
# PyMOL involved from here on.
durations = [1000 // FPS] * len(frame_paths)
durations[-1] = 1500

gif_path = os.path.join(OUTDIR, "rfd3_diffusion.gif")
frames_to_gif(frame_paths, durations, gif_path)
print("GIF saved: %s" % gif_path)

try:
    mp4_path = os.path.join(OUTDIR, "rfd3_diffusion.mp4")
    frames_to_mp4(list(zip(frame_paths, durations)), mp4_path, fps=FPS)
    print("MP4 saved: %s" % mp4_path)
except RuntimeError as e:
    print("Skipped MP4 (%s) -- GIF above is still available." % e)
