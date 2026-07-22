# RFdiffusion3 movie/figure kit — quickstart

Everything you need to load an RFdiffusion3 diffusion trajectory in PyMOL and
turn it into Figure-1-style stills or a looping GIF/MP4 — without adopting the
rest of this repo's (foliage-green) color scheme. One file, no dependency on
`.pymolrc`.

**Use the `_noisy_model_*.cif.gz` trajectory if you have a choice** — that's
what this kit was built and tested against. The `_denoised_model_*.cif.gz`
trajectory also works, but hasn't been verified to follow the same
noise/folded frame-order convention; see the note under Quickstart below.

## Install — pick one

**Option A: download the one file, `run` it from your own `.pymolrc`**

```bash
curl -o ~/rfd3_movie_kit.pml https://raw.githubusercontent.com/SethWoodbury/pymolrc/main/rfd3_movie_kit.pml
```

Then add one line to the end of your own `~/.pymolrc`:

```
run ~/rfd3_movie_kit.pml
```

(Or skip editing `.pymolrc` and just type `run ~/rfd3_movie_kit.pml` at the
`PyMOL>` prompt whenever you want it for that session.)

**Option B: copy-paste**

Open [`rfd3_movie_kit.pml`](rfd3_movie_kit.pml) and paste its contents
straight into your own `.pymolrc`, anywhere after your other custom commands.

Either way, restart PyMOL (or `run` the file) and you have 8 new commands:
`rfd3_movie`, `rfd3_check_frames`, `style_fixed`, `apply_camera`,
`draw_connectors`, `catalytic_sel_from_motif`, `add_custom_bond`,
`color_bb_rfdiffusion3`.

## Quickstart

```
PyMOL> run ~/rfd3_movie_kit.pml
PyMOL> rfd3_movie /path/to/your_noisy_model_2.cif.gz, reverse=1
PyMOL> style_fixed rfd3_traj
PyMOL> mplay
```

Want to try it on a real trajectory first? This repo ships one:
`example_pdbs/ZAPP_p1D1_i14_rfd3_noisy_trajectory.cif.gz`. Grab just that file
plus the kit:

```bash
curl -o ~/rfd3_movie_kit.pml https://raw.githubusercontent.com/SethWoodbury/pymolrc/main/rfd3_movie_kit.pml
curl -o ~/rfd3_example.cif.gz https://raw.githubusercontent.com/SethWoodbury/pymolrc/main/example_pdbs/ZAPP_p1D1_i14_rfd3_noisy_trajectory.cif.gz
```

```
PyMOL> run ~/rfd3_movie_kit.pml
PyMOL> rfd3_movie ~/rfd3_example.cif.gz, reverse=1
PyMOL> style_fixed rfd3_traj
PyMOL> mplay
```

`reverse=1` is worth using by default — on real trajectories seen so far, raw
file state 1 is the folded structure and state N is noise, so `reverse=1` makes
frame 1 show noise and the last frame show the folded design. **`rfd3_movie`
doesn't actually know or care whether your file is a `_noisy_model_` or
`_denoised_model_` trajectory** — it loads and renders either one identically,
and only ever uses the filename to auto-locate the design JSON. So every call
auto-runs `rfd3_check_frames` and prints a line telling you whether your
`reverse=` setting is actually correct for *that specific file* — read that
line rather than assuming `reverse=1` is always right.

## Command reference

| Command | Description |
|---|---|
| `rfd3_movie <traj.cif.gz>` | Load one trajectory. Diffusing protein shows as CA spheres + a faint sidechain cloud in the RFd3 gradient; the fixed catalytic motif/cofactor/metal are auto-detected and held static, colored `fixed_color` (default orange) on fixed-protein carbons. Options: `reverse=1` (see above), `color_scheme=0` (plain cyan instead of gradient), `cloud=0` (CA spheres only, no sidechain cloud), `fixed_sidechain=1` (hide fixed backbone, sidechains only), `fixed_json=` (cross-check against the design JSON's `select_fixed_atoms`, else auto-derived from the trajectory filename). Then `mplay`. |
| `rfd3_check_frames <obj> [, reverse]` | Sanity-check which end of a loaded trajectory is noise vs folded, from the actual geometry (CA radius from centroid of the diffusing chain — noise is a collapsed blob, folded is expanded) rather than assumed from the filename. Runs automatically at the end of every `rfd3_movie` call; run it again by hand any time you want to double-check, e.g. after loading a file some other way, or on a `_denoised_model_` trajectory you're not sure follows the same convention. |
| `style_fixed <obj>` | Solid, opaque ball-and-stick on the fixed theozyme. Fixes the classic "fixed atoms look see-through" problem at the root (disconnected motif fragments render as `nonbonded`/`lines`, which ignore transparency settings entirely — see `docs/RFD3_FIGURES.md` item 2) and sizes real backbone atoms bigger than the placeholder sidechain atoms so the two read clearly apart. |
| `apply_camera <obj> [, mode] [, sel] [, span_states] [, n_states] [, turns]` | One entry point for camera control: `mode="orient"` (default, auto-frame), `mode="zoom"` (crop on `sel`), `mode="view"` (exact `get_view` matrix via `custom_view=`). `span_states=1, n_states=<n>` takes one fixed camera from the folded frame and holds it across the whole trajectory — needed because RFd3 noise is a *collapsed* blob, not an expanded one (the folded frame is the largest state and the one that can clip; see `docs/RFD3_FIGURES.md` item 3). |
| `draw_connectors <obj> [, color] [, max_dist]` | Dashed lines from each fixed motif residue to its nearest scaffold Cα — a stand-in for "this theozyme residue sits here" when there's no covalent bond to draw. Only meaningful on a folded frame. |
| `catalytic_sel_from_motif <clean_obj>, <motif_sel> [, max_match]` | If you also have a real, sequence-assigned output model of the design (not just the placeholder-sequence trajectory), identifies its catalytic residues by coordinate-matching them to the trajectory's fixed motif — far more precise than a distance cutoff off the ligand. Returns a selection string. |
| `add_custom_bond <sel_a>, <sel_b> [, label]` | Force a bond PyMOL won't perceive on its own (e.g. a covalent TS-adduct). Each selection must resolve to exactly one atom; writes topology once, so it holds across every state. |
| `color_bb_rfdiffusion3 [sel]` | The RFdiffusion3 gradient (pink → purple → teal → dark blue) N→C on any selection — what `rfd3_movie` uses internally, callable on its own too. |

## Watch out for

- **Never toggle representations (e.g. `show cartoon`) per-frame on the live
  trajectory object** — it hard-hangs PyMOL, unrecoverable with Ctrl-C. If you
  want a cartoon on just the final folded frame, pull that one state out first
  with `cmd.create(tmp, obj, source_state, 1)` and style the copy instead.
- `reverse=1` only flips the movie *frame*→state mapping. Anything that takes
  a literal state number (`cmd.create`'s `source_state`) still needs the real
  state, not the frame number.

Full write-up of these and more (GIF/MP4 assembly, why noise is collapsed not
expanded, why `.pymolrc` functions are called bare instead of via `cmd.*`,
etc.) is in [`docs/RFD3_FIGURES.md`](docs/RFD3_FIGURES.md) — worth a read
before re-debugging any of it. Full runnable examples (including the
GIF/MP4-assembly step) are in [`examples/`](examples/), built on
`scripts/frames_to_movie.py` for the post-processing.
