# pymolrc

My PyMOL startup config (`.pymolrc`): a foliage-green color scheme plus performance
fixes and custom commands for browsing protein-design outputs.

![preview](preview_original_plus_fixes.png)

## Install

```bash
cp .pymolrc ~/.pymolrc        # back up any existing ~/.pymolrc first
# or symlink:  ln -s "$PWD/.pymolrc" ~/.pymolrc
```

The bundled `show_termini.py` sits next to `.pymolrc` and loads automatically.

## Features

- **Foliage color scheme**, applied per object on load: greens for protein, amethyst
  ligands, soft element colors, black background.
- **Metals as spheres**, detected by element, residue name, or atom name — including
  metals buried inside a ligand residue.
- **Catalytic residues** highlighted from `REMARK 666` (sidechain sticks; whole residue
  for PRO/GLY).
- **Fast batch browsing**: load a glob of designs and page through them with pgup/pgdn.
  `autosolo` loads extras hidden for near-instant startup on large sets.
- **Structural align-all**, background toggle, one-line sequence dump, and residue-gradient
  coloring.

## Commands

| Command | Description |
|---|---|
| `seq [sel]` | print the one-letter sequence of a selection |
| `show_metals [sel]` | show metals as spheres |
| `show_catres [sel]` · `only_catres` | catalytic-residue sticks |
| `style_all [sel]` | (re)apply the full style |
| `color_palette [sel]` | apply the foliage look |
| `bg_white` · `bg_black` | background toggle |
| `align_all [ref]` · `center_all [ref]` | superpose / co-center all objects |
| `autostyle` · `autoalign` · `autosolo` `on\|off` | load-time behavior toggles |
| `publication_ray_trace` | ray-traced figure render |
| `color_bb_rfdiffusion [sel]` | RFdiffusion gradient (dark blue → navaho), N→C |
| `color_bb_rfdiffusion3 [sel]` | RFdiffusion3 gradient (pink → purple → teal → dark blue), N→C |
| `gaussian_mode [sel]` · `gaussian_off [sel]` | GaussView / QM look (glossy ball-and-stick, bond orders, perspective), and restore |
| `gaussian_spin` `on\|off` | gentle continuous spin (GUI) |
| `rfd3_movie <traj.cif.gz>` | RFdiffusion3 diffusion movie from one trajectory: diffusing protein as CA spheres + V-cloud in the RFd3 gradient; fixed motif/cofactor auto-detected, held static (cofactor in your ligand style; `fixed_color` orange on fixed-protein carbons only); options `reverse=1`, `fixed_sidechain=1`, JSON cross-check (`fixed_json=`/auto); auto-runs `rfd3_check_frames` and prints whether `reverse=` actually plays noise → folded for that file; then `mplay` |
| `rfd3_check_frames <obj> [, reverse]` | which end of a loaded RFd3 trajectory is noise vs folded, measured from geometry (CA radius from centroid — noise is collapsed, folded is expanded), not assumed from the `_noisy_model_`/`_denoised_model_` filename (`rfd3_movie` doesn't treat the two differently); reports whether a given `reverse=` setting is correct for that file |
| `style_fixed <obj>` | solid, opaque ball-and-stick on an `rfd3_movie` trajectory's fixed theozyme; fixes the "fixed atoms never look opaque" problem at its root (see `docs/RFD3_FIGURES.md`) and grades stick/sphere size so real backbone reads clearly against placeholder sidechain atoms |
| `apply_camera <obj> [, mode] [, sel] [, custom_view] [, span_states] [, turns]` | one entry point for camera control: `mode="orient"` (auto-frame), `"zoom"` (crop on `sel`), `"view"` (exact `get_view` matrix); `span_states=1` takes one fixed camera from the folded frame for a whole trajectory (RFd3 noise is a collapsed blob, not an expansion — see `docs/RFD3_FIGURES.md`) |
| `draw_connectors <obj> [, color] [, max_dist]` | dashed lines from each fixed motif residue to its nearest scaffold Cα, on a folded frame — a visual stand-in for "this theozyme residue sits here" when there's no covalent bond to draw |
| `catalytic_sel_from_motif <clean_obj>, <motif_sel> [, max_match]` | identify a clean (real-sequence) model's catalytic residues by coordinate-matching them to a trajectory's fixed motif — precise even where a ligand-distance cutoff would over-select |
| `add_custom_bond <sel_a>, <sel_b> [, label]` | force a bond PyMOL won't perceive on its own (e.g. a covalent TS-adduct); writes topology once, holds across every state |

The `color_bb_*` commands default to `chain A` and recolor carbons only (non-carbon atoms
left alone); pass `all_atom=1` to recolor every atom, or `backbone_only=1` to keep
sidechains as they are.

## RFdiffusion3 trajectory figures

**Just want the RFd3 movie/figure commands, not this whole foliage color
scheme?** See **[`RFD3_QUICKSTART.md`](RFD3_QUICKSTART.md)** — a standalone
page for sharing with labmates. It covers grabbing just
[`rfd3_movie_kit.pml`](rfd3_movie_kit.pml) (one file, no other dependency) and
either `run`-ing it from your own `.pymolrc` or copy-pasting it in directly,
plus a command reference and a quickstart example.

**Use the `_noisy_model_*.cif.gz` trajectory, not the `_denoised_model_*.cif.gz`
one, if you have a choice.** `rfd3_movie` loads and renders either one the same
way, and both are supported — but the noise/folded frame-order convention
(which drives the `reverse=` default) has only been verified against noisy
trajectories. `rfd3_check_frames` runs automatically on every `rfd3_movie` call
and prints a warning if `reverse=` looks wrong for whatever file you loaded, so
a denoised trajectory will still work correctly — just double check that
printed warning if you use one.

`docs/RFD3_FIGURES.md` collects the hard-won principles behind building
Figure-1-style stills and storyboard movies from an RFd3 trajectory on top of
the commands above — frame vs. state indexing, why disconnected motif
fragments ignore transparency, why RFd3 noise is a collapsed blob (not an
expanding cloud) and what that means for framing a camera, why you must never
toggle representations per-frame on the live trajectory object, GIF/MP4
assembly pitfalls, and more. Read it before re-debugging any of these.

`examples/make_rfd3_figure.py` and `examples/make_rfd3_gif.py` are complete,
runnable worked examples built on the bundled trajectory in `example_pdbs/`
(`pymol -cq examples/make_rfd3_figure.py`). `scripts/frames_to_movie.py` is
the pure-Python (no PyMOL) post-processing step they hand off to for
GIF/MP4 assembly — usable standalone once you have rendered PNG frames.

## Examples

`example_pdbs/` holds sample structures to try the commands on:
- `ZETA_1__A1_metalloesterase_theozyme.pdb` — a metalloesterase theozyme (Zn + His triad + substrate); good for `gaussian_mode` and `color_bb_*`.
- `ZAPP_p1D1_i14_rfd3_noisy_trajectory.cif.gz` — an RFdiffusion3 diffusion trajectory for `rfd3_movie` and the `examples/` scripts above.

## Reverting

The original pre-fix config is kept locally (git-ignored) as
`pymolrc_backup_2026-06-12.pymolrc`. Restore it with:

```bash
cp ~/pymolrc/pymolrc_backup_2026-06-12.pymolrc ~/.pymolrc
```
