# Command reference

**Not sure what's here?** Type `pymolrc_help` in PyMOL for an index of every command
group, or `help_crystal` / `help_rfd3` for a specific kit's full argument reference.

| Command | Description |
|---|---|
| `theme [list / next / prev / name]` | browse palettes or apply one to this session |
| `rfd4_autoload on/off` | switch annotation inspection for subsequent loads |
| `seq [sel-or-file]` | print the one-letter sequence of a selection/object, or a structure file path (`seq ~/x.pdb` loads it temporarily); prints usage on a bad input. Aliases: `sequence`, `get_seq`, `get_sequence` |
| `show_metals [sel]` | show metals as spheres |
| `show_catres [sel]` · `only_catres` | catalytic-residue sticks |
| `style_all [sel]` | (re)apply the full style |
| `color_palette [sel]` | align objects and apply the legacy foliage cycle and ligand shell |
| `bg_white` · `bg_black` | background toggle |
| `align_all [ref]` · `center_all [ref]` | superpose / co-center all objects |
| `autostyle` · `autoalign` · `autosolo` `on\|off` | load-time behavior toggles |
| `publication_ray_trace` | ray-traced figure render |
| `color_bb_rfdiffusion [sel]` | RFdiffusion gradient (dark blue → navaho), N→C |
| `color_bb_rfdiffusion3 [sel]` | RFdiffusion3 gradient (pink → purple → teal → dark blue), N→C |
| `color_bb_rfdiffusion4_proteina [sel]` | RFD4-Proteina Northern Nights gradient (teal → blue → violet → purple), N→C |
| `gaussian_mode [sel]` · `gaussian_off [sel]` | GaussView / QM look (glossy ball-and-stick, bond orders, perspective), and restore |
| `gaussian_spin` `on\|off` | gentle continuous spin (GUI) |
| `rfd3_movie <traj.cif.gz>` | RFdiffusion3 diffusion movie from one trajectory: diffusing protein as CA spheres + V-cloud, fixed motif/cofactor auto-detected and styled via `style_fixed`; fully configurable — `color_scheme=` takes the default RFd3 gradient, a flat custom color, or a space-separated custom multi-color gradient; `ca_scale=`/`cloud_scale=`/`cloud_transparency=` size knobs; `bonds=`/`unbonds=` lists to force/remove bonds on load; any other keyword (`trim_atoms=`, `color_map=`, `lig_rep=`, ...) passes straight through to `style_fixed`; plus `reverse=1`, `fixed_sidechain=1`, JSON cross-check (`fixed_json=`/auto); auto-runs `rfd3_check_frames` and prints whether `reverse=` actually plays noise → folded for that file; then `mplay` |
| `rfd3_check_frames <obj> [, reverse]` | which end of a loaded RFd3 trajectory is noise vs folded, measured from geometry (CA radius from centroid — noise is collapsed, folded is expanded), not assumed from the `_noisy_model_`/`_denoised_model_` filename (`rfd3_movie` doesn't treat the two differently); reports whether a given `reverse=` setting is correct for that file |
| `style_fixed <obj>` | solid, opaque ball-and-stick on an `rfd3_movie` trajectory's fixed theozyme; fixes the "fixed atoms never look opaque" problem at its root (see `docs/RFD3_FIGURES.md`); configurable — `trim_atoms=` (dict of `{resn: atom_names}` to trim any fully-fixed residue beyond the built-in GLU/HIS, or `{}` to disable trimming), `color_map=` (per-residue/selection color overrides applied over `fixed_color`), `lig_rep=` (`'sticks'` default, or `'lines'`/etc); idempotent, safe to call again any time to restyle without reloading |
| `apply_camera <obj> [, mode] [, sel] [, custom_view] [, span_states] [, turns]` | one entry point for camera control: `mode="orient"` (auto-frame), `"zoom"` (crop on `sel`), `"view"` (exact `get_view` matrix); `span_states=1` takes one fixed camera from the folded frame for a whole trajectory (RFd3 noise is a collapsed blob, not an expansion — see `docs/RFD3_FIGURES.md`) |
| `draw_connectors <obj> [, color] [, max_dist]` | dashed lines from each fixed motif residue to its nearest scaffold Cα, on a folded frame — a visual stand-in for "this theozyme residue sits here" when there's no covalent bond to draw |
| `catalytic_sel_from_motif <clean_obj>, <motif_sel> [, max_match]` | identify a clean (real-sequence) model's catalytic residues by coordinate-matching them to a trajectory's fixed motif — precise even where a ligand-distance cutoff would over-select |
| `add_custom_bond <sel_a>, <sel_b> [, label]` | force a bond PyMOL won't perceive on its own (e.g. a covalent TS-adduct); writes topology once, holds across every state |
| `remove_bond <sel_a>, <sel_b> [, label]` | remove a bond PyMOL perceived that you don't want — the mirror of `add_custom_bond` |
| `xtal <path> [, design=]` | crystallography everything-button: load crystals + sibling MTZ maps, split subunits, superpose all subunits (and the design) onto a reference, style, highlight `REMARK 666` catalytic residues, and (with `focus=active`) zoom to the pocket — see [`crystallography/`](../crystallography/) |
| `xtal_density <sel> [, level] [, carve] [, diff] [, dlevel]` | 2Fo-Fc mesh (blue) + Fo-Fc difference (green +/red −) around a subunit/selection; density tracks the aligned atoms. Tighter: raise `level` / lower `carve`; larger: the reverse |
| `xtal_rms` · `xtal_align` · `xtal_catres` · `xtal_focus` | RMSD table (every subunit vs reference and vs design); re-superpose; highlight catalytic residues; zoom pocket vs all. `help_crystal` for full args |
| `pymolrc_help` · `help_crystal` · `help_rfd3` | print a command index / a kit's full argument reference at the `PyMOL>` prompt |
| `rfd4_load <annotated.cif>, <object>` | inspect coordinate masks on cartoons with small target spheres, reference context, and chain centroids |
| `rfd4_style <object>, <mode>` · `rfd4_info <object>` · `help_rfd4` | switch condition colors or inspect annotations; see the [inspection guide](../RFD4_QUICKSTART.md) |
| `rfd4_centroids <object>, 0` / `rfd4_centroids <object>, 1` | hide/show annotated centroid targets; the choice survives subsequent style calls |
| `rfd4_cache_info` | show lightweight reads, disk/memory cache hits, native parses, and last reader timing |

The `color_bb_*` commands default to `chain A` and recolor carbons only (non-carbon atoms
left alone); pass `all_atom=1` to recolor every atom, or `backbone_only=1` to keep
sidechains as they are.

The RFD4-Proteina gradient uses the manuscript's **Northern Nights**
`molecular_spectrum`: teal `#186C70` → blue `#3B5A8D` → violet `#463D73` →
purple `#8D52AD`. These stops come from `paper/style/theme.json` at revision
`f336faf492d3`, the manuscript revision cited by the
[website stylesheet](https://github.com/RosettaCommons/RFD4-Proteina/blob/main/website/assets/site-colors.css),
and match the [preprint's](https://rfd4-proteina.org/assets/paper/RFD4-Proteina.pdf)
molecular palette. The website's interactive viewers use a flat teal instead.
Colors are interpolated by residue number, as in `color_bb_rfdiffusion`.
The individual colors are also available as `rfd4p_teal`, `rfd4p_blue`,
`rfd4p_violet`, and `rfd4p_purple`.

```pml
color_bb_rfdiffusion4_proteina chain A
color_bb_rfdiffusion4_proteina chain A, backbone_only=1
color_bb_rfdiffusion4_proteina my_design, all_atom=1
```
