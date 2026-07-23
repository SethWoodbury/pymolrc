# Crystallography kit for PyMOL — quickstart

Custom PyMOL commands for working up refined crystal structures: load a folder
of crystals with their maps, split every subunit into its own object, superpose
all subunits onto each other **and** onto a design model, show 2Fo-Fc / Fo-Fc
electron density that follows the aligned atoms, print an RMSD summary, and
highlight the design's REMARK 666 catalytic residues — all from one command.

Self-contained: one file (`crystal_kit.pml`), no dependency on the rest of this
repo or on `.pymolrc`. Copy just this `crystallography/` folder and you have the
kit plus a runnable demo.

## Install — pick one

**Option A: download the one file, `run` it from your own `.pymolrc`**

```bash
curl -o ~/crystal_kit.pml https://raw.githubusercontent.com/SethWoodbury/pymolrc/main/crystallography/crystal_kit.pml
```

Then add one line to the end of your own `~/.pymolrc`:

```
run ~/crystal_kit.pml
```

(Or just type `run ~/crystal_kit.pml` at the `PyMOL>` prompt when you want it.)

**Option B: copy-paste** the contents of `crystal_kit.pml` into your own
`.pymolrc`, anywhere after your other custom commands.

Either way you get seven commands: `xtal`, `xtal_align`, `xtal_density`,
`xtal_rms`, `xtal_catres`, `xtal_focus`, `help_crystal`.

## Quickstart — try it on the bundled demo

This folder ships two real crystals of the same design (`ZAPP_P1D1_Pin5` and
`ZAPP_P1D1_Pin16`, each with its refined `.mtz`) plus the design model
(`ZAPP_P1D1_design.pdb`, which carries six REMARK 666 catalytic residues and a
YYE transition-state analog).

```
PyMOL> run ~/pymolrc/crystallography/crystal_kit.pml
PyMOL> xtal ~/pymolrc/crystallography/example, design=~/pymolrc/crystallography/example/ZAPP_P1D1_design.pdb, focus=active
PyMOL> xtal_rms
PyMOL> xtal_density xtal_Pin5_A
```

That one `xtal` call loads both crystals + their maps, splits each into
per-chain subunits (`xtal_Pin5_A..D`, `xtal_Pin16_A..D`), superposes all eight
subunits onto a reference and onto the design, styles everything (crystals
colored per-crystal, design in green with an amethyst ligand), highlights the
catalytic residues in magenta, and zooms to the active-site pocket.

## The everything-button

```
xtal <path> [, design=] [, split=1] [, align=1] [, to_design=1] [, catres=1]
            [, focus=all] [, surface_metals=1]
```

**`path`** (the only required argument) can be a directory of `.pdb` files, a
glob, or a single `.pdb`. For each crystal it looks for a sibling `<name>.mtz`
and loads its maps. Sensible defaults do the whole workup; every step is also a
standalone command you can re-run with different options:

| Option | Default | Meaning |
|---|---|---|
| `design=` | *(none)* | Path to a design/reference model. Its REMARK 666 residues become the catalytic highlight. If it sits in the same folder as the crystals it is **not** double-loaded as a crystal. |
| `split=` | `1` | Split each crystal into per-chain subunit objects (each keeps its own active-site metals). |
| `align=` | `1` | Superpose all subunits onto a reference (first subunit). |
| `to_design=` | `1` | Also superpose the design onto the reference (by sequence, so catalytic residues register correctly). |
| `catres=` | `1` | Highlight the design's REMARK 666 residues in every subunit. |
| `focus=` | `all` | `active` zooms straight to the catalytic pocket; `all` frames everything. |
| `surface_metals=` | `1` | Keep all metals (useful). `0` hides metals outside the active site. |

## Electron density — what to type, and what it is

```
xtal_density <sel> [, level=1.0] [, carve=1.8] [, diff=1] [, dlevel=3.0] [, surface=0]
```

Draws the **2Fo-Fc** map (blue mesh) around `sel` at `level` sigma — this is
the standard "where is the model supported by data" map — and, with `diff=1`,
the **Fo-Fc difference** map at ±`dlevel` sigma: **green** = positive
(density present, nothing modeled — a missing atom/ligand/water) and **red** =
negative (something modeled with no density under it). Levels are in sigma
(`normalize_ccp4_maps` is turned on for you).

The map is copied into the subunit's aligned frame before contouring, so the
density **tracks the atoms even after alignment** — you can overlay density from
both crystals on the same superposed pocket.

- **Tighter** density: raise `level` (e.g. `1.5`) and/or lower `carve` (`1.2`).
- **Larger / more context**: lower `level` (`0.8`) and/or raise `carve` (`2.5`).
- `surface=1` draws a translucent solid isosurface instead of a mesh.

```
PyMOL> xtal_density xtal_Pin5_A, level=1.5, carve=1.2      # tight
PyMOL> xtal_density xtal_Pin5_A and resi 89+92, carve=2.5  # just two residues, roomy
```

## The rest

| Command | What it does |
|---|---|
| `xtal_align [, reference=] [, to_design=1]` | Re-superpose all subunits onto a reference (default: first subunit); `super` with a `cealign` fallback. |
| `xtal_rms [, reference=] [, design=]` | Print an RMSD table: every subunit vs the reference and vs the design (chain-agnostic — NCS copies A/B/C/D pair up correctly), plus per-crystal Zn/ligand info. |
| `xtal_catres [, design=] [, color=magenta] [, rep=sticks] [, label=0]` | Highlight the design's REMARK 666 catalytic residues in every subunit + the design. Change the color/representation, or `label=1` for `resn+resi` labels. |
| `xtal_focus [, target=active]` | Zoom: `active` = catalytic pocket (default), `all` = everything, or any selection string. |
| `help_crystal` | Print the full command reference at the `PyMOL>` prompt. |

## Customizing

Every command is a thin wrapper over ordinary PyMOL `load`/`create`/`super`/
`isomesh`/`color` calls on plain selections, so you keep full manual override at
the prompt afterward.

**Default palette** (only carbons are recolored — every other element keeps its
CPK color, so N/O/S/metals read normally): crystals get **grey** backbones
(`cx_grey1/2/3`, cycled per crystal) with **goldenrod** catalytic carbons
(`cx_cat_xtal`); the design keeps its **foliage-green** backbone (`cx_green`)
with **aqua-cyan** catalytic carbons (`cx_cat_design`) and an **amethyst** ligand
(`cx_amethyst`). Warm-gold vs cool-cyan is what separates the crystal and design
active sites when they superpose (and it's colorblind-safe — it rides the
blue-yellow axis). All of these are `set_color`'d at the top of
`crystal_kit.pml`; override per call with e.g.
`xtal_catres color=orange, design_color=magenta`.

## Notes

- **The two original PDBs are kept as reference objects.** After splitting, each
  crystal's full, unsplit, **un-aligned** structure stays loaded (`xtal_Pin5`,
  `xtal_Pin16`) at the top of the object panel — but **disabled**, so it's out of
  view until you tick it on. Handy for checking anything against the untouched
  original in its own crystal coordinates. The `_A..D` subunits below it are the
  aligned working copies.
- **Stick and metal sizes follow your `.pymolrc`** (`stick_radius 0.25`,
  `sphere_scale 0.7`, `valence off`, metals in the `elemZn` blue-grey) so kit
  objects match the rest of your figures rather than looking thin.
- **2Fo-Fc vs Fo-Fc** are the two field-standard maps; both come straight from
  the refined `.mtz` (PyMOL auto-generates `<map>.2fofc` and `<map>.fofc` on
  load). If a `.mtz` holds only unmerged data (no map coefficients), `xtal`
  says so rather than failing silently.
- Subunit split includes **only that chain's own active-site metals**, not
  nearby protein or modified residues bleeding in from a neighboring chain — so
  a shared catalytic residue number (e.g. a carbamylated Lys, KCX 16) isn't
  double-counted.
- The design is aligned onto the crystals **by sequence** (`align`, not
  `super`), because a purely structural superposition can lock onto a shifted
  register and put the catalytic residues in the wrong place.
