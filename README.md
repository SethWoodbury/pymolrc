# pymolrc

My PyMOL startup config (`.pymolrc`) — a foliage-green color scheme plus performance
fixes and custom commands for browsing protein-design outputs: metals as spheres,
catalytic-residue highlighting from `REMARK 666`, fast batch loading, structural
align-all, sequence printing, and background toggles. Full changelog below.

## Install
```bash
cp .pymolrc ~/.pymolrc        # back up any existing ~/.pymolrc first
# or symlink:  ln -s "$PWD/.pymolrc" ~/.pymolrc
```

> On load, the config runs the bundled `show_termini.py` (marks protein N- and
> C-termini with colored spheres). It lives in this repo next to `.pymolrc`, so a
> fresh clone gets the command automatically.

---

# ~/.pymolrc — surgical fixes (2026-06-12, final)

Your **original aesthetic is preserved byte-for-byte** (foliage greens, amethyst
ligand, element colors, render/lighting settings, black background, `color_palette`).
Only the functional fixes you asked for were added.

## Backup (local only, not committed)
`pymolrc_backup_2026-06-12.pymolrc` — the ORIGINAL pre-fix config, kept next to this
repo as a revert target. It's git-ignored, so it stays local and isn't pushed to GitHub.

## Revert
```bash
cp ~/pymolrc/pymolrc_backup_2026-06-12.pymolrc ~/.pymolrc
```

## What was changed (and ONLY this)
1. **Lag fix** — removed the `cmd.load`/`cmd.fetch` monkeypatch that ran
   `center('all')` on every load (O(N²)). Replaced with a lean per-object hook.
2. **pgup/pgdn** — removed the `orient()` in `structure_step`; cycling no longer
   resets the camera.
3. **Backbone unstuck** — `color_palette` now `unset cartoon_color` instead of
   pinning it, so `color red, myobj` works on the cartoon.
4. **Metals → spheres**, detected by **element OR residue name OR atom name**
   (atom-name match is guarded with `not polymer.protein` so the alpha-carbon
   "CA" is never mistaken for calcium). Spheres a bit bigger (`sphere_scale 0.7`).
   Works on metals embedded in a ligand residue (e.g. the test Zn lives in `PSZ`).
5. **REMARK 666 catalytic residues** — parsed from the loaded PDB
   (chain = field[9], resi = field[11] on the `MATCH MOTIF` side), shown as sticks
   in a distinct highlight color (`catres_hi`, orange). Only the **sidechain**
   carbons are highlighted; the catalytic **backbone** keeps the protein color —
   EXCEPT for PRO/GLY, whose backbone carbons are highlighted (no usable sidechain).
   Styled PER OBJECT on load (and by `color_palette`). There is **no persistent
   `catres` selection** — one aggregate selection OR'd across many objects did not
   scale to big multi-file/glob loads (looked like "only the first one"), so styling
   is applied per object instead. `show_catres` / `only_catres` also work per object.
6. **bg toggle** — `bg_white` / `bg_black` (black stays the default).
7. **`seq [sel]`** — prints one continuous one-letter sequence of `sele` (or any
   selection); His triad → `HHH`.

## New commands
- View/style: `seq [sel]` · `show_metals [sel]` · `show_catres [sel]` · `only_catres` ·
  `style_all [sel]` · `bg_white` · `bg_black`
- Align: `align_all [ref]` (structural superpose all → ref; `super` + `cealign` fallback) ·
  `center_all [ref]` (instant translate-only co-centering)
- Color: `color_bb_rfdiffusion [selection]` — apply the RFdiffusion gradient (dark blue → navaho)
  N→C across residues. Defaults to `chain A`, carbons only (non-C atoms untouched).
  Options: `all_atom=1` (recolor every element) · `backbone_only=1` (keep sidechains as-is).
- Color: `color_bb_rfdiffusion3 [selection]` — RFdiffusion3 gradient (pink → purple → teal → dark blue)
  N→C, teal auto-centered at the chain midpoint. Same defaults/options as above.
- Toggles: `autostyle on|off` · `autoalign on|off` · `autosolo on|off`
- Unchanged: `color_palette`, `structure_*`, `publication_ray_trace`, …

## Align loaded structures into the same space
After loading a batch/glob of diverse structures:
- `align_all` — structurally superpose every object onto the first (or `align_all <ref>`).
  Each object aligned once to the reference (O(N)); ~0.3 s for 30 structures.
- `center_all` — if you just want them roughly co-located fast (translation only).
- `autoalign on` — superpose each new object onto the first AS it loads (hands-free).

## Faster loading (measured on 30 copies of a1_hit)
The old `center('all')`-on-every-load was O(N²) (the "laggy after a while" bug) — removed.
Per-object styling was also O(N²) (each `show`/`set` scales with #objects loaded). Fixes:
global `sphere_scale`/`stick_radius` settings (no costly per-selection `set`), merged
color ops, and **lazy styling**.
- **Default** (`autostyle on`, `autosolo off`): every object styled+shown on load — ~2.2 s/30.
- **Fast browse** (`autosolo on`): extras load HIDDEN and UNstyled — **~0.23 s/30** — then
  pgup/pgdn shows one at a time and styles each on first view; `structure_show_all` or
  `style_all` styles everything. Best for big globs. Make it permanent by adding a line
  `autosolo on` near the end of `~/.pymolrc`.

## Change the catalytic highlight color
Defaults to orange (`catres_hi`). To change: edit `set_color('catres_hi', [r,g,b])` in the
rc (then `style_all` to re-apply), or color specific residues directly, e.g.
`color yellow, <obj> and resi 118+130+134 and sidechain`.

Preview of your colors + the fixes on `a1_hit.pdb`: `preview_original_plus_fixes.png`.

## Update — lazy loading + `.bashrc` (much faster `random_shuffle`)
Profiling `random_shuffle 10` on the real design dir: the bottleneck was the `.bashrc`
`pymol()` wrapper forcing `-d "color_palette"` on every launch (~2 s of super-align +
clash scan + ligand-shell) which also re-enabled every object (cancelling autosolo).

- **`.bashrc` `pymol()`** no longer forces `color_palette`. **`random_shuffle`** now uses
  lazy browse. Run `source ~/.bashrc` (or open a new shell) to pick up the changes.
- **Foliage coloring moved into auto-styling** — loading any structure now applies the
  full look (greens + amethyst ligand + element colors + ligand sticks + metal spheres +
  catalytic sticks) per object, so `color_palette` is no longer needed just for the look.
- **AF3 clash-removal OFF by default** — `clash_bond_filter on` (and `ligand_center_filter on`)
  to re-enable for AlphaFold3 / diffusion outputs with spurious bonds.
- **`browse_list <listfile>`** (used by `random_shuffle`): loads ALL structures (so they
  all appear in the object list) and shows **ONLY the first**. **Thread-free** — a
  background loader was tried but PyMOL does load+render on ONE thread, so a worker thread
  just freezes the GUI (and raced, leaving the last one shown); true rotate-while-loading
  is impossible. The unavoidable styling pause is made a single predictable startup pause:
  - **`browse_prestyle on`** (DEFAULT): pre-style everything up front — one-time pause
    (~1.5 s for 10 real designs), then pgup/pgdn is instant.
  - **`browse_prestyle off`**: instant start — first shows fast, the rest load light and
    style on first view (~0.2 s each, instant after; grows on very large sets).
  - **`style_all`** (batched: one global show per representation, not per-object) styles
    everything on demand.

Measured on 10 real designs: pre-style-all startup **~1.5 s** then instant scroll; or
instant-start **~0.2 s** with ~0.2 s on first view of each. Only the first is ever shown.

## Optimization round (committee-reviewed) — O(N²) → O(N)
A 3-agent review + benchmarks found styling was O(N²) (each per-object `show`/`color` is
re-evaluated across the whole scene). Fixes:
- **`style_all` fully batched:** ONE global op per color/representation + catalytic residues
  via a **chunked transient combined selection** (object-scoped clauses, residues grouped
  per chain, built→used→deleted; avoids the stale-aggregate-selection bug) + the whole batch
  under `suspend_updates` (try/finally). Single foliage green on the fast path (you view one
  at a time); `color_palette` keeps cycling greens for overlays.
- **Load hook O(N²)→O(N):** a persistent `_KNOWN_OBJS` set → ONE `get_object_list()` per load
  (was 3); a `delete` hook prunes it so delete+reload re-styles correctly; `ref` passed in.
- **REMARK 666 cached** by (mtime,size) → no re-read on reload.
- **Direct globs routed to the fast path:** the `.bashrc` `pymol()` wrapper sends a many-file
  (>4, no flags) glob like `pymol .../*scaffold_5*.pdb` through `browse_list` (load all,
  batch pre-style, show first). 2–3 files or any flag → normal overlay as before.

**Measured on 25 real scaffold_5 designs: 9.8 s → 0.855 s (~11×).** All loaded + pre-styled +
only first shown. Catalytic backbone rule, metals, delete+reload all verified intact.
Rejected by the review (would break things): `defer_builds_mode` (reintroduces pgup lag),
filename-derived object names (breaks on multiplex/dedup), a persistent combined `catres`
selection (the old "only the first" bug), changing `connect_mode` (risks metal/ligand bonds).
PyMOL's parser is single-threaded — per-file parse (~15 ms) is already at its floor.

CPU / memory (honest): PyMOL's parser is single-threaded — one structure's load can't be
split across cores. The free wins used here: lazy loading (don't parse what you don't
view), background I/O prefetch (overlaps NFS latency), RAM caching of seen structures,
and removing redundant heavy work (`color_palette`/clash on every launch). For a big
plain `pymol *.pdb` glob, use `autosolo on` or `browse_list`.
