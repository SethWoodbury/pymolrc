# RFD4 annotated CIFs in PyMOL

Open an annotated input exactly like a normal structure:

```bash
pymol /absolute/path/to/annotated_input.cif
```

Or use `load /absolute/path/to/annotated_input.cif` inside PyMOL. The installed
`~/.pymolrc` automatically detects native RFD4 condition fields and displays the
coordinate-conditioned atoms with the normal element palette. Protein absolute-target
carbons are cyan, substrate/cofactor carbons are amethyst, and internal-group carbons
are orange. Nitrogen, oxygen, sulfur, phosphorus and hydrogen use the named colors
from `.pymolrc`; other elements keep PyMOL's element colors. Coordinate-conditioned
atoms are sticks/small spheres; unconstrained reference atoms are lines, or crosses
when unbonded. Protein carbons in wholly sequence-free residues use a slightly
lighter shade of their primary color. No special loader
command is needed. Ordinary PDB/CIF files retain their existing loader behavior.
The same automatic handling works with `browse_list` and `random_shuffle`.

Chain-centroid targets from the source annotations appear as magenta translucent
spheres with crosses, labeled `Centroid<outputchain>`, one per conditioned chain.
These are `ORI`/`X` visual markers at the supplied target coordinates, not means
computed from the displayed atoms. No markers are created when centroid targets
are absent. Use `rfd4_centroids motif, 0` to hide them or
`rfd4_centroids motif, 1` to show them; the choice is remembered across style calls.
Markers belong to the parent object, so they follow its navigation, transforms,
and saved PyMOL session appearance. They do not change the source CIF.

Metals use the normal metal color and large spheres, inheriting the global
`sphere_scale` setting (normally 0.7), including subsequent changes to that setting.
Nonmetal annotation spheres remain small. Reference-only metals appear as crosses,
not spheres. `context=0` hides all reference atoms. Element colors identify chemistry; use the flat coordinate
mode below to inspect absolute/internal status on every atom, including metals.

Automatic detection requires explicit `_condition_*` or
`_atom_site.condition_*` fields. The `annotation_expsegmin/max` columns alone
are bookkeeping and also occur in saved predictions; they do **not** activate
the input reader. AF3 and generated CIFs with only those leftover columns use
PyMOL's standard loader, retaining ligand/cofactor atoms even when their
occupancy is undefined. Detection depends on content, not the filename.
Use explicit `rfd4_load` for an unusual instruction-only file without condition
fields. This does not alter the source CIF or repair previously loaded objects:
restart PyMOL and reopen them, or run `run ~/pymolrc/rfd4_viewer.py` and reload
the affected structures into new objects.

Review scripts can override load-time styling with `hide everything` followed
by `show sticks`. An isolated metal has no stick to draw. The updated `.pymolrc`
therefore adds spheres for selected metals when `show`/`as` requests sticks or
lines (also nonbonded representations). Restart PyMOL to activate this change;
`show_metals` remains available in an older session. This is separate from loader
detection and does not alter metal coordination or add bonds. Explicit hiding
still works, and `autostyle off` disables the automatic addition.

An annotated input can contain an expandable scaffold placeholder whose coordinates
are NaN. Ordinary PyMOL's mmCIF loader includes it, making the viewing bounds
NaN and the scene blank. The automatic adapter omits instruction placeholders
from the display; it never
changes the annotated source. Native coordinates, charges and bonds are retained.

Self-contained inputs first use a lightweight Biopython CIF tokenizer and NumPy
reader. This route needs no native scientific imports, subprocess, display cache,
or filesystem scan of the RFD4 installation. Its narrow contract requires one
model, no alternate locations, matching author/label identities, known condition
encodings, identity `atom_array_index` bookkeeping when present, explicit local
atom/bond templates, canonical protein components and
simple polymer links. It preserves atom-site formal charges. Unknown schemas,
metals, hydrogens, noncanonical polymers, physical explicit links, and other
unsupported cases fall back to the native reader. Explicit `repo=`/`python=` or
RFD4 runtime environment overrides also select the native path. This is a
display reader only; inference always uses the authoritative native parser.

The fallback native route uses a persistent, verified JSON display cache under
`~/.cache/pymol/rfd4/`. A cached file opens directly in PyMOL, including in a
fresh session, without starting the native scientific Python environment.
The source bytes, reader/condition code, parsing settings, package versions,
and local configuration determine cache validity. Changed input or reader code
causes a fresh native parse. Cached coordinates, bonds, charges, and annotation
records are checksummed and validated before use; corrupt entries are rebuilt.
Code hashes are reused only when their file identity/size/ctime/mtime match.
Native package binary/data files are tracked by their filesystem change stamps.
File renames reuse the same content cache and display the current source path.

On an unsupported, uncached file, the native CPU parser starts once, saves the resulting
display data, and is reused for subsequent uncached files. It exits after 60
seconds without requests. A 32-model in-memory cache also speeds backtracking.
No checkpoint, remote server, or GPU is used. Re-running the helper preserves
loaded annotation metadata and does not stack duplicate loader hooks.

For an input library that needs the native fallback, optionally prepare its display cache in advance
using a few persistent CPU workers:

```bash
python ~/pymolrc/scripts/rfd4_model_cache.py --source-dir /path/to/annotated --jobs 4
```

This reads only the directory's immediate `*.cif` files, keeps all source files
unchanged, and reports cache hits, newly parsed files, elapsed time, and the
reader fingerprint. It is optional: normal opens populate the same cache.
`RFD4_PYMOL_CACHE=/another/path` overrides the cache location. The entire cache is
disposable; removing it simply requires native parsing again. Neither the
cache nor a PyMOL export is a replacement for the annotated source CIFs.

`rfd4_cache_info` reports the actual loading route, including `lightweight_reads`,
and its elapsed time inside PyMOL. This excludes desktop GUI startup/rendering;
filesystem and CPU load can affect timings.

If the old session is already blank, **start a fresh PyMOL session** with the
terminal command above. The old object and its NaN camera bounds are not repaired
in place, and ordinary loading into a populated session deliberately preserves
the current camera. Starting fresh gets both the new loader and a valid view.
Re-running `run ~/pymolrc/rfd4_viewer.py` installs the adapter in an existing
session, but does not by itself repair objects loaded before that update.

Annotation colors remain optional. With exactly one annotated object visible,
you may omit its name: `rfd4_style , hbond`. Otherwise use its actual sidebar name.
`motif` below is an example object name created explicitly with
`rfd4_load /path/input.cif, motif`, not a name automatically assigned to every file.

| Command | Display |
|---|---|
| `rfd4_style motif, coordinates` | Element colors; cyan protein C, amethyst ligand/cofactor C, orange internal C; reference atoms as lines; lighter protein C for free identities |
| `rfd4_style motif, coordinates, element_colors=0` | Flat cyan absolute-target atoms / orange internal-group atoms, including metals |
| `rfd4_style motif, coordinates, context=0` | Hide reference-only context lines |
| `rfd4_style motif, sequence` | Teal fixed identities; gray redesignable identities |
| `rfd4_style motif, unindexed` | Violet unindexed residues |
| `rfd4_style motif, hbond` | Blue donors; salmon acceptors; magenta both roles |
| `rfd4_style motif, burial` | Green buried; orange partially buried; blue exposed |
| `rfd4_style motif, elements` | Normal green protein/amethyst ligand colors; reference atoms as lines; lighter protein C for free identities |
| `rfd4_centroids motif, 0` / `rfd4_centroids motif, 1` | Hide/show the source annotation's chain-centroid markers |
| `rfd4_info motif` | Source path, condition counts, expandable length ranges |
| `rfd4_cache_info` | Session cache hits, native parses, and last reader timing |
| `help_rfd4` | Command reference |

Turn labels off with `rfd4_style motif, coordinates, labels=0`.
`context=1` (default) shows reference-only atoms as element-colored lines in either
`coordinates` or `elements` mode; unbonded atoms, including reference metals, are
crosses. `context=0` hides their representations and labels. A residue receives
the lighter carbon shade only when none of its retained atoms has fixed sequence
identity; mixed masks retain the primary shade and remain inspectable through the
atom-level selections. `element_colors=0` disables the shade to retain exact flat
coordinate-condition colors. The other annotation modes keep their existing flat colors
and display the retained physical structure for inspection.
The default view distinguishes **retained coordinates**, **absolute coordinate
targets**, **internal group geometry**, and **fixed sequence identity**. A retained
backbone is not necessarily constrained. Internal groups preserve relative geometry
without fixing the group's placement; sequence conditioning never implies backbone
conditioning. These are input conditions, not guarantees of exact generated output.

When an annotation view is requested, the `<object>_annotations` group contains named selections such as `motif_absolute`,
`motif_internal`, `motif_reference`, `motif_fixed` (an alias for absolute),
`motif_sequence_free`, `motif_hb_donor`, and `motif_partial`. Other objects keep
their current styling when a new annotated CIF is loaded. An existing object
name is never overwritten or extended: reloading gets a new unique name.
Browsing preserves condition colors. Saved PyMOL sessions preserve their own
colors, representations, visibility, and camera when restored; navigation also
preserves colors explicitly assigned by review scripts using `load_raw`.

Restart PyMOL to activate the startup-hook changes. Re-running the viewer alone
does not update the already loaded startup hooks or old object metadata; reload
the original CIF into a new object to obtain the new coordinate categories.

The source CIF is unchanged. Physical coordinates, formal charges, and native
bond topology are retained; expandable scaffold placeholders, zero-occupancy
atoms, and nonfinite coordinates are omitted from the view. Printed scaffold
length ranges are generation instructions, not extra atoms to display. The
viewer uses the asymmetric unit, first model, and first alternate location.
Residue labels are annotated-CIF IDs, not the original protein numbering or eventual designed sequence positions.
Lines show reference geometry without constrained coordinates, even within a residue
whose amino-acid identity is fixed. Carbon shade describes identity, not coordinate freedom.

H-bond roles and burial colors show **requested conditioning**, not measured
contacts or solvent accessibility. The viewer supports the categories listed
above; other native annotations remain in the source file. It does not edit
conditions. Save screenshots or a PyMOL session for viewing, but keep the original
annotated CIF for RFD4: exporting a structure from PyMOL can include the visual
centroid markers and does not preserve the RFD4 annotation tables. Such exports
are not replacements for the original input. Named selections and styling survive a `.pse`; the
`rfd4_style`/`rfd4_info` Python metadata is restored by loading the source again.

Keep `rfd4_viewer.py` together with `scripts/rfd4_cif_export.py`,
`scripts/rfd4_model_cache.py`, `scripts/rfd4_cif_detect.py`, and
`scripts/rfd4_cif_fast.py`; copying the whole repository preserves this layout.
The lightweight reader requires NumPy and Biopython in **PyMOL's Python
environment**. Install them with that environment's package manager, for example
`python -m pip install numpy biopython` using its Python executable. Missing
lightweight dependencies select the native fallback.

The fallback requires an installed RFD4 environment with its native parsing
dependencies. Its defaults are `~/git/RFD4-Proteina-dev/.venv/bin/python` and
cluster `digs`.
For another installation, pass `repo=/path/to/RFD4, python=/path/to/python` to
`rfd4_load`, or set `RFD4_REPO`, `RFD4_PYTHON`, and `RFD4_CLUSTER` before starting
PyMOL. A valid native cached load checks installation provenance but does not
launch that environment. Supported lightweight loads need no RFD4 installation.
