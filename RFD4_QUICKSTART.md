# RFD4 annotated CIFs in PyMOL

Open an annotated input exactly like a normal structure:

```bash
pymol /absolute/path/to/annotated_input.cif
```

Or use `load /absolute/path/to/annotated_input.cif` inside PyMOL. The installed
`~/.pymolrc` automatically detects native RFD4 condition fields and displays the
physical structure as sticks with your usual element colors. No special loader
command is needed. Ordinary PDB/CIF files retain their existing loader behavior.
The same automatic handling works with `browse_list` and `random_shuffle`.

The native input contains an expandable scaffold placeholder whose coordinates
are NaN. Ordinary PyMOL's mmCIF loader includes it, making the viewing bounds
NaN and the scene blank. The automatic adapter instead uses the native RFD4/
AtomWorks parser and omits instruction placeholders from the display; it never
changes the annotated source. Native coordinates, charges and bonds are retained.

Annotated files now use a persistent, verified JSON display cache under
`~/.cache/pymol/rfd4/`. A cached file opens directly in PyMOL, including in a
fresh session, without starting the native scientific Python environment.
The source bytes, reader/condition code, parsing settings, package versions,
and local configuration determine cache validity. Changed input or reader code
causes a fresh native parse. Cached coordinates, bonds, charges, and annotation
records are checksummed and validated before use; corrupt entries are rebuilt.
Code hashes are reused only when their file identity/size/ctime/mtime match.
Native package binary/data files are tracked by their filesystem change stamps.
File renames reuse the same content cache and display the current source path.

On an uncached file, the native CPU parser starts once, saves the resulting
display data, and is reused for subsequent uncached files. It exits after 60
seconds without requests. A 32-model in-memory cache also speeds backtracking.
No checkpoint, remote server, or GPU is used. Re-running the helper preserves
loaded annotation metadata and does not stack duplicate loader hooks.

For a new large input library, optionally prepare its display cache in advance
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

The measured fresh-session load for the 216-atom example dropped from about
9.7 seconds to 0.043 seconds after preparing its cache. These are model-loading
measurements, excluding desktop GUI startup and network filesystem variability.

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
| `rfd4_style motif, coordinates` | Cyan coordinate targets; gray reference coordinates |
| `rfd4_style motif, sequence` | Teal fixed identities; gray redesignable identities |
| `rfd4_style motif, unindexed` | Violet unindexed residues |
| `rfd4_style motif, hbond` | Blue donors; salmon acceptors; magenta both roles |
| `rfd4_style motif, burial` | Green buried; orange partially buried; blue exposed |
| `rfd4_style motif, elements` | Conventional element colors |
| `rfd4_info motif` | Source path, condition counts, expandable length ranges |
| `rfd4_cache_info` | Session cache hits, native parses, and last reader timing |
| `help_rfd4` | Command reference |

Turn labels off with `rfd4_style motif, coordinates, labels=0`.
When an annotation view is requested, the `<object>_annotations` group contains named selections such as `motif_fixed`,
`motif_sequence_free`, `motif_hb_donor`, and `motif_partial`. Other objects keep
their current styling when a new annotated CIF is loaded. An existing object
name is never overwritten or extended: reloading gets a new unique name.

The source CIF is unchanged. Physical coordinates, formal charges, and native
bond topology are retained; expandable scaffold placeholders, zero-occupancy
atoms, and nonfinite coordinates are omitted from the view. Printed scaffold
length ranges are generation instructions, not extra atoms to display. The
viewer uses the asymmetric unit, first model, and first alternate location.
Residue labels are source IDs, not the eventual designed sequence positions.
Gray atoms can be present as reference geometry without fixed coordinates.

H-bond roles and burial colors show **requested conditioning**, not measured
contacts or solvent accessibility. The viewer supports the categories listed
above; other native annotations remain in the source file. It does not edit
conditions. Save screenshots or a PyMOL session for viewing, but keep the original
annotated CIF for RFD4: exporting a structure from PyMOL does not preserve the
RFD4 annotation tables. Named selections and styling survive a `.pse`; the
`rfd4_style`/`rfd4_info` Python metadata is restored by loading the source again.

Defaults: `~/git/RFD4-Proteina-dev/.venv/bin/python`, cluster `digs`.
For another installation, pass `repo=/path/to/RFD4, python=/path/to/python` to
`rfd4_load`, or set `RFD4_REPO`, `RFD4_PYTHON`, and `RFD4_CLUSTER` before starting
PyMOL. Copy this whole directory if moving machines: the helper also needs
`scripts/rfd4_cif_export.py`, `scripts/rfd4_model_cache.py`, and an installed
RFD4 environment. A valid cached load checks installation provenance but does
not launch that environment.
