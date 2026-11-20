# Inspect RFdiffusion4-Proteina CIFs in PyMOL

```sh
pymolrfd4 prediction.cif annotated.cif -d 'set grid_mode, 1'
```

The launcher enables annotation inspection for that PyMOL process, including later `load` and `browse_list` calls. Ordinary `pymol` takes the standard PyMOL loading path. To inspect one file in either session, enter `rfd4_load annotated.cif, motif`. To change routing for future opens, use `rfd4_autoload on` or `rfd4_autoload off`. An explicitly inspected file opens in a new object, even if you already loaded it normally.

The reader recognizes declared condition fields in `.cif`, `.mmcif`, and their gzip forms. AF3 predictions and CIFs with only scaffold bookkeeping use ordinary PyMOL loading. Some instruction CIFs contain placeholders with NaN coordinates. When native loading gives an invalid camera, the startup hook frames the finite atoms. Inspection mode also omits the placeholders and displays their segment instructions.

## Read the coordinate view

| Appearance | Meaning |
|---|---|
| Ochre cartoon, gold beads | Absolute coordinates are fixed |
| Lilac cartoon and beads | The group’s internal geometry is fixed. Its overall placement is free |
| Blue cartoon | Reference backbone with free coordinates |
| Thin lines or crosses | Reference ligand atoms and disconnected fragments |
| Labeled translucent pink sphere | A chain-centroid target supplied in the CIF |

A cartoon follows the mask on its alpha carbon. A residue can contain a mixture of fixed and free atoms. Use its beads or named selections for an atom-level readout. The stick color reports chemistry, including bright orange phosphorus, blue nitrogen, red oxygen, and yellow sulfur in the default palette. In chemical coloring modes, protein carbons are lighter when every retained atom in the residue has free sequence identity.

On proteins with at least 40 absolute-target alpha carbons, continuous, fully fixed residues use the ochre cartoon alone. Partial residues, internal groups, termini, ligands, and isolated targets keep their atom detail. Nonmetal beads default to `fixed_scale=0.16`. Metal spheres inherit PyMOL’s global `sphere_scale`. A new inspection starts with residue labels off.

```pml
rfd4_style motif, coordinates, markers=all
rfd4_style motif, coordinates, fixed_scale=0.12, fixed_color=yellow
rfd4_style motif, coordinates, fixed_cartoon=wheat, internal_color=violet
rfd4_style motif, coordinates, overview_min=100
rfd4_style motif, coordinates, context=lines
rfd4_style motif, coordinates, context=0, labels=1
```

`context=lines` (or `1`) restores atom-level reference lines. `context=0` hides reference context. `element_colors=0` applies the status color to every target atom, including heteroatoms. Use `rfd4_style motif, legacy` for the original cyan/orange coloring, soft elements, reference lines, and 0.25-scale nonmetal beads.

## Other annotations

| Command | View |
|---|---|
| `rfd4_style motif, sequence` | Teal fixed identities and gray free identities |
| `rfd4_style motif, unindexed` | Violet unindexed atoms |
| `rfd4_style motif, hbond` | Blue donors, salmon acceptors, magenta for both |
| `rfd4_style motif, burial` | Green buried, orange intermediate, blue exposed |
| `rfd4_style motif, elements` | Molecular palette with chemical element colors |
| `rfd4_centroids motif, 0` or `1` | Hide or show the supplied centroid targets |
| `rfd4_info motif` | Print source, condition counts, and expandable segment ranges |
| `rfd4_cache_info` | Print the reader route, read time, and cache hits |

You can omit the object name when exactly one annotated object is visible, for example `rfd4_style , sequence`. The annotation group starts closed and disabled so it takes no grid cell. Expand `<object>_annotations` for selections such as `_absolute`, `_fixed` (the absolute alias), `_internal`, `_reference`, and `_sequence_free`. Absolute, internal, and reference masks are available on opening. Other named selections are created on demand or by explicit `rfd4_load`.

Coordinate freedom and sequence identity are independent. H-bond and burial values describe requested conditions. Centroid positions come directly from the annotations, and a target at the origin is valid. Markers follow the parent object’s transforms. Saved PyMOL sessions preserve appearances, named selections, and the camera. Reload the annotated source to use its Python metadata in a new session.

## Reader and cache

The inspection view retains finite, occupied physical atoms, formal charges, and bond orders from the asymmetric unit, first model, and first alternate location. Expandable segments and zero-occupancy atoms are excluded. Keep the original CIF for downstream use. PyMOL exports can contain display-only centroid markers and omit annotation tables.

Self-contained single-model CIFs with matching author and label identifiers and local chemical templates use an in-process NumPy/Biopython reader. It supports probability and source-bookkeeping columns in completed outputs, locally defined single-atom ions and waters, simple peptide links, and unambiguous covalent disulfides. Ambiguous links, alternate locations, unsupported chemistry, or unfamiliar conditions use the native reader. The fast reader checks identities and condition masks before display.

The native CPU reader uses `~/git/RFD4-Proteina-dev/.venv/bin/python` and cluster `digs` by default. Set `RFD4_REPO`, `RFD4_PYTHON`, or `RFD4_CLUSTER` to force that runtime, or pass `repo=` and `python=` to `rfd4_load`. Verified native display models are cached by source content and reader fingerprint in `~/.cache/pymol/rfd4`. Set `RFD4_PYMOL_CACHE` to change the location. A process-local cache holds up to 32 native models. The native worker exits after a minute of inactivity.

To warm the native display cache for a directory, run the following command in a suitable CPU allocation.

```sh
python scripts/rfd4_model_cache.py --source-dir /path/to/annotated --jobs 4
```

To measure file opening on your installation, run the benchmark once with each launcher. It reports time from `cmd.load` through styling, excluding process startup and rendering.

```sh
RFD4_TIMING_SOURCE=/path/to/annotated.cif pymol -cq -r tests/benchmark_rfd4_loading.py
RFD4_TIMING_SOURCE=/path/to/annotated.cif pymolrfd4 -cq -r tests/benchmark_rfd4_loading.py
```
