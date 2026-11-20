# Molecular palettes

Type `theme next` or `theme prev` in PyMOL to cycle through the palettes. `theme list` marks the active one. Selecting a palette changes the open structures and sets the colors for subsequent loads in the session. All palettes share the original matte lighting, orthographic camera, and black outlines on ray renders. `bg_white` and `bg_black` switch the background.

![Actual black viewport and white ray renders of Mineral and Fern, including an N-to-C gradient and a close-up of the same ligand and zinc ions](molecular-palettes.png)

| Theme | Sidechain carbon | Backbone carbon | Ligand carbon |
|---|---|---|---|
| **`mineral` (default)** | Taupe `#776D5F` | Warm stone `#8C8474` | Lagoon `#36A2AF` |
| `fern` | Sage `#608169` | Fern `#587A5D` | Orchid `#BA79B6` |
| `classic` | Sage `sc05` | Forest `bb05` | Amethyst `lig_amethyst3` |
| `fjord` | Mint `#8CCBC3` | Teal `#559EAA` | Lilac `#BE97DD` |
| `clay` | Peach `#E4BAA4` | Terracotta `#BD8579` | Aqua `#7BBCC5` |
| `iris` | Lavender `#B8ACE0` | Violet `#9283C2` | Jade `#7BC7A5` |
| `foliage` | Sage `sc05` | Forest `bb05` | Amethyst `lig_amethyst3` |

The molecular themes use stock element colors, including white hydrogen, blue nitrogen, red oxygen, bright orange phosphorus (`#FF8000`), and yellow sulfur. Foliage restores the historical soft elements and uniform gray-blue metals. The old custom nitrogen is also available as `legacy_nitrogen`. In Mineral and Fern, REMARK 666 residues use lighter protein-colored sticks so that phosphorus is easy to pick out.

**Ligands.** A single carbon color makes the same molecule easy to follow across models. Use `ligand_colors by_type` to assign accents by chemical component name (the PyMOL `resn` field). Copies of a component get the same accent across chains, objects, and later loads in that session. Two molecules both named `LIG` share an accent. `ligand_colors list` prints the assignments. Mineral cycles lagoon, orchid, and olive. With more components the accents repeat. `ligand_colors uniform` restores the theme’s ligand carbon. Water, metals, and atoms belonging to polymer residues keep their own colors.

**Gradients.** `theme mineral, gradient=on` uses warm umber `#786F62` at N, stone `#8C8474` at the midpoint, and pale taupe `#A09886` at C. Fern runs from `#456950` through `#587A5D` to `#77977A`. Each chain starts its own ramp in residue-number and insertion-code order. Only the cartoon changes. Molecular sticks keep their atom colors, and RFD4 condition cartoons keep their constraint colors. Use `rfd4_style , elements` to see the molecular palette on an inspected file. `theme mineral, gradient=off` returns to a flat cartoon.

Set `PYMOLRC_THEME` to a palette name and `PYMOLRC_GRADIENT` to `on` or `off` to choose your startup look. For example, `PYMOLRC_THEME=fern PYMOLRC_GRADIENT=on pymol model.cif` opens a single session with a green gradient. Export the same variables in your shell profile to make that choice persistent.

The pocket comparison used the same public model, camera, and matte lighting for each candidate. In the first round, green ligand carbons approached gray sidechains under a simulated red–green color-vision deficiency. Moving the ligand toward cyan and warming the protein improved their separation. The zinc spheres retain PyMOL’s native lavender. [Color-vision simulation method](https://pubmed.ncbi.nlm.nih.gov/19834201/) · [Original palette gallery](viewer-guide.png)

The named `sc01`–`sc10`, `bb01`–`bb10`, amethyst, Seaborn, and RFdiffusion palettes are available to existing scripts. `color_palette` applies the original per-object foliage cycle, alignment, and ligand shell. The three `color_bb_rfdiffusion*` commands apply residue gradients. [Gradient examples](COMMANDS.md) · [Original foliage render](foliage-original.png)

## Brand colors

Use these names with PyMOL’s `color` command, for example `color ipd_logo_blue, chain A and elem C`.

**[NVIDIA](https://www.nvidia.com/content/dam/en-zz/Solutions/about-us/documents/NVIDIA-Trademark-and-Logo-Usage-Guidelines.pdf)**

| Color name | HEX |
|---|---|
| `nvidia_green` | `#76B900` |
| `nvidia_black` | `#000000` |
| `nvidia_dark_gray` | `#1E1E1E` |
| `nvidia_white` | `#FFFFFF` |

Black is also specified in the [NVIDIA logo guidelines](https://www.nvidia.com/content/dam/en-zz/Solutions/about-us/NVIDIA-Brand-Guidelines-for-NVIDIA-Partner-Network-v03-3-5-19.pdf).

**[IPD](https://www.ipd.uw.edu/wp-content/uploads/2025/03/IPD-logo-UW-color-2025.png)**

These values are sampled from solid fills in IPD’s current official logo,
not a published color specification. Existing `ipd_blue` (`#4196B5`) and
`ipd_purple` (`#B41F85`) remain available with their original values.

| Color name | HEX |
|---|---|
| `ipd_logo_blue` | `#4D97B7` |
| `ipd_logo_magenta` | `#B11186` |
| `ipd_uw_purple` | `#4B2E83` |

**[UW](https://www.washington.edu/brand/brand-elements/colors/)**

Uses the published digital HEX values; some RGB labels on UW’s page differ.
`uw_purple` / `uw_gold` are aliases for Husky Purple / Husky Gold.

| Color name | HEX |
|---|---|
| `uw_purple` | `#32006E` |
| `uw_husky_purple` | `#32006E` |
| `uw_spirit_purple` | `#4B2E83` |
| `uw_gold` | `#B7A57A` |
| `uw_husky_gold` | `#B7A57A` |
| `uw_web_gold` | `#E8E3D3` |
| `uw_heritage_gold` | `#85754D` |
| `uw_spirit_gold` | `#FFC700` |
| `uw_accent_green` | `#AADB1E` |
| `uw_accent_teal` | `#2AD2C9` |
| `uw_accent_pink` | `#E93CAC` |
| `uw_accent_lavender` | `#C5B4E3` |

**[Michigan](https://brand.umich.edu/design-resources/colors/)**

| Color name | HEX |
|---|---|
| `umich_maize` | `#FFCB05` |
| `umich_blue` | `#00274C` |
| `umich_tappan_red` | `#9A3324` |
| `umich_ross_orange` | `#D86018` |
| `umich_rackham_green` | `#75988D` |
| `umich_wave_field_green` | `#A5A508` |
| `umich_taubman_teal` | `#00B2A9` |
| `umich_arboretum_blue` | `#2F65A7` |
| `umich_a2_amethyst` | `#702082` |
| `umich_matthaei_violet` | `#575294` |
| `umich_peony_pink` | `#E01F7C` |
| `umich_umma_tan` | `#CFC096` |
| `umich_burton_tower_beige` | `#9B9A6D` |
| `umich_angell_hall_ash` | `#989C97` |
| `umich_law_quad_stone` | `#655A52` |
| `umich_puma_black` | `#131516` |

Sources checked September 10, 2026.
