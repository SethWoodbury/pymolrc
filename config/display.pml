# Environment
cmd.space('rgb')

# Original matte cartoon lighting and black outlines on ray-traced images.
set specular, 0
unset shininess
set ray_shadow, off
set ray_trace_mode, 1
set ray_trace_disco_factor, 1
set ray_trace_gain, 0.1
set ray_trace_color, black
set ray_texture, 0
set antialias, 2
set power, 0.2
set ambient, 0.4
set direct, 0.45
set reflect, 1
set reflect_power, 0
set cartoon_sampling, 10

# Structure display
set orthoscopic, 1
set valence, off
set cartoon_flat_sheets, off
set cartoon_fancy_sheets, 1
set cartoon_fancy_helices, 1
set cartoon_transparency, 0.0
set cartoon_gap_cutoff, 0
set ribbon_width, 8
set line_width, 3

# Performance (faster building/loading, no visual change)
set async_builds, 1
# Global stick/sphere sizing so the per-object auto-styling never calls the costly
# per-selection `set` (that was the main load bottleneck). stick_radius 0.25 is the
# PyMOL default, sphere_scale 0.7 is the metal-sphere size. Override per-object anytime.
set stick_radius, 0.25
set sphere_scale, 0.7

# --- Foliage-forward protein carbon palettes (cycling per object) ---
set_color sc01 = [165,206,154]
set_color sc02 = [155,198,147]
set_color sc03 = [146,191,140]
set_color sc04 = [137,184,133]
set_color sc05 = [128,177,126]
set_color sc06 = [120,170,120]
set_color sc07 = [112,163,114]
set_color sc08 = [104,156,108]
set_color sc09 = [97,149,102]
set_color sc10 = [90,142,97]

set_color bb01 = [120,172,115]
set_color bb02 = [108,164,109]
set_color bb03 = [97,155,99]
set_color bb04 = [89,146,94]
set_color bb05 = [81,138,89]
set_color bb06 = [74,130,84]
set_color bb07 = [67,122,80]
set_color bb08 = [61,114,75]
set_color bb09 = [55,106,71]
set_color bb10 = [50,99,67]

# --- Ligand carbon options (refined amethyst family) ---
set_color lig_amethyst1 = [201,170,221]
set_color lig_amethyst2 = [185,151,208]
set_color lig_amethyst3 = [168,132,194]
set_color lig_amethyst4 = [152,114,180]
set_color lig_amethyst5 = [136,97,166]

# --- Element colors (soft-but-distinct) ---
set_color elemN = [67,120,209]
set_color elemO = [214,86,98]
set_color elemS = [212,177,89]
set_color elemP = [188,128,78]
set_color elemH = [237,239,234]
set_color elemZn = [109,124,142]

# --- Legacy colors (kept for manual use) ---
set_color legacy_nitrogen = [2,118,253]
set_color pymol_gray, [51,51,51]
set_color pymol_black, [34,34,34]
set_color good_yellow, [250,199,44]
set_color good_teal, [40,176,193]
set_color good_green, [170,195,47]
set_color good_pink, [236,114,164]
set_color good_peach, [249,145,120]
set_color good_blue, [68,153,231]
set_color good_gray, [220,220,220]
set_color good_red, [228,74,62]
set_color ipd_blue, [65,150,181]
set_color ipd_purple, [180,31,133]

# --- Brand colors (manual use; digital HEX values, checked 2026-09-10) ---

# NVIDIA
# Sources: https://www.nvidia.com/content/dam/en-zz/Solutions/about-us/documents/NVIDIA-Trademark-and-Logo-Usage-Guidelines.pdf
# Black: https://www.nvidia.com/content/dam/en-zz/Solutions/about-us/NVIDIA-Brand-Guidelines-for-NVIDIA-Partner-Network-v03-3-5-19.pdf
set_color nvidia_green, [118,185,0]  # 76B900
set_color nvidia_black, [0,0,0]  # 000000
set_color nvidia_dark_gray, [30,30,30]  # 1E1E1E
set_color nvidia_white, [255,255,255]  # FFFFFF

# Institute for Protein Design
# RGB sampled from the solid fills of the current official IPD logo:
# https://www.ipd.uw.edu/wp-content/uploads/2025/03/IPD-logo-UW-color-2025.png
# Existing ipd_blue / ipd_purple above retain their legacy values.
set_color ipd_logo_blue, [77,151,183]  # 4D97B7
set_color ipd_logo_magenta, [177,17,134]  # B11186
set_color ipd_uw_purple, [75,46,131]  # 4B2E83

# University of Washington
# https://www.washington.edu/brand/brand-elements/colors/
# Use published HEX values (some RGB labels on the page disagree).
# uw_purple / uw_gold alias Husky Purple / Husky Gold.
set_color uw_purple, [50,0,110]  # 32006E
set_color uw_husky_purple, [50,0,110]  # 32006E
set_color uw_spirit_purple, [75,46,131]  # 4B2E83
set_color uw_gold, [183,165,122]  # B7A57A
set_color uw_husky_gold, [183,165,122]  # B7A57A
set_color uw_web_gold, [232,227,211]  # E8E3D3
set_color uw_heritage_gold, [133,117,77]  # 85754D
set_color uw_spirit_gold, [255,199,0]  # FFC700
set_color uw_accent_green, [170,219,30]  # AADB1E
set_color uw_accent_teal, [42,210,201]  # 2AD2C9
set_color uw_accent_pink, [233,60,172]  # E93CAC
set_color uw_accent_lavender, [197,180,227]  # C5B4E3

# University of Michigan
# https://brand.umich.edu/design-resources/colors/
set_color umich_maize, [255,203,5]  # FFCB05
set_color umich_blue, [0,39,76]  # 00274C
set_color umich_tappan_red, [154,51,36]  # 9A3324
set_color umich_ross_orange, [216,96,24]  # D86018
set_color umich_rackham_green, [117,152,141]  # 75988D
set_color umich_wave_field_green, [165,165,8]  # A5A508
set_color umich_taubman_teal, [0,178,169]  # 00B2A9
set_color umich_arboretum_blue, [47,101,167]  # 2F65A7
set_color umich_a2_amethyst, [112,32,130]  # 702082
set_color umich_matthaei_violet, [87,82,148]  # 575294
set_color umich_peony_pink, [224,31,124]  # E01F7C
set_color umich_umma_tan, [207,192,150]  # CFC096
set_color umich_burton_tower_beige, [155,154,109]  # 9B9A6D
set_color umich_angell_hall_ash, [152,156,151]  # 989C97
set_color umich_law_quad_stone, [101,90,82]  # 655A52
set_color umich_puma_black, [19,21,22]  # 131516

# --- Seaborn palette ---
set_color sblightgreen, [153,204,153]
set_color sblightpurple, [204,153,204]
set_color sblightorange, [255,219,153]
set_color sblightblue, [153,153,255]
set_color sblightred, [255,153,153]
set_color sblightyellow, [255,255,153]

# --- RFdiffusion gradient ---
# RFd_darkblue, RFd_purple, RFd_pink, and paper_teal are defined by the RFD3 kit.
set_color RFd_blue, [102,134,197]
set_color RFd_lightblue, [149,150,198]
set_color RFd_melon, [255,198,178]
set_color RFd_navaho, [255,224,172]

# --- RFD4-Proteina: Northern Nights molecular_spectrum, N->C ---
# paper/style/theme.json at f336faf492d3, the revision cited by
# https://rfd4-proteina.org/assets/site-colors.css (checked 2026-11-02).
set_color rfd4p_teal, [24,108,112]  # 186C70
set_color rfd4p_blue, [59,90,141]  # 3B5A8D
set_color rfd4p_violet, [70,61,115]  # 463D73
set_color rfd4p_purple, [141,82,173]  # 8D52AD

# Molecular themes. Hex values are listed beside their names in docs/COLORS.md.
set_color mineral_stone, [140,132,116]
set_color mineral_taupe, [119,109,95]
set_color mineral_site, [159,145,119]
set_color mineral_lagoon, [54,162,175]
set_color fern_leaf, [88,122,93]
set_color fern_sage, [96,129,105]
set_color fern_site, [119,153,122]
set_color fern_orchid, [186,121,182]
set_color ligand_olive, [139,151,80]
set_color fjord_mint, [140,203,195]
set_color fjord_teal, [85,158,170]
set_color fjord_lilac, [190,151,221]
set_color clay_peach, [228,186,164]
set_color clay_terracotta, [189,133,121]
set_color clay_aqua, [123,188,197]
set_color iris_lavender, [184,172,224]
set_color iris_violet, [146,131,194]
set_color iris_jade, [123,199,165]

# Background


bg_color black
