def pymolrc_help(_self=cmd):
    """Top-level index of every command group in this pymolrc config."""
    print("""
=== pymolrc command groups ===  (repo: github.com/SethWoodbury/pymolrc)
Type the group-help command for full per-command args (mandatory vs optional):

  Crystallography (crystals + maps + alignment) ....... help_crystal
  RFdiffusion3 movies / figures ....................... help_rfd3
  RFdiffusion4 annotations (pymolrfd4) ................ help_rfd4

Crystallography : xtal  xtal_align  xtal_density  xtal_rms  xtal_catres  xtal_focus
RFdiffusion3    : rfd3_movie  rfd3_check_frames  style_fixed  apply_camera
                  draw_connectors  catalytic_sel_from_motif  add_custom_bond  remove_bond
Look & color    : theme  ligand_colors  color_palette  autostyle  style_all  publication_ray_trace
                  color_bb_rfdiffusion  color_bb_rfdiffusion3
                  color_bb_rfdiffusion4_proteina
                  gaussian_mode  gaussian_off  gaussian_spin  bg_white  bg_black
Align & arrange : align_all  center_all  autoalign  autosolo
Select & inspect: seq  show_catres  only_catres  show_metals
Browse many     : browse_list  browse_next  browse_prev  browse_step  browse_prestyle
                  (arrow-key structure nav: structure_next/prev/step/show_all)
Bond fixes      : add_custom_bond  remove_bond  clash_bond_filter  ligand_center_filter

Full feature list + install/revert notes: README.md
""")


cmd.extend("pymolrc_help", pymolrc_help)
cmd.extend("help_pymolrc", pymolrc_help)
