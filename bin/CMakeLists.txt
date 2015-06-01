# Copy files from here to the binary directory.
set( FILES_TO_INSTALL
  alignment_chunks
  alignment_locs
  alignment_summary_stats
  alignments_maligner_to_malign_viz  
  alignments_maligner_to_soma
  alignments_soma_to_maligner
  convert_maps
  filter_alignments
  get_maps_by_id
  get_boundary_frags
  sample_alignments
  sample_maps_file
  soma_consensus_maps
  soma_query_maps
  split_maps_file
  summarize_maps_file
  smooth_maps_file
  make_insilico_map
  generate_random_maps
  parse_nucmer_coords
  parse_somacbcb_locs
)

install(FILES ${FILES_TO_INSTALL}
        DESTINATION ${MALIGNER_BIN_DIR}
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
  