# Copy files from here to the binary directory.
set( FILES_TO_INSTALL
  convert_maps
  filter_alignments
  sample_maps_file
  split_maps_file
  summarize_maps_file
  sample_alignments
  alignment_locs
  alignment_chunks
  alignment_summary_stats
  soma_consensus_maps
)

install(FILES ${FILES_TO_INSTALL}
        DESTINATION ${MALIGNER_BIN_DIR}
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
