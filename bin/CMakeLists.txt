# Copy files from here to the binary directory.
set( FILES_TO_INSTALL
  convert_maps filter_alignments sample_maps_file
  split_maps_file summarize_maps_file alignment_locs alignment_chunks alignmnet_summary_stats
)

install(FILES ${FILES_TO_INSTALL}
        DESTINATION ${MALIGNER_BIN_DIR}
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
