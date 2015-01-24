# #include headers from this folder
# include_directories(${UTIL_SOURCE_DIR}/include)

include_directories(
  ${MALIGNER_SOURCE_DIR}/src/common
)

# source files
set(KMER_MATCH_SOURCES
  "error_model.cpp"
  "map.cpp"
  "map_reader.cpp"
  "map_chunk.cpp"
  "map_chunk_db.cpp"
  "map_frag.cpp"
  "ref_alignment.cpp"
)

# Build Library
add_library(kmer_match STATIC ${KMER_MATCH_SOURCES})

# install directory
install(TARGETS kmer_match ARCHIVE DESTINATION "${MALIGNER_LIB_DIR}")