# #include headers from this folder
# include_directories(${UTIL_SOURCE_DIR}/include)

# source files
set(KMER_MATCH_SOURCES "error_model.cpp" "map.cpp" "map_chunk.cpp" "map_chunk_db.cpp" "map_frag.cpp" "ref_alignment.cpp")

# Build Library
add_library(kmer_match STATIC ${KMER_MATCH_SOURCES})

# install directory
install(TARGETS kmer_match 
  ARCHIVE DESTINATION "${MALIGNER_SOURCE_DIR}/lib" )