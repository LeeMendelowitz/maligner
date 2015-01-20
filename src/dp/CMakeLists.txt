# #include headers from this folder
# include_directories(${UTIL_SOURCE_DIR}/include)

# source files
set(DP_SOURCES
  "align.cc"
  "globals.cc"
  "ScoreCell.cc"
  # "ScoreMatrix.cc"
  "alignment.cc"
  "chunk.cc"
  "matched_chunk.cc"
  "map_data.cc"
)

# Build Library
add_library(dp STATIC ${DP_SOURCES})

# install directory
install(TARGETS dp ARCHIVE DESTINATION "${MALIGNER_LIB_DIR}")

# CAN ATTEMPT TO BUILD THE PYTHON LIBRARY WITH HELP FROM THIS:
# http://stackoverflow.com/questions/5851016/boost-python-and-cmake-link-and-load-errors
# Best to move the malignpy glue code to a difference directory, outside of dp.
# Also boost python may have trouble using the std::shared_ptr. Please test tomorrow.