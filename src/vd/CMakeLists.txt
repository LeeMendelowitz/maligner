# #include headers from this folder
# include_directories(${UTIL_SOURCE_DIR}/include)

include_directories(
  ${MALIGNER_SOURCE_DIR}/src/common
  ${MALIGNER_SOURCE_DIR}/src/dp
)

# source files
set(VD_SOURCES
  "score_matrix_vd.cpp"
  "score_matrix_profile.cpp"
)

# Build Library
add_library(vd STATIC ${VD_SOURCES})

# install directory
# install(TARGETS dp ARCHIVE DESTINATION "${MALIGNER_LIB_DIR}")
