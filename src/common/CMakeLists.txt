# #include headers from this folder
# include_directories(${UTIL_SOURCE_DIR}/include)

# source files
set(COMMON_SRC
  "timer.cpp"
  "map.cpp"
  "map_reader.cpp"
  "common_math.cpp"
  )

# Build Library
add_library(common STATIC ${COMMON_SRC})

# install directory
# install(TARGETS common ARCHIVE DESTINATION "${MALIGNER_LIB_DIR}")