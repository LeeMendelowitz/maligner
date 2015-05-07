# #include headers from this folder
# include_directories(${UTIL_SOURCE_DIR}/include)

# source files
set(COMMON_SRC
  "timer.cpp"
  "map.cpp"
  "map_data.cpp"
  "map_reader.cpp"
  "common_math.cpp"
  )

# Build Library
add_library(common STATIC ${COMMON_SRC})

# install directory
# install(TARGETS common ARCHIVE DESTINATION "${MALIGNER_LIB_DIR}")

# CAN ATTEMPT TO BUILD THE PYTHON LIBRARY WITH HELP FROM THIS:
# http://stackoverflow.com/questions/5851016/boost-python-and-cmake-link-and-load-errors
# Best to move the malignpy glue code to a difference directory, outside of dp.
# Also boost python may have trouble using the std::shared_ptr. Please test tomorrow.