# WORK IN PROGRESS!

cmake_minimum_required(VERSION 2.8)
project(MALIGNER)

# compiler
SET(CMAKE_CXX_COMPILER "g++")
set(CMAKE_CXX_FLAGS "-std=c++11")


#########################################
# LMM 12/17: Commenting out Boost and Python stuff for now.

# Boost
# set(Boost_USE_STATIC_LIBS   ON)
# set(Boost_USE_MULTITHREADED ON)
# find_package(Boost COMPONENTS
#                 system
#                 thread
#                 python
#              REQUIRED)
# include_directories(${Boost_INCLUDE_DIRS})
# link_directories(${Boost_LIBRARY_DIRS})

# # Python
# find_package(PythonLibs REQUIRED)
# include_directories(${PYTHON_INCLUDE_DIRS})

# add_library(test SHARED
#         main.cpp
#     )
# # need to use .so suffix instead of .dylib, as python only looks for .so
# set_target_properties( test
#     PROPERTIES
#     IMPORT_SUFFIX ".so"
# )
# target_link_libraries( test
#         ${Boost_LIBRARIES}
#         ${PYTHON_LIBRARIES}
#     )

# various project names
# go on adding multiple names for different projects
SET(PROJECT_NAMES kmer_match)

foreach(PROJECT_NAME ${PROJECT_NAMES})

  # add sub-directories
  add_subdirectory("${CMAKE_SOURCE_DIR}/${PROJECT_NAME}"
    "${CMAKE_CURRENT_BINARY_DIR}/out/${PROJECT_NAME}")

endforeach(PROJECT_NAME ${PROJECT_NAMES})

