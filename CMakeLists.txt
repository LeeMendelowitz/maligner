cmake_minimum_required(VERSION 2.8.11)
project(MALIGNER)

# For CMAKE help, check out this resource:
# http://www.cmake.org/Wiki/CMake

message(STATUS "CMAKE_GENERATOR : ${CMAKE_GENERATOR}")
message(STATUS "CMAKE_BUILD_TYPE : ${CMAKE_BUILD_TYPE}")

# Compiler flags
set(CMAKE_CXX_COMPILER "g++")
set(CMAKE_CXX_FLAGS_RELEASE "-std=c++11 -O3")
set(CMAKE_CXX_FLAGS_DEBUG "-std=c++11 -O0 -g")
# set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS_DEBUG})

#########################################
# LMM 12/17: Commenting out Boost and Python stuff for now. Investigate further.
# src from here: http://stackoverflow.com/questions/5851016/boost-python-and-cmake-link-and-load-errors
#
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

set(MALIGNER_LIB_DIR ${MALIGNER_BINARY_DIR}/lib)
set(MALIGNER_BIN_DIR ${MALIGNER_BINARY_DIR}/bin)
message(STATUS "MALIGNER_LIB_DIR : ${MALIGNER_LIB_DIR}")
message(STATUS "MALIGNER_BIN_DIR : ${MALIGNER_BIN_DIR}")

include_directories(${MALIGNER_BINARY_DIR}/src/dp ${MALIGNER_BINARY_DIR}/src/kmer_match)

add_subdirectory (src)

