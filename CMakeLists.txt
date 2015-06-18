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
 

set(MALIGNER_LIB_DIR ${MALIGNER_BINARY_DIR}/lib)
set(MALIGNER_BIN_DIR ${MALIGNER_BINARY_DIR}/bin)
message(STATUS "MALIGNER_LIB_DIR : ${MALIGNER_LIB_DIR}")
message(STATUS "MALIGNER_BIN_DIR : ${MALIGNER_BIN_DIR}")
message(STATUS "CMAKE_ROOT : ${CMAKE_ROOT}")

include_directories(${MALIGNER_BINARY_DIR}/src/dp ${MALIGNER_BINARY_DIR}/src/kmer_match)
add_subdirectory(src)
add_subdirectory(bin)
add_subdirectory(lib)