cmake_minimum_required(VERSION 2.8.11)
project(MALIGNER)

# For CMAKE help, check out this resource:
# http://www.cmake.org/Wiki/CMake

message(STATUS "CMAKE_GENERATOR : ${CMAKE_GENERATOR}")
message(STATUS "CMAKE_BUILD_TYPE : ${CMAKE_BUILD_TYPE}")

# Compiler flags
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
        message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

message(STATUS "Using CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")

set(MALIGNER_LIB_DIR ${MALIGNER_BINARY_DIR}/lib)
set(MALIGNER_BIN_DIR ${MALIGNER_BINARY_DIR}/bin)
message(STATUS "MALIGNER_LIB_DIR : ${MALIGNER_LIB_DIR}")
message(STATUS "MALIGNER_BIN_DIR : ${MALIGNER_BIN_DIR}")
message(STATUS "CMAKE_ROOT : ${CMAKE_ROOT}")

include_directories(${MALIGNER_BINARY_DIR}/src/dp ${MALIGNER_BINARY_DIR}/src/kmer_match)
add_subdirectory(src)
add_subdirectory(bin)
add_subdirectory(lib)