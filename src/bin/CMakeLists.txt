# #include headers from this folder

include_directories(${MALIGNER_SOURCE_DIR}/src/kmer_match)

# source files
message(STATUS "making maligner.cpp....")
message(STATUS "CMAKE_CURRENT_SOURCE_DIR : ${CMAKE_CURRENT_SOURCE_DIR}")
message(STATUS "MALIGNER_SOURCE_DIR : ${MALIGNER_SOURCE_DIR}")
message(STATUS "INCLUDE_DIRECTORIES : ${INCLUDE_DIRECTORIES}")

# Build Library
set(maligner_SRCS "maligner.cpp")
add_executable(maligner ${maligner_SRCS})
target_link_libraries(maligner kmer_match)

# install directory
install(TARGETS maligner DESTINATION "${MALIGNER_BIN_DIR}")