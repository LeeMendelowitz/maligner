# #include headers from this folder

include_directories(
  ${MALIGNER_SOURCE_DIR}/src/kmer_match
  ${MALIGNER_SOURCE_DIR}/src/dp
  ${MALIGNER_SOURCE_DIR}/src/common
)

# source files
message(STATUS "making maligner.cpp....")
message(STATUS "CMAKE_CURRENT_SOURCE_DIR : ${CMAKE_CURRENT_SOURCE_DIR}")
message(STATUS "MALIGNER_SOURCE_DIR : ${MALIGNER_SOURCE_DIR}")
message(STATUS "INCLUDE_DIRECTORIES : ${INCLUDE_DIRECTORIES}")

# build maligner
set(maligner_SRCS "maligner.cpp")
add_executable(maligner ${maligner_SRCS})
target_link_libraries(maligner kmer_match common)

# build maligner_dp
set(maligner_dp_SRCS "maligner_dp.cpp")
add_executable(maligner_dp ${maligner_dp_SRCS})
target_link_libraries(maligner_dp kmer_match dp common)

# install directory
install(TARGETS maligner maligner_dp DESTINATION "${MALIGNER_BIN_DIR}")
