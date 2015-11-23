# #include headers from this folder

include_directories(
  ${MALIGNER_SOURCE_DIR}/src/ix
  ${MALIGNER_SOURCE_DIR}/src/dp
  ${MALIGNER_SOURCE_DIR}/src/common
  ${MALIGNER_SOURCE_DIR}/src/vd
)

# source files
message(STATUS "making maligner.cpp....")
message(STATUS "CMAKE_CURRENT_SOURCE_DIR : ${CMAKE_CURRENT_SOURCE_DIR}")
message(STATUS "MALIGNER_SOURCE_DIR : ${MALIGNER_SOURCE_DIR}")
message(STATUS "INCLUDE_DIRECTORIES : ${INCLUDE_DIRECTORIES}")

# build maligner
set(maligner_ix_SRCS "maligner_ix.cpp")
add_executable(maligner_ix ${maligner_ix_SRCS})
target_link_libraries(maligner_ix ix dp common)

# build maligner_dp
set(maligner_dp_SRCS "maligner_dp.cpp")
add_executable(maligner_dp ${maligner_dp_SRCS})
target_link_libraries(maligner_dp vd ix dp common)

# build maligner_vd
set(maligner_vd_SRCS "maligner_vd.cpp")
add_executable(maligner_vd ${maligner_vd_SRCS})
target_link_libraries(maligner_vd vd ix dp common)


# install directory
install(TARGETS maligner_ix maligner_dp maligner_vd DESTINATION "${MALIGNER_BIN_DIR}")
