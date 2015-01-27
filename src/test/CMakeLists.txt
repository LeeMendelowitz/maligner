# #include headers
include_directories(
  ${MALIGNER_SOURCE_DIR}/src/kmer_match
  ${MALIGNER_SOURCE_DIR}/src/dp
  ${MALIGNER_SOURCE_DIR}/src/common
  ${MALIGNER_SOURCE_DIR}/src/bin
)


set(test_scorematrix_src "test_scorematrix.cpp")
add_executable(test_scorematrix ${test_scorematrix_src})
target_link_libraries(test_scorematrix dp common)

set(test_scorematrix_ijkl_src "test_scorematrix_ijkl.cpp")
add_executable(test_scorematrix_ijkl ${test_scorematrix_ijkl_src})
target_link_libraries(test_scorematrix_ijkl dp common)

set(test_scorematrix_shape_src "test_scorematrix_shape.cpp")
add_executable(test_scorematrix_shape ${test_scorematrix_shape_src})
target_link_libraries(test_scorematrix_shape dp common)


add_executable(test_maligner_dp "test_maligner_dp.cpp")
target_link_libraries(test_maligner_dp dp kmer_match common)

add_executable(test_maligner_dp_row_vs_column "test_maligner_dp_row_vs_column.cpp")
target_link_libraries(test_maligner_dp_row_vs_column dp kmer_match common)

add_executable(test_maligner_dp_sizing_penalty "test_maligner_dp_sizing_penalty.cpp")
target_link_libraries(test_maligner_dp_sizing_penalty dp kmer_match common)

add_executable(test_maligner_dp_max_miss "test_maligner_dp_max_miss.cpp")
target_link_libraries(test_maligner_dp_max_miss dp kmer_match common)

add_executable(test_maligner_dp_get_all "test_maligner_dp_get_all.cpp")
target_link_libraries(test_maligner_dp_get_all dp kmer_match common)

add_executable(test_PartialSumsC "test_PartialSumsC.cpp")
target_link_libraries(test_PartialSumsC dp kmer_match common)

# install directory
install(TARGETS
  test_scorematrix
  test_scorematrix_ijkl
  test_scorematrix_shape
  test_maligner_dp
  test_maligner_dp_get_all
  test_maligner_dp_row_vs_column
  test_maligner_dp_sizing_penalty
  test_maligner_dp_max_miss
  test_PartialSumsC
  DESTINATION "${MALIGNER_BIN_DIR}")
