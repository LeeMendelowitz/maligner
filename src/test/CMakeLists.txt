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

add_executable(test_maligner_dp "test_maligner_dp.cpp")
target_link_libraries(test_maligner_dp dp kmer_match common)

add_executable(test_maligner_dp_get_all "test_maligner_dp_get_all.cpp")
target_link_libraries(test_maligner_dp_get_all dp kmer_match common)

# install directory
install(TARGETS test_scorematrix test_scorematrix_ijkl test_maligner_dp test_maligner_dp_get_all DESTINATION "${MALIGNER_BIN_DIR}")
