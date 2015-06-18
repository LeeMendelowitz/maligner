# #include headers
include_directories(
  ${MALIGNER_SOURCE_DIR}/src/ix
  ${MALIGNER_SOURCE_DIR}/src/dp
  ${MALIGNER_SOURCE_DIR}/src/common
  ${MALIGNER_SOURCE_DIR}/src/bin
)


# set(test_scorematrix_src "test_scorematrix.cpp")
# add_executable(test_scorematrix ${test_scorematrix_src})
# target_link_libraries(test_scorematrix dp common)

# set(test_scorematrix_ijkl_src "test_scorematrix_ijkl.cpp")
# add_executable(test_scorematrix_ijkl ${test_scorematrix_ijkl_src})
# target_link_libraries(test_scorematrix_ijkl dp common)

# set(test_scorematrix_shape_src "test_scorematrix_shape.cpp")
# add_executable(test_scorematrix_shape ${test_scorematrix_shape_src})
# target_link_libraries(test_scorematrix_shape dp common)


# add_executable(test_maligner_dp "test_maligner_dp.cpp")
# target_link_libraries(test_maligner_dp dp ix common)

# add_executable(test_maligner_dp_row_vs_column "test_maligner_dp_row_vs_column.cpp")
# target_link_libraries(test_maligner_dp_row_vs_column dp ix common)

# add_executable(test_maligner_dp_sizing_penalty "test_maligner_dp_sizing_penalty.cpp")
# target_link_libraries(test_maligner_dp_sizing_penalty dp ix common)

# add_executable(test_maligner_dp_max_miss "test_maligner_dp_max_miss.cpp")
# target_link_libraries(test_maligner_dp_max_miss dp ix common)

# add_executable(test_maligner_dp_get_all "test_maligner_dp_get_all.cpp")
# target_link_libraries(test_maligner_dp_get_all dp ix common)


add_executable(test_mad "test_mad.cpp")
target_link_libraries(test_mad common)

add_executable(test_map_sampler "test_map_sampler.cpp")
target_link_libraries(test_map_sampler common)

add_executable(test_partial_sums "test_PartialSumsC.cpp")
target_link_libraries(test_partial_sums dp ix common)

# install directory
# install(TARGETS
#   test_scorematrix
#   test_scorematrix_ijkl
#   test_scorematrix_shape
#   test_maligner_dp
#   test_maligner_dp_get_all
#   test_maligner_dp_row_vs_column
#   test_maligner_dp_sizing_penalty
#   test_maligner_dp_max_miss
#   test_PartialSumsC
#   test_mad
#   DESTINATION "${MALIGNER_BIN_DIR}")
install(TARGETS
  test_map_sampler
  test_partial_sums
  DESTINATION "${MALIGNER_BIN_DIR}/test")
