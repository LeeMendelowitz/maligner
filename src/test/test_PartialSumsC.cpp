#include "align.h"

#include <vector>
#include <iostream>
#include <cstdlib>

using maligner_dp::PartialSums;

int main() {
  
  std::vector<int> iv {1, 2, 3, 4, 5};
  size_t max_misses = 2;

  PartialSums ps_forward(iv, max_misses, PartialSums::forward_tag());
  for (size_t i = 0; i < iv.size(); i++) {
    for (size_t j = 0; j < max_misses + 1; j++) {
      std::cout << "i: " << i << " j: " << j << " " << ps_forward(i,j) << "\n";
    }
  }

  PartialSums ps_rev(iv, max_misses, PartialSums::reverse_tag());
  for (size_t i = 0; i < iv.size(); i++) {
    for (size_t j = 0; j < max_misses + 1; j++) {
      std::cout << "i: " << i << " j: " << j << " " << ps_rev(i,j) << "\n";
    }
  }

  return EXIT_SUCCESS;
}