#include "align.h"

#include <vector>
#include <iostream>
#include <cstdlib>

using maligner_dp::PartialSums;


template< typename T>
std::ostream& operator<<(std::ostream& os, std::vector<T> vec) {
  size_t n = vec.size();
  os << "[ ";
  for (size_t i = 0; i < n; i++) {
    os << vec[i];
    if (i != n-1) {
      os << ", ";
    }
  }
  os << " ]";
  return os;
}

int main() {
  
  std::vector<int> iv {1, 2, 3, 4, 5};
  std::vector<int> iv_rev(iv.rbegin(), iv.rend());
  size_t max_misses = 2;

  std::cout << "frags: " << iv << "\n";
  std::cout << "max_misses: " << max_misses << "\n";
  std::cout << "\n";
  std::cout << "Partial Sums Forward:\n";
  PartialSums ps_forward(iv, max_misses, PartialSums::forward_tag());
  std::cout << ps_forward << "\n";
  std::cout << "-------------------------------------\n\n\n";

  std::cout << "Partial Sums Reverse:\n";
  PartialSums ps_rev(iv, max_misses, PartialSums::reverse_tag());
  std::cout << ps_rev << "\n";
  std::cout << "-------------------------------------\n\n\n";

  std::cout << "Partial Sums Reverse 2:\n";
  PartialSums ps_rev2(iv_rev, max_misses, PartialSums::forward_tag());
  std::cout << ps_rev2 << "\n";
  std::cout << "-------------------------------------\n\n\n";


  return EXIT_SUCCESS;
}