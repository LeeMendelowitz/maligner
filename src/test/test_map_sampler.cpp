// Test whether the shape of the score matrix is important for the iteration.
// ScoreMatrix is still column order.
// For example: try 100 x 1000000 vs 10000 x 10000
//
// A ScoreCell is 32 bytes, so a column of 100 is 32kb, a scorematrix of 10000 is 3200kb.
// A scorematrix of 3200kb is likely to have more cache misses?

#include <iostream>
#include <algorithm>
#include <vector>
#include <string>
#include <unordered_map>
#include <sstream>
#include <fstream>
#include <chrono>
#include <getopt.h>



// dp includes
#include "common_math.h"

// common includes
#include "timer.h"

#include "map_sampler.h"

// #define DEBUG

template< class T >
std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec) {

  const size_t N = vec.size();

  for(size_t i = 0; i < N; i++) {
    os << vec[i] << " ";
  }

  os << "\n";

  return os;

}

int main(int argc, char* argv[]) {

  using lmm_utils::Timer;
  using namespace std;
  typedef std::vector<int> IntVec;

  IntVec v1 {1, 2, 3, 4, 5};
  std::cout << v1 << "\n";
  std::cout << "permute: " << permute(v1) << "\n";
  std::cout << "permute: " << permute(v1) << "\n";
  std::cout << "permute: " << permute(v1) << "\n";
  std::cout << "sample: " << sample_with_replacement(v1, 2) << "\n";
  std::cout << "sample: " << sample_with_replacement(v1, 3) << "\n";
  std::cout << "sample: " << sample_with_replacement(v1, 4) << "\n";


  return EXIT_SUCCESS;

}
