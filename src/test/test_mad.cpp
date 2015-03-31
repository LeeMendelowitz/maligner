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

// #define DEBUG

int main(int argc, char* argv[]) {

  using lmm_utils::Timer;
  using namespace std;
  typedef std::vector<double> DoubleVec;

  DoubleVec v1 {1.0, 2.0, 3.0};
  std::cout << "median v1: " << median(v1) << "\n";
  std::cout << "mad v1: " << mad(v1) << "\n";

  DoubleVec v2 {1.0};
  std::cout << "median v2: " << median(v2) << "\n";
  std::cout << "mad v2: " << mad(v2) << "\n";

  DoubleVec v3 {4.0, 1.0}; 
  std::cout << "median v3: " << median(v3) << "\n";
  std::cout << "mad v3: " << mad(v3) << "\n";

  return EXIT_SUCCESS;

}
