#include <iostream>
#include <algorithm>
#include <vector>
#include <string>
#include <unordered_map>
#include <sstream>
#include <fstream>
#include <chrono>
#include <getopt.h>
#include <limits>

// kmer_match includes
// #include "map.h"
// #include "map_reader.h"
// #include "map_frag.h"
// #include "map_frag_db.h"
// #include "error_model.h"
// #include "map_chunk.h"
// #include "ref_alignment.h"
// #include "map_chunk_db.h"

// dp includes
#include "align.h"
#include "utils.h"
#include "ScoreMatrix.h"
#include "ScoreCell.h"

// common includes
#include "timer.h"

int main(int argc, char* argv[]) {

  using lmm_utils::Timer;
  using namespace maligner_dp;
  using namespace std;

  Timer timer;
  
  typedef ScoreMatrix<row_order_tag> ScoreMatrixType;
 
  std::cout << "sizeof(ScoreCell): " << sizeof(ScoreCell) << "\n"
            << "sizeof(ScoreMatrix): " << sizeof(ScoreMatrixType) << "\n";

  ScoreCell cell1, cell2;

  std::cout << "Cell1: " << ScoreCellFullOutput(&cell1) << "\n";
  std::cout << "Cell2: " << ScoreCellFullOutput(&cell2) << "\n";
  std::cout << "Cell1==Cell2: " << (cell1 == cell2) << "\n";

  const double INF = std::numeric_limits<double>::infinity();
  std::cout << "INF == INF: " << (INF==INF) << "\n";
  std::cout << "-INF == -INF: " << (-INF==-INF) << "\n";
  return EXIT_SUCCESS;

}
