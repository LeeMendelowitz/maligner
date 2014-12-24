#include <iostream>
#include <algorithm>
#include <vector>
#include <string>
#include <unordered_map>
#include <sstream>
#include <fstream>
#include <chrono>
#include <getopt.h>

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

// common includes
#include "timer.h"

int main(int argc, char* argv[]) {

  using lmm_utils::Timer;
  using namespace maligner_dp;
  using namespace std;

  Timer timer;
  
  size_t num_cols = 50000;
  size_t num_rows = 50;
  const int num_trials = 100;

  timer.start();
  ScoreMatrix sm(num_rows, num_cols);
  
  timer.end();

  cerr << "Made score matrix: " << timer << "\n";

  timer.start();
  ScoreCellVec& sm_data = sm.getData();
  for(int i = 0; i < sm_data.size(); i++) {
    sm_data[i].score_ = i;
  }
  timer.end();

  cerr << "Set score values..." <<  timer << "\n";


  //////////////////////////////////////////////////////////////
  {
    cerr << "Iterating straight across...";
    double sum = 0.0;    
    timer.start();
    for(int itrial = 0; itrial < num_trials; itrial++) {
      for(size_t i = 0; i < sm_data.size(); i++) {
        sum += sm_data[i].score_;
      }
    }
    timer.end();
    cerr << "sum=" << sum << "\n";
    cerr << "done iterator straight across: " << timer << "\n";
  }
  //////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////
  {
    cerr << "Iterating down columns";
    double sum = 0.0;
    timer.start();
    for(int itrial = 0; itrial < num_trials; itrial++) {
      for(size_t col = 0; col < num_cols; col++) {
        for(size_t row = 0; row < num_rows; row++) {
          const ScoreCell* cell = sm.getCell(row, col);
          sum += cell->score_;
        }
      }
    }
    timer.end();
    cerr << "sum=" << sum << "\n";
    cerr << "done iterator down columns: " << timer << "\n";
  }

  //////////////////////////////////////////////////////////////
  {
    cerr << "Iterating across rows";
    double sum = 0.0;
    timer.start();
    for(int itrial = 0; itrial < num_trials; itrial++) {
      for(size_t row = 0; row < num_rows; row++) {
        for(size_t col = 0; col < num_cols; col++) {
          const ScoreCell* cell = sm.getCell(row, col);
          sum += cell->score_;
        }
      }
    }
    timer.end();
    cerr << "sum=" << sum << "\n";
    cerr << "done iterator across rows: " << timer << "\n";
  }

  return EXIT_SUCCESS;

}
