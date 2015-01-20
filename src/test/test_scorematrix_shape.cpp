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

// #define DEBUG

int main(int argc, char* argv[]) {

  using lmm_utils::Timer;
  using namespace maligner_dp;
  using namespace std;

  Timer timer;

  cerr << "Size of score cell:  " << sizeof(ScoreCell) << " bytes\n";

  size_t num_rows = 100;
  size_t num_cols = 50000;

  size_t row_delta = 2;
  size_t col_delta = 5;

  // num_cols = 10;
  // num_rows = 3;
  const int num_trials = 100;


  cerr << "*******************************************************\n";
  //////////////////////////////////////////////////////////////
  {

    timer.start();
    ScoreMatrix<column_order_tag> sm(num_rows, num_cols);
    timer.end();
    cerr << "Made column order score matrix: " << timer << "\n";

    timer.start();
    ScoreCellVec& sm_data = sm.getData();
    for(int i = 0; i < sm_data.size(); i++) {
      sm_data[i].score_ = 1.0;
    }
    timer.end();
    cerr << "Set score values..." <<  timer << "\n";

    cerr << "col row col row using get cell\n"
       << "(" << sm.getNumRows() << "x" << sm.getNumCols() <<") size=" << sm.getSize() << "\n";
    double sum = 0.0;
    timer.start();
    for(int itrial = 0; itrial < num_trials; itrial++) {

      for(size_t col = 0; col < num_cols; col++) {
        int max_col_delta = min(col_delta, col);
        for(size_t row = 0; row < num_rows; row++) {

          int max_row_delta = min(row_delta, row);

            #ifdef DEBUG
            std::cerr << " col: " << col
                      << " row: " << row 
                      << " max_col_delta: " << max_col_delta 
                      << " max_row_delta: " << max_row_delta
                      << "\n";
            #endif
          for(int k = 1; k <= max_col_delta; k++) {
          for(int l = 1; l <= max_row_delta; l++) {

              #ifdef DEBUG
              std::cerr << " col: " << col
                        << " row: " << row 
                        << " k: " << k
                        << " l: " << l
                        << " row_sel: " << row - l
                        << " col_sel: " << col - k << "\n";
              #endif

              const ScoreCell* cell = sm.getCell(row - l, col - k);
              sum += cell->score_;
            }
          }
        }
      }

    }
    timer.end();
    cerr << "sum=" << sum << "\n";
    cerr << "done " << timer << "\n";
  }
  cerr << "*******************************************************\n";

  //////////////////////////////////////////////////////////////
  {

    timer.start();
    ScoreMatrix<row_order_tag> sm(num_rows, num_cols);
    timer.end();
    cerr << "Made row order score matrix: " << timer << "\n";

    timer.start();
    ScoreCellVec& sm_data = sm.getData();
    for(int i = 0; i < sm_data.size(); i++) {
      sm_data[i].score_ = 1.0;
    }
    timer.end();
    cerr << "Set score values..." <<  timer << "\n";

    cerr << "row col row col using get cell\n"
       << "(" << sm.getNumRows() << "x" << sm.getNumCols() <<") size=" << sm.getSize() << "\n";
    double sum = 0.0;
    timer.start();
    for(int itrial = 0; itrial < num_trials; itrial++) {

      for(size_t row = 0; row < num_rows; row++) {

        int max_row_delta = min(row_delta, row);

        for(size_t col = 0; col < num_cols; col++) {

          int max_col_delta = min(col_delta, col);

            #ifdef DEBUG
            std::cerr << " col: " << col
                      << " row: " << row 
                      << " max_col_delta: " << max_col_delta 
                      << " max_row_delta: " << max_row_delta
                      << "\n";
            #endif

          for(int l = 1; l <= max_row_delta; l++) {         
          for(int k = 1; k <= max_col_delta; k++) {
          

              #ifdef DEBUG
              std::cerr << " col: " << col
                        << " row: " << row 
                        << " k: " << k
                        << " l: " << l
                        << " row_sel: " << row - l
                        << " col_sel: " << col - k << "\n";
              #endif

              const ScoreCell* cell = sm.getCell(row - l, col - k);
              sum += cell->score_;
            }
          }
        }
      }

    }
    timer.end();
    cerr << "sum=" << sum << "\n";
    cerr << "done " << timer << "\n";
  }
  cerr << "*******************************************************\n";


  return EXIT_SUCCESS;

}
