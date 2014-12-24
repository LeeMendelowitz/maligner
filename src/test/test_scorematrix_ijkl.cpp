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
  size_t num_rows = 100;
  size_t row_delta = 3;
  size_t col_delta = 3;

  // num_cols = 10;
  // num_rows = 3;
  const int num_trials = 100;

  timer.start();
  ScoreMatrix sm(num_rows, num_cols);
  
  timer.end();

  cerr << "Size of score cell:  " << sizeof(ScoreCell) << " bytes\n";
  cerr << "Made score matrix: " << timer << "\n";

  timer.start();
  ScoreCellVec& sm_data = sm.getData();
  for(int i = 0; i < sm_data.size(); i++) {
    sm_data[i].score_ = i;
  }
  timer.end();

  cerr << "Set score values..." <<  timer << "\n";


  ///////////////////////////////////////////////////////////
  // Try to go straight across without using get cell
  {
    cerr << "col row col row"
         << "(" << num_rows << "x" << num_cols <<")";
    double sum = 0.0;
    timer.start();
    for(int itrial = 0; itrial < num_trials; itrial++) {

      for(size_t col = 0; col < num_cols; col++) {
        size_t col_offset = col*num_rows;

        int max_col_delta = min(col_delta, col);
        for(size_t row = 0; row < num_rows; row++) {

          int max_row_delta = min(row_delta, row);

            #if DEBUG > 0
            std::cerr << " col: " << col
                      << " row: " << row 
                      << " max_col_delta: " << max_col_delta 
                      << " max_row_delta: " << max_row_delta
                      << "\n";
            #endif
          for(int k = 1; k <= max_col_delta; k++) {
            size_t col_offset_delta = k*num_rows;
          for(int l = 1; l <= max_row_delta; l++) {

              #if DEBUG > 0
              std::cerr << " col: " << col
                        << " row: " << row 
                        << " k: " << k
                        << " l: " << l
                        << " row_sel: " << row - l
                        << " col_sel: " << col - k << "\n";
              #endif

              const ScoreCell* cell = &sm_data[col_offset - col_offset_delta + row - l];         
              sum += cell->score_;
            }
          }
        }
      }

    }
    timer.end();
    cerr << "sum=" << sum << "\n";
    cerr << "done: " << timer << "\n";
  }
  cerr << "*******************************************************\n";
  //////////////////////////////////////////////////////////////
  {
    cerr << "col row col row"
         << "(" << num_rows << "x" << num_cols <<")";
    double sum = 0.0;
    timer.start();
    for(int itrial = 0; itrial < num_trials; itrial++) {

      for(size_t col = 0; col < num_cols; col++) {
        int max_col_delta = min(col_delta, col);
        for(size_t row = 0; row < num_rows; row++) {

          int max_row_delta = min(row_delta, row);

            #if DEBUG > 0
            std::cerr << " col: " << col
                      << " row: " << row 
                      << " max_col_delta: " << max_col_delta 
                      << " max_row_delta: " << max_row_delta
                      << "\n";
            #endif
          for(int k = 1; k <= max_col_delta; k++) {
          for(int l = 1; l <= max_row_delta; l++) {

              #if DEBUG > 0
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
    cerr << "done iterator down columns: " << timer << "\n";
  }
  cerr << "*******************************************************\n";

  {
    cerr << "col row row col" 
         << "(" << num_rows << "x" << num_cols <<")";
    double sum = 0.0;
    timer.start();
    for(int itrial = 0; itrial < num_trials; itrial++) {

      for(size_t col = 0; col < num_cols; col++) {
        int max_col_delta = min(col_delta, col);
        for(size_t row = 0; row < num_rows; row++) {

          int max_row_delta = min(row_delta, row);

            #if DEBUG > 0
            std::cerr << " col: " << col
                      << " row: " << row 
                      << " max_col_delta: " << max_col_delta 
                      << " max_row_delta: " << max_row_delta
                      << "\n";
            #endif

          for(int l = 1; l <= max_row_delta; l++) {
          for(int k = 1; k <= max_col_delta; k++) {
            

              #if DEBUG > 0
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
    cerr << "column row row column: " << timer << "\n";
  }
  cerr << "*******************************************************\n";

  //////////////////////////////////////////////////////////////
  {
    cerr << "(row col row col) "
         << "(" << num_rows << "x" << num_cols <<")";
    double sum = 0.0;
    timer.start();
    for(int itrial = 0; itrial < num_trials; itrial++) {

      for(size_t row = 0; row < num_rows; row++) {
        int max_row_delta = min(row_delta, row);
        for(size_t col = 0; col < num_cols; col++) {
          int max_col_delta = min(col_delta, col);

          #if DEBUG > 0
            std::cerr << " col: " << col
                      << " row: " << row 
                      << " max_col_delta: " << max_col_delta 
                      << " max_row_delta: " << max_row_delta
                      << "\n";
          #endif

          for(int k = 1; k <= max_row_delta; k++) {
            for(int l = 1; l <= max_col_delta; l++) {

              #if DEBUG > 0
              std::cerr << " col: " << col
                        << " row: " << row 
                        << " k: " << k
                        << " l: " << l
                        << " row_sel: " << row - k
                        << " col_sel: " << col - l << "\n";
              #endif

              const ScoreCell* cell = sm.getCell(row - k, col - l);
              sum += cell->score_;
            }
          }
        }
      }

    }
    timer.end();
    cerr << "sum=" << sum << "\n";
    cerr << "done iterator across rows: " << timer << "\n";
  }
  cerr << "*******************************************************\n";

  return EXIT_SUCCESS;

}
