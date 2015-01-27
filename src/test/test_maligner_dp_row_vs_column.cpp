// Run maligner and perform alignments in both row mode and column mode.
// Assert that the alignments matrixces are filled out the same.

#include <iostream>
#include <algorithm>
#include <vector>
#include <string>
#include <unordered_map>
#include <sstream>
#include <fstream>
#include <chrono>
#include <getopt.h>
#include <exception>

// kmer_match includes
#include "map.h"
#include "map_reader.h"
#include "map_frag.h"
#include "map_frag_db.h"
#include "error_model.h"
#include "map_chunk.h"
#include "ref_alignment.h"
#include "map_chunk_db.h"

// dp includes
#include "map_data.h"
#include "alignment.h"
#include "align.h"
#include "utils.h"
#include "ScoreMatrix.h"

// common includes
#include "timer.h"
#include "common_defs.h"

using std::string;
using std::unordered_map;
using namespace maligner_maps;
// using namespace std;

#include "maligner_dp_includes.h"

// Wrapper around all of the Map structures we need to store
// in order to perform DP alignments.
using namespace kmer_match;
using namespace maligner_dp;
// using kmer_match::Map;

// using maligner_dp::MapData;
// using maligner_dp::PartialSums;
// using maligner_dp::ScoreMatrix
// using maligner_dp::make_partial_sums;
// using maligner_dp::AlignOpts;
using lmm_utils::Timer;

class ScoreMatrixMismatchError: public std::exception
{
  virtual const char* what() const throw()
  {
    return "ScoreMatrices do not match!";
  }
};

template<class SM1, class SM2>
void check_score_matrix(SM1& sm1, SM2& sm2) {
  if (! (sm1 == sm2)) {
    throw ScoreMatrixMismatchError();
  }
}


typedef ScoreMatrix<row_order_tag> RowScoreMatrix;
typedef ScoreMatrix<column_order_tag> ColumnScoreMatrix;


typedef unordered_map<string, RefMapWrapper> RefMapWrapperDB;


void print_partial_sums(const PartialSums& ps, size_t num_frags, int max_misses) {
  for(size_t i = 0; i < num_frags; i++) {
    for(int miss = 0; miss < max_misses + 1; miss++) {
      std::cerr << "i: " << i << " miss: " << miss << " " << ps(i, miss) << "\n";
    }
  }
}

int main(int argc, char* argv[]) {

  maligner_dp::opt::program_name = argv[0];
  parse_args(argc, argv);
  Timer timer;

  // cerr << "................................................\n"
  //      << "MALIGNER settings:\n"
  //      << "\tquery_maps_file: " << maligner_dp::opt::query_maps_file << "\n"
  //      << "\tref_maps_file: " << maligner_dp::opt::ref_maps_file << "\n"
  //      << "\tmax. consecutive unmatched sites: " << maligner_dp::opt::max_unmatched_sites << "\n"
  //      << "\tmax. unmatched rate: " << maligner_dp::opt::max_unmatched_rate << "\n"
  //      << "\trelative_error: " << maligner_dp::opt::rel_error << "\n"
  //      << "\tabsolute_error: " << maligner_dp::opt::min_abs_error << "\n"
  //      << "\tminimum query frags: " << maligner_dp::opt::min_frag << "\n"
  //      << "\tmax matches per query: " << maligner_dp::opt::max_match << "\n"
  //      << "................................................\n\n";

  AlignOpts align_opts(maligner_dp::opt::query_miss_penalty,
                       maligner_dp::opt::ref_miss_penalty,
                       maligner_dp::opt::query_max_misses,
                       maligner_dp::opt::ref_max_misses,
                       maligner_dp::opt::sd_rate,
                       maligner_dp::opt::min_sd,
                       maligner_dp::opt::max_chunk_sizing_error,
                       maligner_dp::opt::ref_max_miss_rate,
                       maligner_dp::opt::query_max_miss_rate,
                       maligner_dp::opt::alignments_per_reference,
                       maligner_dp::opt::alignments_per_reference,
                       maligner_dp::opt::min_alignment_spacing,
                       maligner_dp::opt::neighbor_delta,
                       maligner_dp::opt::query_is_bounded,
                       maligner_dp::opt::ref_is_bounded);
  
  cerr << align_opts << "\n";
  
  // Build a database of reference maps. 
  MapVec ref_maps(read_maps(maligner_dp::opt::ref_maps_file));
  cerr << "Read " << ref_maps.size() << " reference maps.\n";

  // Store reference maps in an unordered map.
  RefMapWrapperDB map_db;
  for(auto i = ref_maps.begin(); i != ref_maps.end(); i++) {
    map_db.insert( RefMapWrapperDB::value_type(i->name_, RefMapWrapper(*i,
      maligner_dp::opt::ref_max_misses,
      maligner_dp::opt::sd_rate,
      maligner_dp::opt::min_sd)) );
  }

 cerr << "Wrapped " << map_db.size() << " reference maps.\n";

 RowScoreMatrix sm_row, sm_row2, sm_row3, sm_row4, sm_row5;
 ColumnScoreMatrix sm_col, sm_col2, sm_col3, sm_col4;

 MapReader query_map_reader(maligner_dp::opt::query_maps_file);
 Map query_map;
 AlignmentVec alns;

 while(query_map_reader.next(query_map)) {

    QueryMapWrapper qmw(query_map,align_opts.query_max_misses);

    // print_partial_sums(qmw.ps_forward_, qmw.m_.frags_.size(), align_opts.query_max_misses);

    const size_t num_query_frags = query_map.frags_.size();

    Timer query_timer;
    query_timer.start();

    for(auto ref_map_iter = map_db.begin();
        ref_map_iter != map_db.end();
        ref_map_iter++) {

      RefMapWrapper& rmw = ref_map_iter->second;
      const size_t num_ref_frags = rmw.m_.frags_.size();

      AlignTask<RowScoreMatrix, Chi2SizingPenalty> task_row(&qmw.md_, &rmw.md_,
        &qmw.m_.frags_, &rmw.m_.frags_, 
        &qmw.ps_forward_, &rmw.ps_, &rmw.sd_inv_,
        0,
        &sm_row, &alns,
        true, //is_forward
        align_opts
      );

      AlignTask<ColumnScoreMatrix, Chi2SizingPenalty> task_col(&qmw.md_, &rmw.md_,
        &qmw.m_.frags_, &rmw.m_.frags_, 
        &qmw.ps_forward_, &rmw.ps_, &rmw.sd_inv_,
        0,
        &sm_col, &alns,
        true, //is_forward
        align_opts
      );

      AlignTask<RowScoreMatrix, Chi2SizingPenalty> task_row2(&qmw.md_, &rmw.md_,
        &qmw.m_.frags_, &rmw.m_.frags_, 
        &qmw.ps_forward_, &rmw.ps_, &rmw.sd_inv_,
        0,
        &sm_row2, &alns,
        true, //is_forward
        align_opts
      );

      AlignTask<ColumnScoreMatrix, Chi2SizingPenalty> task_col2(&qmw.md_, &rmw.md_,
        &qmw.m_.frags_, &rmw.m_.frags_, 
        &qmw.ps_forward_, &rmw.ps_, &rmw.sd_inv_,
        0,
        &sm_col2, &alns,
        true, //is_forward
        align_opts
      );

      AlignTask<RowScoreMatrix, Chi2SizingPenalty> task_row3(&qmw.md_, &rmw.md_,
        &qmw.m_.frags_, &rmw.m_.frags_, 
        &qmw.ps_forward_, &rmw.ps_, &rmw.sd_inv_,
        0,
        &sm_row3, &alns,
        true, //is_forward
        align_opts
      );

      AlignTask<ColumnScoreMatrix, Chi2SizingPenalty> task_col3(&qmw.md_, &rmw.md_,
        &qmw.m_.frags_, &rmw.m_.frags_, 
        &qmw.ps_forward_, &rmw.ps_, &rmw.sd_inv_,
        0,
        &sm_col3, &alns,
        true, //is_forward
        align_opts
      );

      AlignTask<RowScoreMatrix, Chi2SizingPenalty> task_row4(&qmw.md_, &rmw.md_,
        &qmw.m_.frags_, &rmw.m_.frags_, 
        &qmw.ps_forward_, &rmw.ps_, &rmw.sd_inv_,
        0,
        &sm_row4, &alns,
        true, //is_forward
        align_opts
      );

      AlignTask<RowScoreMatrix, Chi2SizingPenalty> task_row5(&qmw.md_, &rmw.md_,
        &qmw.m_.frags_, &rmw.m_.frags_, 
        &qmw.ps_forward_, &rmw.ps_, &rmw.sd_inv_,
        0,
        &sm_row5, &alns,
        true, //is_forward
        align_opts
      );

      AlignTask<ColumnScoreMatrix, Chi2SizingPenalty> task_col4(&qmw.md_, &rmw.md_,
        &qmw.m_.frags_, &rmw.m_.frags_, 
        &qmw.ps_forward_, &rmw.ps_, &rmw.sd_inv_,
        0,
        &sm_col4, &alns,
        true, //is_forward
        align_opts
      );



      std::cout << "Aligning " << query_map.name_ << " to " << rmw.m_.name_ << "\n";


      ////////////////////////////////////////////////////
      // ALIGN BY ROW
      timer.start();
      fill_score_matrix_using_partials(task_row);
      timer.end();
      std::cout << "fill_score_matrix_using_partials, row order: " << timer << "\n";

      timer.start();
      fill_score_matrix_using_partials_with_breaks(task_row2);
      timer.end();
      std::cout << "fill_score_matrix_using_partials_with_breaks, row order: " << timer << "\n";
      
      timer.start();
      fill_score_matrix_with_breaks(task_row3);
      timer.end();
      std::cout << "fill_score_matrix_with_breaks, row order: " << timer << "\n";

      timer.start();
      fill_score_matrix_using_partials_with_breaks_hardcode_penalty(task_row4);
      timer.end();
      std::cout << "fill_score_matrix_using_partials_with_breaks_hardcode_penalty, row order: " << timer << "\n";

      timer.start();
      fill_score_matrix_using_partials_with_breaks_hardcode_penalty_ijlk(task_row5);
      timer.end();
      std::cout << "fill_score_matrix_using_partials_with_breaks_hardcode_penalty_ijlk, row order: " << timer << "\n";




      //////////////////////////////////////////////////////
      // ALIGN BY COLUMN

      timer.start();
      fill_score_matrix_using_partials(task_col);
      timer.end();
      std::cout << "fill_score_matrix_using_partials, col order: " << timer << "\n";

      timer.start();
      fill_score_matrix_using_partials_with_breaks(task_col2);
      timer.end();
      std::cout << "fill_score_matrix_using_partials_with_breaks, col order: " << timer << "\n";

      timer.start();
      fill_score_matrix_with_breaks(task_col3);
      timer.end();
      std::cout << "fill_score_matrix_with_breaks, col order: " << timer << "\n";
      
      timer.start();
      fill_score_matrix_using_partials_with_breaks_hardcode_penalty(task_col4);
      timer.end();
      std::cout << "fill_score_matrix_using_partials_with_breaks_hardcode_penalty, col order: " << timer << "\n";



      std::cout << "Row and column matrices match:\n"
        << "\ta: " << (sm_col == sm_col2) << "\n"
        << "\tb: " << (sm_col == sm_col3) << "\n"
        << "\tb2: " << (sm_col == sm_col4) << "\n"
        << "\tc: " << (sm_row == sm_row2) << "\n"
        << "\td: " << (sm_row == sm_row3) << "\n"
        << "\tf: " << (sm_row == sm_row4) << "\n"
        << "\tg: " << (sm_row == sm_row5) << "\n"
        << "\te: " << (sm_row == sm_col) << "\n";

      check_score_matrix(sm_col, sm_col2);
      check_score_matrix(sm_col, sm_col3);
      check_score_matrix(sm_col, sm_col4);

      check_score_matrix(sm_row, sm_row2);
      check_score_matrix(sm_row, sm_row3);
      check_score_matrix(sm_row, sm_row4);
      check_score_matrix(sm_row, sm_row5);

      check_score_matrix(sm_col, sm_row);


    }
    query_timer.end();

    alns.clear();
    
 }
 

  return EXIT_SUCCESS;

}
