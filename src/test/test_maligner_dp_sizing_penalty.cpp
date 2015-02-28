// Run maligner and perform alignments in both row mode and column mode,
// and with different sizing penalty functions.
//
// Compare the run times!

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

typedef ScoreMatrix<row_order_tag> RowScoreMatrix;
typedef ScoreMatrix<column_order_tag> ColumnScoreMatrix;

typedef unordered_map<string, RefMapWrapper> RefMapWrapperDB;

////////////////////////////////////////////////////////////////////////
// Different methods of alignment to test


struct FillWithPartialsBreaks {

 template<class AlignTaskType>
 void operator()(AlignTaskType& align_task) {

  // std::cout << "FillWithPartialsBreaks!\n";
  fill_score_matrix_using_partials_with_breaks(align_task);
 }

};

struct FillWithPartialsBreaksHardcode {

 template<class AlignTaskType>
 void operator()(AlignTaskType& align_task) {
  // std::cout << "FillWithHardcode!\n";
  fill_score_matrix_using_partials_with_breaks_hardcode_penalty(align_task);
 }

};

struct FillWithPartialsBreaksHardcodeMaxMiss {

 template<class AlignTaskType>
 void operator()(AlignTaskType& align_task) {
  // std::cout << "FillWithHardcode!\n";
  fill_score_matrix_using_partials_with_breaks_hardcode_penalty_max_miss(align_task);
 }

};

struct FillWithPartialsBreaksHardcode_ijlk {

 template<class AlignTaskType>
 void operator()(AlignTaskType& align_task) {
  fill_score_matrix_using_partials_with_breaks_hardcode_penalty_ijlk(align_task);
 }

};



////////////////////////////////////////



template<class ScoreMatrixType, class SizingPenaltyType, class AlignFunction>
void align_to_all_ref(RefMapWrapperDB& ref_map_db,
                      QueryMapWrapper& qmw,
                      ScoreMatrixType* sm,
                      AlignmentVec& alns,
                      SizingPenaltyType,
                      AlignOpts& align_opts,
                      AlignFunction& align_func)
{
  // Align the query to all references, but only in the forward direction
  const size_t num_query_frags = qmw.num_frags();

  for(RefMapWrapperDB::iterator ref_map_iter = ref_map_db.begin();
          ref_map_iter != ref_map_db.end();
          ref_map_iter++) {

      RefMapWrapper& rmw = ref_map_iter->second;
      const size_t num_ref_frags = rmw.num_frags();

      sm->resize(num_query_frags + 1, num_ref_frags + 1);

      // std::cout << "num_query_frags:  " << num_query_frags << " num_ref_frags: " << num_ref_frags << "\n";

      AlignTask<ScoreMatrixType, SizingPenaltyType> align_task (
        &qmw.md_,
        &rmw.md_,
        &qmw.m_.frags_,
        &rmw.m_.frags_, 
        &qmw.ps_forward_,
        &rmw.ps_,
        &rmw.sd_inv_,
        sm, &alns,
        true, //is_forward
        align_opts
      );

      align_func(align_task);

  }

}



int main(int argc, char* argv[]) {

  maligner_dp::opt::program_name = argv[0];
  parse_args(argc, argv);
  Timer timer;

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
                       maligner_dp::opt::min_alignment_spacing,
                       maligner_dp::opt::neighbor_delta,
                       maligner_dp::opt::query_is_bounded,
                       maligner_dp::opt::ref_is_bounded);

  // Override the max_chunk_sizing_error, because we need the ScoreMatrix to be densely populated
  // to compare the effect of using Chi2 Sizing error vs no sizing error. The sizing error can bound
  // the search space when max_chunk_sizing_error is set.
  //align_opts.max_chunk_sizing_error = std::numeric_limits<double>::max();

  cerr << align_opts << "\n";

  // Build a database of reference maps. 
  MapVec ref_maps(read_maps(maligner_dp::opt::ref_maps_file));
  cerr << "Read " << ref_maps.size() << " reference maps.\n";

  // Store reference maps in an unordered map.
  RefMapWrapperDB ref_map_db;
  for(auto i = ref_maps.begin(); i != ref_maps.end(); i++) {
    ref_map_db.insert( RefMapWrapperDB::value_type(i->name_,
      RefMapWrapper(*i, maligner_dp::opt::ref_max_misses,
        maligner_dp::opt::sd_rate,
        maligner_dp::opt::min_sd)) );
  }

  cerr << "Wrapped " << ref_map_db.size() << " reference maps.\n";

  RowScoreMatrix sm_row;
  ColumnScoreMatrix sm_col;

  MapReader query_map_reader(maligner_dp::opt::query_maps_file);
  Map query_map;
  AlignmentVec alns;

  FillWithPartialsBreaks fill_with_partials_breaks;
  FillWithPartialsBreaksHardcode fill_hardcode;
  FillWithPartialsBreaksHardcode_ijlk fill_hardcode_ijlk;
  FillWithPartialsBreaksHardcodeMaxMiss fill_hardcode_max_miss;

  while(query_map_reader.next(query_map)) {

    QueryMapWrapper qmw(query_map, align_opts.query_max_misses);

    Timer query_timer;

    std::cout << "******************************************************" << std::endl;


    query_timer.start();
    align_to_all_ref(ref_map_db, qmw, &sm_row, alns, Chi2SizingPenalty(), align_opts, fill_with_partials_breaks);
    query_timer.end();
    std::cout << qmw.m_.name_ << " Row, chi2 sizing error: " << query_timer << std::endl;

    query_timer.start();
    align_to_all_ref(ref_map_db, qmw, &sm_row, alns, NoSizingPenalty(), align_opts, fill_hardcode);
    query_timer.end();
    std::cout << qmw.m_.name_ << " Row, hardcode sizing error: " << query_timer << std::endl;

    query_timer.start();
    align_to_all_ref(ref_map_db, qmw, &sm_row, alns, NoSizingPenalty(), align_opts, fill_hardcode_ijlk);
    query_timer.end();
    std::cout << qmw.m_.name_ << " Row, hardcode_ijlk sizing error: " << query_timer << std::endl;

    query_timer.start();
    align_to_all_ref(ref_map_db, qmw, &sm_row, alns, NoSizingPenalty(), align_opts, fill_hardcode_max_miss);
    query_timer.end();
    std::cout << qmw.m_.name_ << " Row, hardcode_max_miss sizing error: " << query_timer << std::endl;


    query_timer.start();
    align_to_all_ref(ref_map_db, qmw, &sm_col, alns, Chi2SizingPenalty(), align_opts, fill_with_partials_breaks);
    query_timer.end();
    std::cout << qmw.m_.name_ << " Column, chi2 sizing error: " << query_timer << std::endl;

    query_timer.start();
    align_to_all_ref(ref_map_db, qmw, &sm_col, alns, NoSizingPenalty(), align_opts, fill_hardcode);
    query_timer.end();
    std::cout << qmw.m_.name_ << " Column, hardcode sizing error: " << query_timer << std::endl;


    // query_timer.start();
    // align_to_all_ref(ref_map_db, qmw, &sm_row, alns, NoSizingPenalty(), align_opts, fill_with_partials_breaks);
    // query_timer.end();
    // std::cout << qmw.m_.name_ << " Row, no sizing error: " << query_timer << std::endl;

    // query_timer.start();
    // align_to_all_ref(ref_map_db, qmw, &sm_col, alns, NoSizingPenalty(), align_opts, fill_with_partials_breaks);
    // query_timer.end();
    // std::cout << qmw.m_.name_ << " Column, no sizing error: " << query_timer << std::endl;



    alns.clear();
      
  }
 

  return EXIT_SUCCESS;

}
