// Run maligner and perform alignments in both row mode and column mode.
// Assert that the alignments are the same. Compare alignment times.

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

struct MapWrapper {

  MapWrapper(const Map& m, int num_missed_sites) :
    m_(m),
    md_(m_.name_, m_.frags_.size()),
    ps_(make_partial_sums(m_.frags_, num_missed_sites))
  {

  }

  Map m_;
  MapData md_;
  PartialSums ps_;

};

typedef unordered_map<string, MapWrapper> MapWrapperDB;

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
  MapWrapperDB map_db;
  for(auto i = ref_maps.begin(); i != ref_maps.end(); i++) {
    map_db.insert( MapWrapperDB::value_type(i->name_, MapWrapper(*i, maligner_dp::opt::ref_max_misses)) );
  }

 cerr << "Wrapped " << map_db.size() << " reference maps.\n";

 RowScoreMatrix sm_row, sm_row2;
 ColumnScoreMatrix sm_col, sm_col2;

 MapReader query_map_reader(maligner_dp::opt::query_maps_file);
 Map query_map;
 AlignmentVec alns;

 while(query_map_reader.next(query_map)) {

    MapWrapper qmw(query_map,align_opts.query_max_misses);

    const size_t num_query_frags = query_map.frags_.size();

    Timer query_timer;
    query_timer.start();

    for(MapWrapperDB::iterator ref_map_iter = map_db.begin();
        ref_map_iter != map_db.end();
        ref_map_iter++) {

      MapWrapper& rmw = ref_map_iter->second;
      const size_t num_ref_frags = rmw.m_.frags_.size();

      timer.start();

      sm_row.resize(num_query_frags + 1, num_ref_frags + 1);
      sm_row2.resize(num_query_frags + 1, num_ref_frags + 1);
      sm_col.resize(num_query_frags + 1, num_ref_frags + 1);
      sm_col2.resize(num_query_frags + 1, num_ref_frags + 1);

      timer.end();

      AlignTask<RowScoreMatrix> task_row(&qmw.md_, &rmw.md_,
        &qmw.m_.frags_, &rmw.m_.frags_, 
        &qmw.ps_, &rmw.ps_,
        0,
        &sm_row, &alns,
        true, //is_forward
        align_opts
      );

      AlignTask<ColumnScoreMatrix> task_col(&qmw.md_, &rmw.md_,
        &qmw.m_.frags_, &rmw.m_.frags_, 
        &qmw.ps_, &rmw.ps_,
        0,
        &sm_col, &alns,
        true, //is_forward
        align_opts
      );

      AlignTask<RowScoreMatrix> task_row2(&qmw.md_, &rmw.md_,
        &qmw.m_.frags_, &rmw.m_.frags_, 
        &qmw.ps_, &rmw.ps_,
        0,
        &sm_row2, &alns,
        true, //is_forward
        align_opts
      );

      AlignTask<ColumnScoreMatrix> task_col2(&qmw.md_, &rmw.md_,
        &qmw.m_.frags_, &rmw.m_.frags_, 
        &qmw.ps_, &rmw.ps_,
        0,
        &sm_col2, &alns,
        true, //is_forward
        align_opts
      );

      std::cout << "Aligning " << query_map.name_ << " to " << rmw.m_.name_ << "\n";


      //////////////////////////////////////////////////////
      // ALIGN BY ROW
      timer.start();
      fill_score_matrix_using_partials(task_row);
      timer.end();
      std::cout << "fill_score_matrix_using_partials, row order: " << timer << "\n";

      timer.start();
      fill_score_matrix_using_partials_with_breaks(task_row2);
      timer.end();
      std::cout << "fill_score_matrix_using_partials_with_breaks, row order: " << timer << "\n";

      // print_filled_by_row(std::cout, sm_row);
      // std::cout << "percent filled last row\t" << double(sm_row.countFilledByRow(sm_row.getNumRows()-1))/sm_row.getNumCols() << "\n";
      // std::cout << "percent filled\t" << sm_row.percentFilled() << "\n";

      // std::cout << sm_row << "\n";

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

      // print_filled_by_row(std::cout, sm_col);
      // std::cout << "percent filled last row\t" << double(sm_col.countFilledByRow(sm_col.getNumRows()-1))/sm_col.getNumCols() << "\n";
      // std::cout << "percent filled\t" << sm_col.percentFilled() << "\n";   
      // std::cout << "done.\n";
      // std::cout << sm_col << "\n";

      std::cout << "Row and column matrices match:\n"
        << "\t0: " << (sm_col == sm_col2) << "\n"
        << "\t1: " << (sm_row == sm_col) << "\n"
        << "\t2: " << (sm_row2 == sm_col2) << "\n"
        << "\t3: " << (sm_row == sm_row2) << "\n";

    }
    query_timer.end();

    alns.clear();
    
 }
 

  return EXIT_SUCCESS;

}
