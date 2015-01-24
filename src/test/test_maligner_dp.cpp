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

typedef ScoreMatrix<row_order_tag> ScoreMatrixType;
typedef AlignTask<ScoreMatrixType, Chi2SizingPenalty> AlignTaskType;

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
  
  // Build a database of reference maps. 
  MapVec ref_maps(read_maps(maligner_dp::opt::ref_maps_file));
  cerr << "Read " << ref_maps.size() << " reference maps.\n";

  // Store reference maps in an unordered map.
  MapWrapperDB map_db;
  for(auto i = ref_maps.begin(); i != ref_maps.end(); i++) {
    map_db.insert( MapWrapperDB::value_type(i->name_, MapWrapper(*i, maligner_dp::opt::ref_max_misses)) );
  }

 cerr << "Wrapped " << map_db.size() << " reference maps.\n";

 // Generate a single ScoreMatrix to use throughout this program.
 ScoreMatrixType sm;
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
      sm.resize(num_query_frags + 1, num_ref_frags + 1);
      timer.end();

      std::cout << "Resize: ("
        << sm.getNumRows() << "x" << sm.getNumCols()
        << ") capacity: " << sm.getCapacity()
        << " size: " << sm.getSize() 
        << " size-to-capacity: " << double(sm.getSize())/double(sm.getCapacity())
        << " " << timer << "\n";
      // Temporary loop to see speed without resizing
      // for(int i =0; i <5; i++) {
        AlignTaskType task(&qmw.md_, &rmw.md_,
          &qmw.m_.frags_, &rmw.m_.frags_, 
          &qmw.ps_, &rmw.ps_,
          0,
          &sm, &alns,
          true, //is_forward
          align_opts
        );

        std::cout << "Aligning " << query_map.name_ << " to " << rmw.m_.name_ << "\n";
        // make_best_alignments_using_partials(task);

        timer.start();
        fill_score_matrix_using_partials(task);
        timer.end();
        std::cout << "fill_score_matrix_using_partials: " << timer << "\n";
        print_filled_by_row(std::cout, sm);
        std::cout << "percent filled last row\t" << double(sm.countFilledByRow(sm.getNumRows()-1))/sm.getNumCols() << "\n";
        std::cout << "percent filled\t" << sm.percentFilled() << "\n";

        // timer.start();
        // fill_score_matrix_using_partials_no_size_penalty(task);
        // timer.end();
        // std::cout << "fill_score_matrix_using_partials_no_size_penalty: " << timer << "\n";
        // std::cout << "percent filled\t" << sm.percentFilled() << "\n";

        // timer.start();
        // fill_score_matrix_using_partials_with_size_penalty_class(task);
        // timer.end();
        // std::cout << "fill_score_matrix_using_partials_with_size_penalty_class: " << timer << "\n";
        // std::cout << "percent filled\t" << sm.percentFilled() << "\n";
        
        // timer.start();
        // fill_score_matrix_using_partials_with_cell_queue(task);
        // timer.end();
        // std::cout << "fill_score_matrix_using_partials_with_cell_queue: " << timer << "\n";
        // print_filled_by_row(std::cout, sm);
        // std::cout << "percent filled last row\t" << double(sm.countFilledByRow(sm.getNumRows()-1))/sm.getNumCols() << "\n";
        // std::cout << "percent filled\t" << sm.percentFilled() << "\n"; 

        timer.start();
        fill_score_matrix_using_partials_with_cell_mark(task);
        timer.end();
        std::cout << "fill_score_matrix_using_partials_with_cell_mark: " << timer << "\n";
        print_filled_by_row(std::cout, sm);
        std::cout << "percent filled last row\t" << double(sm.countFilledByRow(sm.getNumRows()-1))/sm.getNumCols() << "\n";
        std::cout << "percent filled\t" << sm.percentFilled() << "\n";   

        size_t white_count = sm.countColor(ScoreCellColor::WHITE);
        size_t green_count =sm.countColor(ScoreCellColor::GREEN);
        std::cout << "White: " << white_count << " (" << double(white_count)/sm.getSize() << ") "
                  << " Green: " << green_count << " (" << double(green_count)/sm.getSize() << ")" 
                  << "\n";     

        std::cout << "done.\n";
      // }

    }
    query_timer.end();

    std::cout << "Found " << alns.size() << " alignments. " << query_timer << "\n";

    alns.clear();
    
 }
 

  return EXIT_SUCCESS;

}
