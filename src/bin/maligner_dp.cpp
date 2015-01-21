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
typedef AlignTask<ScoreMatrixType> AlignTaskType;

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

  using maligner_dp::Alignment;
  maligner_dp::opt::program_name = argv[0];
  parse_args(argc, argv);

  print_args(std::cerr);

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
 AlignmentVec alns_forward, alns_reverse, all_alignments;


 std::cout << AlignmentHeader();
 while(query_map_reader.next(query_map)) {

    all_alignments.clear();
    alns_forward.clear();
    alns_reverse.clear();

    const IntVec& query_frags_forward = query_map.frags_;
    const IntVec query_frags_reverse(query_frags_forward.rbegin(), query_frags_forward.rend());
    PartialSums qps_forward = make_partial_sums(query_frags_forward, align_opts.query_max_misses);
    PartialSums qps_reverse = make_partial_sums(query_frags_reverse, align_opts.query_max_misses);


    MapWrapper qmw(query_map,align_opts.query_max_misses);

    const size_t num_query_frags = query_map.frags_.size();

    Timer query_timer;
    query_timer.start();

    for(MapWrapperDB::iterator ref_map_iter = map_db.begin();
        ref_map_iter != map_db.end();
        ref_map_iter++) {


      MapWrapper& rmw = ref_map_iter->second;
      const size_t num_ref_frags = rmw.m_.frags_.size();
      // sm.resize(num_query_frags + 1, num_ref_frags + 1);

      // const IntVec* p_frags_forward = &query_frags_forward;
      // const IntVec* p_frags_reverse = &query_frags_reverse;

      AlignTaskType task_forward(&qmw.md_, &rmw.md_,
        &query_frags_forward, &rmw.m_.frags_, 
        &qps_forward, &rmw.ps_,
        0, // ref_offset
        &sm, &all_alignments,
        true, // is_forward
        align_opts
      );

      AlignTaskType task_reverse(&qmw.md_, &rmw.md_,
        &query_frags_reverse, &rmw.m_.frags_, 
        &qps_reverse, &rmw.ps_,
        0, // ref_offset
        &sm, &all_alignments,
        false, // is_forward
        align_opts
      );

      // std::cerr << "Aligning " << query_map.name_ << " to " << rmw.m_.name_ << "\n";
      Timer timer;

      // Align Forward
      {
        timer.start();
        // Alignment aln_forward = make_best_alignment_using_partials(task_forward);
        int num_alignments = make_best_alignments_using_partials(task_forward);
        timer.end();
        // std::cerr << "Num alignments forward: " << num_alignments
        //           << " " << timer << "\n";
      }

      // Align Reverse
      {
        timer.start();
        int num_alignments = make_best_alignments_using_partials(task_reverse);
        timer.end();
        // std::cerr << "Num alignments reverse: " << num_alignments
        //           << " " << timer << "\n";
      }
    }

    // Print the best alignments
    std::sort(all_alignments.begin(), all_alignments.end(), AlignmentRescaledScoreComp());
    const int max_ind = std::min(int(all_alignments.size()), opt::max_alignments);
    for(int i = 0; i < max_ind; i++) {
      print_alignment(std::cout, all_alignments[i]);
    }
    
    query_timer.end();    

    std::cerr << "Aligned query " << query_map.name_ << " to all references.\n"
              << query_timer << "\n";
    std::cerr << "*****************************************\n";

 }
 

  return EXIT_SUCCESS;

}
