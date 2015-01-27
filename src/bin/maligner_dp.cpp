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
using namespace maligner_maps;


using lmm_utils::Timer;

typedef ScoreMatrix<row_order_tag> ScoreMatrixType;
typedef AlignTask<ScoreMatrixType, Chi2SizingPenalty> AlignTaskType;



typedef unordered_map<string, RefMapWrapper> RefMapDB;

int main(int argc, char* argv[]) {

  using maligner_dp::Alignment;
  maligner_dp::opt::program_name = argv[0];
  parse_args(argc, argv);

  print_args(std::cerr);

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
                       maligner_dp::opt::alignments_per_reference,
                       maligner_dp::opt::min_alignment_spacing,
                       maligner_dp::opt::neighbor_delta,
                       maligner_dp::opt::query_is_bounded,
                       maligner_dp::opt::ref_is_bounded);

  // Build a database of reference maps. 
  MapVec ref_maps(read_maps(maligner_dp::opt::ref_maps_file));
  cerr << "Read " << ref_maps.size() << " reference maps.\n";

  // Store reference maps in an unordered map.
  RefMapDB ref_map_db;
  for(auto i = ref_maps.begin(); i != ref_maps.end(); i++) {
    ref_map_db.insert( RefMapDB::value_type(i->name_,
        RefMapWrapper(*i, maligner_dp::opt::ref_max_misses,
                          maligner_dp::opt::sd_rate,
                          maligner_dp::opt::min_sd)) );
  }

 cerr << "Wrapped " << ref_map_db.size() << " reference maps.\n";

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
    PartialSums qps_forward(query_frags_forward, align_opts.query_max_misses);
    PartialSums qps_reverse(query_frags_reverse, align_opts.query_max_misses);

    const MapWrapper qmw(query_map, align_opts.query_max_misses);

    const size_t num_query_frags = query_map.frags_.size();

    Timer query_timer;
    query_timer.start();

    for(auto ref_map_iter = ref_map_db.begin();
        ref_map_iter != ref_map_db.end();
        ref_map_iter++) {


      const RefMapWrapper& rmw = ref_map_iter->second;
      const size_t num_ref_frags = rmw.m_.frags_.size();
      // sm.resize(num_query_frags + 1, num_ref_frags + 1);

      // const IntVec* p_frags_forward = &query_frags_forward;
      // const IntVec* p_frags_reverse = &query_frags_reverse;

      AlignTaskType task_forward(const_cast<MapData*>(&qmw.md_),
        const_cast<MapData*>(&rmw.md_),
        &query_frags_forward, &rmw.m_.frags_, 
        &qps_forward, &rmw.ps_,
        &rmw.sd_inv_,
        0, // ref_offset
        &sm, &all_alignments,
        true, // is_forward
        align_opts
      );

      AlignTaskType task_reverse(const_cast<MapData*>(&qmw.md_),
        const_cast<MapData*>(&rmw.md_),
        &query_frags_reverse, &rmw.m_.frags_, 
        &qps_reverse, &rmw.ps_,
        &rmw.sd_inv_,
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
