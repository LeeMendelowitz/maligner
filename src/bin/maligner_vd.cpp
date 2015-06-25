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
#include "map_sampler.h"

// dp includes  
#include "map_data.h"
#include "align.h"
#include "utils.h"
#include "ScoreMatrix.h"

// vd includes
#include "score_matrix_vd.h"
#include "score_matrix_vd_db.h"
#include "score_matrix_profile.h"

// common includes
#include "timer.h"
#include "common_defs.h"

using std::string;
using std::unordered_map;
using std::cerr;
using std::cout;
// using namespace std;


#define PACKAGE_NAME "maligner vd"
#include "maligner_vd_includes.h"



// Wrapper around all of the Map structures we need to store
// in order to perform DP alignments.
// using namespace kmer_match;
using namespace maligner_vd;
using namespace maligner_dp;
using namespace maligner_maps;
using maligner_dp::Alignment;
using maligner_dp::AlignmentHeader;
using maligner_dp::ScoreMatrix;
using maligner_dp::row_order_tag;

using lmm_utils::Timer;

typedef ScoreMatrix<row_order_tag> ScoreMatrixType;
typedef AlignTask<ScoreMatrixType, Chi2SizingPenalty> AlignTaskType;
typedef maligner_vd::RefScoreMatrixVD<ScoreMatrixType> RefScoreMatrixVDType;
typedef maligner_vd::RefScoreMatrixVDDB<RefScoreMatrixVDType> RefScoreMatrixDB;
typedef std::vector<RefScoreMatrixVDType> RefScoreMatrixVDVec;



int main(int argc, char* argv[]) {

  using maligner_dp::Alignment;
  maligner_vd::opt::program_name = argv[0];
  parse_args(argc, argv);

  print_args(std::cerr);

  Timer timer;

  AlignOpts align_opts(maligner_vd::opt::query_miss_penalty,
                       maligner_vd::opt::ref_miss_penalty,
                       maligner_vd::opt::query_max_misses,
                       maligner_vd::opt::ref_max_misses,
                       maligner_vd::opt::sd_rate,
                       maligner_vd::opt::min_sd,
                       maligner_vd::opt::max_chunk_sizing_error,
                       maligner_vd::opt::ref_max_miss_rate,
                       maligner_vd::opt::query_max_miss_rate,
                       maligner_vd::opt::alignments_per_reference,
                       maligner_vd::opt::min_alignment_spacing,
                       maligner_vd::opt::neighbor_delta,
                       maligner_vd::opt::query_is_bounded, // Perhaps this should be part of the MapData instead of AlignOpts
                       maligner_vd::opt::ref_is_bounded, // Perhaps this should be part of the MapData instead of AlignOpts
                       maligner_vd::opt::query_rescaling);

  // Build a database of reference maps. 
  MapVec ref_maps(read_maps(maligner_vd::opt::ref_maps_file));
  cerr << "Read " << ref_maps.size() << " reference maps.\n";

  RefScoreMatrixVDVec ref_score_matrix_vd_vec;
  RefScoreMatrixDB ref_score_matrix_db;
  for(auto& rm : ref_maps ) {

    RefMapWrapper rmw(rm, maligner_vd::opt::reference_is_circular, 
                          maligner_vd::opt::ref_max_misses,
                          maligner_vd::opt::sd_rate,
                          maligner_vd::opt::min_sd);

    ref_score_matrix_db.add_ref_map(std::move(rmw));

  }

  cerr << "Wrapped " << ref_score_matrix_db.size() << " reference maps.\n";

  MapReader query_map_reader(maligner_vd::opt::query_maps_file);
  Map query_map;
  AlignmentVec alns_forward, alns_reverse, all_alignments;

  ////////////////////////////////////////////////////////
  // Open output files for the different alignment types.
  // auto fout_rf_qf = std::ofstream("rf_qf.aln");
  // auto fout_rf_qr = std::ofstream("rf_qr.aln");
  // auto fout_rr_qf = std::ofstream("rr_qf.aln");
  // auto fout_rr_qr = std::ofstream("rr_qr.aln");

  // fout_rf_qf << AlignmentHeader();
  // fout_rf_qr << AlignmentHeader();
  // fout_rr_qf << AlignmentHeader();
  // fout_rr_qr << AlignmentHeader();

  std::cout << maligner_vd::ScoreMatrixRecordHeader() << "\n";

  while(query_map_reader.next(query_map)) {

    all_alignments.clear();
    alns_forward.clear();
    alns_reverse.clear();

    const IntVec& query_frags_forward = query_map.frags_;

    if(query_map.frags_.size() < maligner_vd::opt::min_query_frags) {

      if(maligner_vd::opt::verbose) {
        std::cerr << "Skipping map " << query_map.name_ << " with " 
                  << query_map.frags_.size() << " fragments.\n";
      }

      continue;
    }

    if(query_map.frags_.size() > maligner_vd::opt::max_query_frags) {

      if(maligner_vd::opt::verbose) {
        std::cerr << "Skipping map " << query_map.name_ << " with " 
                  << query_map.frags_.size() << " fragments.\n";
      }

      continue;
    }

    const QueryMapWrapper qmw(query_map, align_opts.query_max_misses);
    const size_t num_query_frags = query_map.frags_.size();

    Timer query_timer;
    query_timer.start();

    ref_score_matrix_db.aln_to_forward_refs(qmw, align_opts);
    ref_score_matrix_db.aln_to_reverse_refs(qmw, align_opts);


    ref_score_matrix_db.compute_query_prefix_mscores(maligner_vd::opt::max_alignments_mad, maligner_vd::opt::min_mad, qmw.get_name());
    ref_score_matrix_db.compute_query_suffix_mscores(maligner_vd::opt::max_alignments_mad, maligner_vd::opt::min_mad, qmw.get_name());

    // std::cerr << "-------------------------------\n";
    // std::cerr << "FORWARD PROFILE:\n";

    // {
    //   maligner_vd::ScoreMatrixProfile profile = ref_score_matrix_db.get_score_matrix_profile_rf_qf(qmw.get_name());
    //   for (auto& rec : profile) {
    //     std::cerr << rec << "\n";
    //   }
    // }

    std::cout << ref_score_matrix_db.get_score_matrix_profile_rf_qf(qmw.get_name())
              << ref_score_matrix_db.get_score_matrix_profile_rf_qr(qmw.get_name())
              << ref_score_matrix_db.get_score_matrix_profile_rr_qf(qmw.get_name())
              << ref_score_matrix_db.get_score_matrix_profile_rr_qr(qmw.get_name());


    std::cerr << "-------------------------------\n";

    query_timer.end();    

    std::cerr << "done aligning query: " << query_map.name_ << " "
              << query_timer << "\n";

    std::cerr << "*****************************************\n";
  
 }

  std::cerr << "maligner_vd done.\n";

  // fout_rf_qf.close();
  // fout_rf_qr.close();
  // fout_rr_qf.close();
  // fout_rr_qr.close();

  return EXIT_SUCCESS;

}
