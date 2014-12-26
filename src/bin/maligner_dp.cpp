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

using std::string;
using std::unordered_map;
// using namespace std;

//
// Getopt
//
#define PACKAGE_NAME "maligner dp"
#define PACKAGE_VERSION "0.2"
#define SUBPROGRAM ""


static const char *VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Lee Mendelowitz. (lmendelo@umiacs.umd.edu) \n"
"\n";

static const int NUM_POSITION_ARGS = 2;

static const char *USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... QUERY_MAPS_FILE REFERENCE_MAPS_FILE\n"
"Align the maps in the QUERY_MAPS_FILE to the maps in the REFERENCE_MAPS_FILE using dynamic programming.\n"
"\n"
"      -h, --help                       display this help and exit\n"
"      -v, --version                    display the version and exit\n"
"      -u, --unmatched VAL              Maximum number of consecutive unmatched sites in reference\n"
"      -r, --rel_error VAL              Maximum allowed relative sizing error. Default: 0.05\n"
"      -a, --abs_error VAL              Minimum sizing error setting. Default: 1000 bp\n"
"      -m, --max_match VAL              Maximum number of matches per query. Default: 100\n"
"      --min_frag VAL                   Skip query maps with less than min-frag interior fragments. Default: 5\n"
"      --max_unmatched_rate  VAL        Maximum unmatched site rate of an alignment. Default: 0.50\n"
;

namespace maligner_dp {
  namespace opt
  {

      static int max_unmatched_sites = 2;
      static double rel_error = 0.05;
      static int min_abs_error = 1000;
      static int min_frag = 5;
      static int max_match = 100;
      static double max_unmatched_rate = 0.50;
      static string query_maps_file;
      static string ref_maps_file;
      string program_name;

      static double query_miss_penalty = 30.0;
      static double ref_miss_penalty = 3.0;
      static int query_max_misses = 2;
      static int ref_max_misses = 5;
      static double sd_rate = 0.10;
      static double min_sd = 500.0;
      static double max_chunk_sizing_error = 4.0;
      static int alignments_per_reference = 1;
      static int min_alignment_spacing = 10;
      static int neighbor_delta = 0;
      static bool query_is_bounded = false;
      static bool ref_is_bounded = false;

  }
}

static const char* shortopts = "u:r:a:m:hv";
enum {OPT_MINFRAG = 1, OPT_MAX_UNMATCHED_RATE};

static const struct option longopts[] = {
    { "unmatched", required_argument, NULL, 'u' },
    { "max_unmatched_rate", required_argument, NULL, OPT_MAX_UNMATCHED_RATE},
    { "rel_error", required_argument, NULL, 'r' },
    { "abs_error", required_argument, NULL, 'a' },
    { "min_frag", required_argument,  NULL, OPT_MINFRAG},
    { "max_match", required_argument, NULL, 'm'},
    { "help",     no_argument,       NULL, 'h' },
    { "version",  no_argument,       NULL, 'v'},
    { NULL, 0, NULL, 0 }
};

void parse_args(int argc, char** argv)
{

    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) 
        {
            case 'u': arg >> maligner_dp::opt::max_unmatched_sites; break;
            case 'r': arg >> maligner_dp::opt::rel_error; break;
            case 'a': arg >> maligner_dp::opt::min_abs_error; break;
            case 'm': arg >> maligner_dp::opt::max_match; break;
            case OPT_MINFRAG: arg >> maligner_dp::opt::min_frag; break;
            case OPT_MAX_UNMATCHED_RATE: arg >> maligner_dp::opt::max_unmatched_rate; break;
            case 'h':
            {
                std::cout << USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
                break;
            }
            case 'v':
            {
                std::cout << VERSION_MESSAGE;
                exit(EXIT_SUCCESS);
                break;
            }
        }
    }

    if (argc - optind < NUM_POSITION_ARGS) 
    {
        std::cerr << SUBPROGRAM ": missing arguments\n";
        die = true;
    } 
    else if (argc - optind > NUM_POSITION_ARGS) 
    {
        std::cerr << SUBPROGRAM ": too many arguments\n";
        die = true;
    }

    if(maligner_dp::opt::max_unmatched_sites < 0)
    {
        std::cerr << "The number of unmatched sites must be greater than or equal to zero.\n";
        die = true;
    }

    if(maligner_dp::opt::rel_error < 0) {
      std::cerr << "The relative error cannot be negative.\n";
      die = true;
    }

    if(maligner_dp::opt::min_abs_error < 0) {
      std::cerr << "The minimum absolute sizing error cannot be negative\n";
      die = true;
    }

    if(maligner_dp::opt::max_match < 1) {
      std::cerr << "The max matches cannot be less than 1\n";
      die = true;
    }

    if(maligner_dp::opt::max_unmatched_rate > 1.0 || maligner_dp::opt::max_unmatched_rate < 0.0) {
      std::cerr << "Invalid max unmatched rate: " << maligner_dp::opt::max_unmatched_rate << "\n";
      die = true;
    }

    if (die) 
    {
        std::cout << "\n" << USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }

    // Parse the query maps file and reference maps file
    maligner_dp::opt::query_maps_file = argv[optind++];
    maligner_dp::opt::ref_maps_file = argv[optind++];

}


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
 ScoreMatrix sm;
 MapReader query_map_reader(maligner_dp::opt::query_maps_file);
 Map query_map;
 AlignmentVec alns_forward, alns_reverse;


 std::cout << AlignmentHeader();
 while(query_map_reader.next(query_map)) {


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
      sm.resize(num_query_frags + 1, num_ref_frags + 1);

      // const IntVec* p_frags_forward = &query_frags_forward;
      // const IntVec* p_frags_reverse = &query_frags_reverse;

      AlignTask task_forward(&qmw.md_, &rmw.md_,
        &query_frags_forward, &rmw.m_.frags_, 
        &qps_forward, &rmw.ps_,
        0, // ref_offset
        &sm, &alns_forward, align_opts
      );

      AlignTask task_reverse(&qmw.md_, &rmw.md_,
        &query_frags_reverse, &rmw.m_.frags_, 
        &qps_reverse, &rmw.ps_,
        0, // ref_offset
        &sm, &alns_reverse, align_opts
      );

      // std::cerr << "Aligning " << query_map.name_ << " to " << rmw.m_.name_ << "\n";
      Timer timer;

      // Align Forward
      timer.start();
      Alignment aln_forward = make_best_alignment_using_partials(task_forward);
      timer.end();
      // std::cerr << "Alignments forward is valid: " << aln_forward.is_valid
      //           << " " << timer << "\n";

      if(aln_forward.is_valid) {
        // alns_forward.push_back(std::move(aln_forward));
        alns_forward.push_back(aln_forward);
        // std::cout << aln_forward << "\n";
      }

      // Align Reverse
      // TO DO: We need to have a means of building reverse alignments with the 
      // correct labeling of indices. Right now we are ignoring this.
      timer.start();
      Alignment aln_reverse = make_best_alignment_using_partials(task_reverse);
      timer.end();
      // std::cerr << "Alignments reverse is valid: " << aln_reverse.is_valid
      //           << " " << timer << "\n";

      if(aln_reverse.is_valid) {
        // alns_reverse.push_back(std::move(aln_reverse));
        alns_reverse.push_back(aln_reverse);
        // std::cout << aln_reverse << "\n";
      }

    }

    // Of all the forward alignments and reverse alignments, print only the best.
    // cerr << "forward scores:\n";
    // for(auto& v: alns_forward) {
    //   cerr << "aln score: " << v.total_score << "\n";
    // }
    // cerr << "reverse scores:\n";
    // for(auto& v: alns_reverse) {
    //   cerr << "aln score: " << v.total_score << "\n";
    // }

    std::sort(alns_forward.begin(), alns_forward.end(), AlignmentScoreComp());

    // cout << "Best forward: ";
    if( alns_forward.empty() ) {
      // cout << "None\n";
    } else {
      print_alignment(std::cout, alns_forward.front(), true);

      // cout << alns_forward[0] << "\n";
    }

    std::sort(alns_reverse.begin(), alns_reverse.end(), AlignmentScoreComp());

    // cout << "Best reverse: ";
    if ( alns_reverse.empty() ) {
      // cout << "None\n";
    } else {
      print_alignment(std::cout, alns_reverse.front(), false);
      // cout << alns_reverse[0] << "\n";
    }
    
    query_timer.end();    

    std::cerr << "Aligned query " << query_map.name_ << "to all references. "
              << query_timer << "\n";
    std::cerr << "*****************************************\n";
 }
 

  return EXIT_SUCCESS;

}
