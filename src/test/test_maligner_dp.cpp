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
"Align the maps in the QUERY_MAPS_FILE to the maps in the REFERENCE_MAPS_FILE\n"
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
      static double min_sd = 1000.0;
      static double max_chunk_sizing_error = 16.0;
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
        AlignTask task(&qmw.md_, &rmw.md_,
          &qmw.m_.frags_, &rmw.m_.frags_, 
          &qmw.ps_, &rmw.ps_,
          0,
          &sm, &alns, align_opts
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
