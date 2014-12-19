#include <iostream>
#include <algorithm>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <chrono>
#include <getopt.h>

#include "map.h"
#include "map_reader.h"
#include "map_frag.h"
#include "map_frag_db.h"
#include "error_model.h"
#include "map_chunk.h"
#include "ref_alignment.h"
#include "map_chunk_db.h"


using namespace std;

//
// Getopt
//
#define PACKAGE_NAME "maligner"
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
            case 'u': arg >> opt::max_unmatched_sites; break;
            case 'r': arg >> opt::rel_error; break;
            case 'a': arg >> opt::min_abs_error; break;
            case 'm': arg >> opt::max_match; break;
            case OPT_MINFRAG: arg >> opt::min_frag; break;
            case OPT_MAX_UNMATCHED_RATE: arg >> opt::max_unmatched_rate; break;
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

    if(opt::max_unmatched_sites < 0)
    {
        std::cerr << "The number of unmatched sites must be greater than or equal to zero.\n";
        die = true;
    }

    if(opt::rel_error < 0) {
      std::cerr << "The relative error cannot be negative.\n";
      die = true;
    }

    if(opt::min_abs_error < 0) {
      std::cerr << "The minimum absolute sizing error cannot be negative\n";
      die = true;
    }

    if(opt::max_match < 1) {
      std::cerr << "The max matches cannot be less than 1\n";
      die = true;
    }

    if(opt::max_unmatched_rate > 1.0 || opt::max_unmatched_rate < 0.0) {
      std::cerr << "Invalid max unmatched rate: " << opt::max_unmatched_rate << "\n";
      die = true;
    }

    if (die) 
    {
        std::cout << "\n" << USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }

    // Parse the query maps file and reference maps file
    opt::query_maps_file = argv[optind++];
    opt::ref_maps_file = argv[optind++];

}


int main(int argc, char* argv[]) {

  using namespace kmer_match;

  opt::program_name = argv[0];
  parse_args(argc, argv);

  cerr << "................................................\n"
       << "MALIGNER settings:\n"
       << "\tquery_maps_file: " << opt::query_maps_file << "\n"
       << "\tref_maps_file: " << opt::ref_maps_file << "\n"
       << "\tmax. consecutive unmatched sites: " << opt::max_unmatched_sites << "\n"
       << "\tmax. unmatched rate: " << opt::max_unmatched_rate << "\n"
       << "\trelative_error: " << opt::rel_error << "\n"
       << "\tabsolute_error: " << opt::min_abs_error << "\n"
       << "\tminimum query frags: " << opt::min_frag << "\n"
       << "\tmax matches per query: " << opt::max_match << "\n"
       << "................................................\n\n";

  ErrorModel error_model(opt::rel_error, opt::min_abs_error);

  auto start_time = chrono::steady_clock::now();

  MapVec ref_maps(read_maps(opt::ref_maps_file));
  MapVec query_maps(read_maps(opt::query_maps_file));

  cerr << "Read " << ref_maps.size() << " reference maps.\n";
  cerr << "Read " << query_maps.size() << " query maps.\n";

  MapChunkDB chunkDB(ref_maps, opt::max_unmatched_sites);
  cerr << "Made MapChunkDB with " << chunkDB.map_chunks_.size() << " chunks.\n";

  size_t aln_count = 0;
  size_t map_with_aln = 0;
  size_t counter = 0;
  size_t interval_report = 100;
  for(auto qi = query_maps.begin(); qi != query_maps.end(); qi++, counter++) {

    if(counter == interval_report) {
      std::cerr << ".";
      counter = 0;
    }

    if(int(qi->frags_.size()) < opt::min_frag) {
      continue;
    }


    const IntVec& frags = qi->get_frags();
    IntPairVec bounds = error_model.compute_bounds(frags.begin() + 1, frags.end() - 1);

    // Determine the maximum number of unmatched sites allowed
    const double mur = opt::max_unmatched_rate;
    size_t max_unmatched = opt::max_unmatched_rate > 0.0 ?
      size_t((mur/(1-mur)) * bounds.size()) : std::numeric_limits<size_t>::max();


    RefAlignmentVec alns = chunkDB.get_compatible_alignments_best(bounds, max_unmatched);
    aln_count += alns.size();

    // std::cerr << qi->name_ << "\n";
    // std::cerr << frags << "\n";
    // std::cerr << "max_unmatched: " << max_unmatched << "\n";
    // std::cerr << "num_frags: " << bounds.size() << "\n";
    // std::cerr << "num aln: " << " " << alns.size() << std::endl;

    if(!alns.empty()) {
      ++map_with_aln;
    }

    // Sort the alignments in ascending order of miss_rate
    std::sort(alns.begin(), alns.end(), ReferenceAlignmentMissRateSort());
    int num_written = 0;
    for(RefAlignmentVec::const_iterator ai = alns.begin();
        ai != alns.end();
        ai++, num_written++) {

       if(ai->miss_rate_ > opt::max_unmatched_rate) continue;
       if(num_written >= opt::max_match) break;

       cout << qi->name_ << " " << *ai << "\n";
    }
  }
  
  cerr << "\ntotal alignments: " << aln_count << "\n";
  cerr << map_with_aln << " maps with alignments.\n";

  auto end_time = chrono::steady_clock::now();
  cerr << chrono::duration <double, milli> (end_time - start_time).count() << " ms" << endl;

  return EXIT_SUCCESS;

}
