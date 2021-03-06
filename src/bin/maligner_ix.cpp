#include <iostream>
#include <algorithm>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <chrono>
#include <getopt.h>
#include <unordered_map>
#include <stdexcept>
#include <limits>

#include "map.h"
#include "map_reader.h"
#include "map_frag.h"
#include "map_frag_db.h"
#include "error_model.h"
#include "map_chunk.h"
#include "ref_alignment.h"
#include "map_chunk_db.h"
#include "matched_chunk.h"
#include "map_wrapper_base.h"
#include "alignment.h"
#include "scorer.h"
#include "bitcover.h"

#include "common_defs.h"

using namespace std;
using namespace maligner_maps;
using namespace kmer_match;

using maligner_dp::Alignment;
using maligner_dp::AlignmentVec;
using maligner_common::Scorer;

//
// Getopt
//
#define PACKAGE_NAME "maligner_ix"
#define SUBPROGRAM ""


static const char *VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by " AUTHOR " " AUTHOR_EMAIL "\n\n";

static const int NUM_POSITION_ARGS = 2;

static const char *USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... QUERY_MAPS_FILE REFERENCE_MAPS_FILE\n"
"Align the maps in the QUERY_MAPS_FILE to the maps in the REFERENCE_MAPS_FILE\n"
"using an index constructed on the reference.\n"
"\n"
"      -h, --help                       display this help and exit\n"
"      -v, --version                    display the version and exit\n"
"      -u, --unmatched VAL              Maximum number of consecutive unmatched sites in reference\n"
"      --rel-error VAL                  Maximum allowed relative sizing error. Default: 0.05\n"
"      -a, --abs-error VAL              Minimum sizing error setting. Default: 1000 bp\n"
"      -m, --max-match VAL              Maximum number of matches per query. Default: 100\n"
"      --min-frag VAL                   Skip query maps with less than min-frag interior fragments. Default: 5\n"
"      --ref-is-circular                Reference map(s) are circular. Default: false\n"
"      --max-unmatched-rate  VAL        Maximum unmatched site rate of an alignment. Default: 0.50\n"
"      --max-score-per-inner-chunk VAL  Maximum score per inner chunk for alignment to be reported. Default: inf\n"
"\n\n"
"Scoring Function Arguments:\n"
"      -q,--query-miss-penalty          Query unmatched site penalty. Default: 18.0\n"
"      -r,--ref-miss-penalty            Reference unmatched site penalty. Default: 3.0\n"
"      --sd-rate                        Standard deviation rate. Default: 0.05\n"
"      --min-sd                         Minimum standard deviation (bp) Default: 500\n"
;


namespace opt
{
    static int max_unmatched_sites = 2;
    static double rel_error = 0.05;
    static int min_abs_error = 1000;
    static int min_frag = 5;
    static int max_match = 100;
    static double max_unmatched_rate = 0.50;
    static bool ref_is_circular = false;
    static bool ref_is_bounded = false;
    static double query_miss_penalty = 18.0;
    static double ref_miss_penalty = 3.0;
    static double sd_rate = 0.05;
    static double min_sd = 500.0;
    static double max_score_per_inner_chunk = std::numeric_limits<double>::infinity();
    static string query_maps_file;
    static string ref_maps_file;
    string program_name;
}

static const char* shortopts = "u:r:a:m:q:r:hv";
enum {OPT_MIN_FRAG = 1, OPT_MAX_UNMATCHED_RATE, OPT_REF_IS_CIRCULAR, OPT_SD_RATE, OPT_MIN_SD,
  OPT_REL_ERROR, OPT_MAX_SCORE_PER_INNER_CHUNK};

static const struct option longopts[] = {
    { "unmatched", required_argument, NULL, 'u' },
    { "max-unmatched-rate", required_argument, NULL, OPT_MAX_UNMATCHED_RATE},
    { "rel-error", required_argument, NULL, OPT_REL_ERROR },
    { "abs-error", required_argument, NULL, 'a' },
    { "min-frag", required_argument,  NULL, OPT_MIN_FRAG},
    { "max-match", required_argument, NULL, 'm'},
    { "query-miss-penalty", required_argument, NULL, 'q'},
    { "ref-miss-penalty", required_argument, NULL, 'r'},
    { "min-sd", required_argument, NULL, OPT_MIN_SD},
    { "sd-rate", required_argument, NULL, OPT_SD_RATE},
    { "ref-is-circular", no_argument, NULL, OPT_REF_IS_CIRCULAR},
    { "max-score-per-inner-chunk", required_argument, NULL, OPT_MAX_SCORE_PER_INNER_CHUNK},
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
            case OPT_REL_ERROR: arg >> opt::rel_error; break;
            case 'a': arg >> opt::min_abs_error; break;
            case 'm': arg >> opt::max_match; break;
            case 'q': arg >> opt::query_miss_penalty; break;
            case 'r': arg >> opt::ref_miss_penalty; break;
            case OPT_SD_RATE: arg >> opt::sd_rate; break;
            case OPT_MIN_SD: arg >> opt::min_sd; break;
            case OPT_MIN_FRAG: arg >> opt::min_frag; break;
            case OPT_MAX_UNMATCHED_RATE: arg >> opt::max_unmatched_rate; break;
            case OPT_MAX_SCORE_PER_INNER_CHUNK: arg >> opt::max_score_per_inner_chunk; break;
            case OPT_REF_IS_CIRCULAR:
              opt::ref_is_circular = true;
              opt::ref_is_bounded = true;
              break;
            case 'h':
            {
                std::cerr << USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
                break;
            }
            case 'v':
            {
                std::cerr << VERSION_MESSAGE;
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
        std::cerr << "\n" << USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }

    // Parse the query maps file and reference maps file
    opt::query_maps_file = argv[optind++];
    opt::ref_maps_file = argv[optind++];

}

typedef std::unordered_map<std::string, MapWrapper> MapDB;

// Select the non overlapping set of alignments. Alignment should already be sorted in the desired order.
maligner_dp::AlignmentVec sift_alignments(maligner_dp::AlignmentVec& alns, MapDB& ref_map_db) {

  maligner_dp::AlignmentVec alns_out;
  alns_out.reserve(alns.size());

  // Reset the bitcovers of the reference maps.
  for(auto i = ref_map_db.begin(); i != ref_map_db.end(); i++) {
    MapWrapper& ref_map = i->second;
    ref_map.bit_cover_.reset();
  }

  for(auto i = alns.begin(); i != alns.end(); i++) {

    maligner_dp::Alignment aln = *i;
    MapWrapper& ref_map = ref_map_db.at(aln.ref_map_data.map_name_);

    if(ref_map.bit_cover_.is_covered(aln.ref_start, aln.ref_end)) {
      continue;
    }

    // Make sure that the alignment does not start past the point of doubling
    // for circular maps
    if(aln.ref_start > ref_map.num_frags()) {
      continue;
    }

    ref_map.bit_cover_.cover(aln.ref_start, aln.ref_end);
    alns_out.push_back(std::move(aln));

  }

  return alns_out;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////
// Convert alignments to the common format used by maligner dp.
// maligner_ix handles reverse alignments by reversing the reference (since we index the reverse of the reference)
// maligner_dp handles reverse alignments by reversing the query and aligning it to the forward of the reference.
// Therefore, we must carefully convert indices and orientation in this function to end up with a maligner_dp::Alignment.
maligner_dp::AlignmentVec convert_alignments(
  const RefAlignmentVec& ref_alignments,
  const MapWrapper& query,
  MapDB& ref_map_db,
  const Scorer& scorer) {

  using maligner_dp::Score;
  using maligner_dp::MatchedChunk;
  using maligner_dp::MatchedChunkVec;
  using maligner_dp::Chunk;
  using maligner_dp::ChunkVec;
  using maligner_dp::MapData;
  using maligner_dp::AlignmentRescaledScoreComp;

  maligner_dp::AlignmentVec alns;

  const MapData& query_map_data = query.map_data_;
  const Map& query_map = query.map_;

  bool is_circular = false;
  bool is_bounded = false;

  Score zero_score;

  //////////////////////////////////////////////////////////////////////
  // Our approach is to represent the query_chunks the same way that the query map would
  // be presented in maligner_dp for forward and reverse alignments.
  // This means that query_chunks_forward and query_chunks_reverse both must have increasing indices,
  // but query_chunks_reverse has the chunk sizes in the *reverse* direction.
  ChunkVec query_chunks_forward;
  ChunkVec query_chunks_reverse;
  const size_t num_query_frags = query_map.frags_.size();
  for(size_t i = 1; i < num_query_frags - 1; i++) {

    const bool is_boundary_chunk_query = (i == 0) || (i == num_query_frags - 1);

    query_chunks_forward.emplace_back(i, i + 1, query_map.frags_[i], is_boundary_chunk_query);

    query_chunks_reverse.emplace_back(i, i + 1, query_map.frags_[i], is_boundary_chunk_query);
    query_chunks_reverse.back().reverse_coords(num_query_frags);

  }
  std::reverse(query_chunks_reverse.begin(), query_chunks_reverse.end());

  
  ///////////////////////////////////////////////////////////////////////////////////
  // In maligner_ix, ref alignments are oriented forward with respect to the query.
  // In maligner_dp, MatchedChunks must be constructed, oriented with respect to forward strand of reference.
  // This means if the alignment is reverse, we need to reverse the reference chunks,
  // and use the reverse representation of the query chunks.
  for(auto i = ref_alignments.begin(); i != ref_alignments.end(); i++) {
    
    const ReferenceAlignment& ref_alignment = *i;
    const Map * p_ref_map = ref_alignment.get_map();
    const MapWrapper& ref = ref_map_db.find(p_ref_map->name_)->second;
    const MapData& ref_map_data = ref.map_data_;
    const bool ref_is_circular = ref.is_circular();

    /////////////////////////////////////////////////////////////////////////////////////////
    // Extract the reference chunks as maligner_dp chunks.
    ChunkVec ref_chunks;
    for(auto rc = ref_alignment.chunks_.begin(); rc != ref_alignment.chunks_.end(); rc++) {

      const MapChunk* p_chunk = *rc;

      bool ref_chunk_is_boundary = ref_is_circular && !opt::ref_is_bounded && 
        ( (p_chunk->start_ == 0) || (p_chunk->end_ == ref.num_frags()) );

      ref_chunks.emplace_back(p_chunk->start_, p_chunk->end_, p_chunk->size_, ref_chunk_is_boundary);

    }

    // Orient the reference chunks forward if necessary
    ChunkVec * p_query_chunks = &query_chunks_forward;
    if(ref_alignment.is_reverse()) {
      std::reverse(ref_chunks.begin(), ref_chunks.end());
      p_query_chunks = &query_chunks_reverse;
    }


    // std::cerr << "query:\n\t" << *p_query_chunks << "\n"
    //           << "ref:\n\t" << ref_chunks << "\n"
    //           << "is_forward: " << ref_alignment.is_forward()
    //           << "\n";

    /////////////////////////
    // Build matched chunks
    const ChunkVec& query_chunks = *p_query_chunks;
    if (query_chunks.size() != ref_chunks.size()) {
      std::cerr << "query_frags: " << query_map.frags_.size() << " query_chunks: " << query_chunks.size() << " ref_chunks: " << ref_chunks.size() << std::endl;
      throw std::runtime_error("query chunks does not match ref_chunks size.");
    }

    MatchedChunkVec matched_chunks;
    Score total_score;
    const size_t num_matched_chunks = query_chunks.size();
    for(size_t i = 0; i < num_matched_chunks; i++) {
      Score score = scorer.compute_score(query_chunks[i], ref_chunks[i]);
      total_score += score;
      matched_chunks.emplace_back(query_chunks[i], ref_chunks[i], score);
    }

    /////////////////////////
    // Construct Alignment
    bool aln_is_forward = ref_alignment.is_forward();
    maligner_dp::Alignment aln(matched_chunks, total_score, query_map_data, ref_map_data, aln_is_forward);
    if (!aln_is_forward) aln.flip_query_coords();
    aln.add_alignment_locs(query.ix_to_locs_, ref.ix_to_locs_);
    alns.push_back(std::move(aln));

  }

  // Sort the alignments by total rescaled score
  std::sort(alns.begin(), alns.end(), AlignmentRescaledScoreComp());

  alns = sift_alignments(alns, ref_map_db);

  return alns;

}


/////////////////////////////////////////
// Read reference maps.
MapWrapperVec read_ref_maps(std::string& maps_file, bool is_circular, bool is_bounded) {
  
  MapReader ref_map_reader(maps_file);
  MapVec ref_maps_temp = ref_map_reader.read_all_maps();
  MapWrapperVec ref_maps;

  for(size_t i = 0; i < ref_maps_temp.size(); i++) {
    Map& ref_map = ref_maps_temp[i];
    ref_maps.emplace_back(std::move(ref_map), is_circular, is_bounded);
  }

  return ref_maps;
}


int main(int argc, char* argv[]) {

  using namespace kmer_match;
  using maligner_dp::AlignmentHeader;

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

  /////////////////////////////////////////////////////////////
  // Read Reference Maps

  MapWrapperVec ref_maps = read_ref_maps(opt::ref_maps_file, opt::ref_is_circular, opt::ref_is_bounded);
 
  // Store a vector of pointers to the underlying Maps.
  MapWrapperPVec p_ref_maps(ref_maps.size());
  for(size_t i = 0; i < ref_maps.size(); i++) {
    p_ref_maps[i] = &ref_maps[i];
  }

  // Store reference maps in an unordered map.
  MapDB ref_map_db;
  for(auto i = ref_maps.begin(); i != ref_maps.end(); i++) {
    ref_map_db.insert( MapDB::value_type(i->map_.name_, *i) );
  }

  cerr << "Read " << ref_maps.size() << " reference maps.\n";

  ////////////////////////////////////////////////////////////////
  // Build Chunk Database
  MapChunkDB chunkDB(p_ref_maps, opt::max_unmatched_sites);
  cerr << "Made MapChunkDB with " << chunkDB.map_chunks_.size() << " chunks.\n";

  // Write alignment header
  std::cout << AlignmentHeader();

  size_t aln_count = 0;
  size_t map_with_aln = 0;
  size_t counter = 0;
  size_t interval_report = 100;
  

  Scorer scorer(opt::query_miss_penalty, opt::ref_miss_penalty, opt::min_sd, opt::sd_rate);

  // Iterate through query maps file.
  Map query_map;
  MapReader query_map_reader(opt::query_maps_file);
  for(; query_map_reader.next(query_map); counter++) {

    if(counter == interval_report) {
      std::cerr << ".";
      counter = 0;
    }

    MapWrapper query = MapWrapper(std::move(query_map), false, false);

    if(int(query.num_frags()) < opt::min_frag) {
      continue;
    }

    std::cerr << "Aligning " << query.get_name() << "\n";

    const FragVec& frags = query.get_frags();

    // Note: We ignore boundary fragments
    IntPairVec bounds = error_model.compute_bounds(frags.begin() + 1, frags.end() - 1);

    // Determine the maximum number of unmatched sites allowed
    const double mur = opt::max_unmatched_rate;
    size_t max_unmatched = mur > 0.0 ?
      size_t((mur/(1-mur)) * bounds.size()) : std::numeric_limits<size_t>::max();


    RefAlignmentVec ref_alns = chunkDB.get_compatible_alignments_best(bounds, max_unmatched);

    // std::cerr << qi->name_ << "\n";
    // std::cerr << frags << "\n";
    // std::cerr << "max_unmatched: " << max_unmatched << "\n";
    // std::cerr << "num_frags: " << bounds.size() << "\n";
    // std::cerr << "num aln: " << " " << ref_alns.size() << std::endl;


    // Sort the alignments in ascending order of miss_rate
    std::sort(ref_alns.begin(), ref_alns.end(), ReferenceAlignmentMissRateSort());
    AlignmentVec alns = convert_alignments(ref_alns, query, ref_map_db, scorer);
    cerr << "Converted " << alns.size() << " alignments.\n";
    int query_aln_count = 0;
    bool map_has_alignment = false;
    for(auto a = alns.begin(); a != alns.end(); a++) {

      if(query_aln_count == opt::max_match) break;

      const maligner_dp::Alignment& aln = *a;
      bool have_aln = false;

      if(aln.score_per_inner_chunk <= opt::max_score_per_inner_chunk) {
        maligner_dp::print_alignment(std::cout, aln);
        ++aln_count;
        ++query_aln_count;
        have_aln = true;
        map_has_alignment = true;
      }

    }

    if(map_has_alignment) {
      ++map_with_aln;
    }

  }
  
  cerr << "\ntotal alignments: " << aln_count << "\n";
  cerr << map_with_aln << " maps with alignments.\n";

  auto end_time = chrono::steady_clock::now();
  cerr << chrono::duration <double, milli> (end_time - start_time).count() << " ms" << endl;


  std::cerr << "maligner done.\n";
  return EXIT_SUCCESS;

}
