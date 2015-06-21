//
// Getopt
//
#include <limits>



static const char *VERSION_MESSAGE = "Version " PACKAGE_VERSION "\n"
"Written by " AUTHOR "(" AUTHOR_EMAIL ") \n"
"\n";

static const int NUM_POSITION_ARGS = 2;

static const char *USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " [OPTION] ... QUERY_MAPS_FILE REFERENCE_MAPS_FILE\n"
"\n"
" Align the maps in the QUERY_MAPS_FILE to the maps in the REFERENCE_MAPS_FILE\n"
" using dynamic programming.\n"
"\n"
" Alignment options:\n"
"      --reference-is-circular          Treat reference maps as circular. Default: false\n"
"      --num-permutation-trials         Number of trials for the permutation test.\n"
"                                           Default: 0\n"
"      --no-query-rescaling             Default: perform query rescaling\n"
"\n"
" Scoring parameters:\n"
"      -q,--query-miss-penalty          Query unmatched site penalty\n"
"      -r,--ref-miss-penalty            Reference unmatched site penalty\n"
"      --sd-rate                        Standard deviation rate.\n"
"      --min-sd                         Minimum standard deviation (bp)\n"
"\n"
" Search space parameters:\n"
"      --query-max-misses               Query max. consecutive unmatched sites\n"
"      --ref-max-misses                 Reference max. consecutive unmatched sites\n"
"      --query-max-miss-rate            Max. rate of unmatched sites in the query\n"
"      --ref-max-miss-rate              Max. rate of unmatched sites in the reference\n"
"      --max-chunk-sizing-error         Max. chunk sizing error score for bounding search\n"
"                                          space. Default: Inf\n"
"      --max-alignments-per-reference   Max. alignments to report per reference\n"
"      --max-alignments                 Max. number of alignments to output\n"
"\n"
" Alignment filters:\n"
"      --max-score-per-inner-chunk      Report alignments with a score per inner chunk less than this\n"
"                                           threshold. Default: Inf\n"
"\n"
" General arguments:\n"
"      -h, --help                       display this help and exit\n"
"      -v, --version                    display the version and exit\n"
"      --verbose                        Verbose output\n";



namespace maligner_dp {
  namespace opt
  {

      static string query_maps_file;
      static string ref_maps_file;
      string program_name;

      static double query_miss_penalty = 18.0;
      static double ref_miss_penalty = 3.0;
      static int query_max_misses = 2;
      static int ref_max_misses = 5;
      static double sd_rate = 0.05;
      static double min_sd = 500.0;
      static double max_chunk_sizing_error = std::numeric_limits<double>::infinity();
      static double ref_max_miss_rate = 0.50;
      static double query_max_miss_rate = 0.25;
      static int alignments_per_reference = 100;
      static int min_alignment_spacing = 10;
      static int neighbor_delta = 0;
      static bool query_is_bounded = false;
      static bool ref_is_bounded = false;
      static int max_alignments = 100;
      static int max_alignments_mad = 100; // Max alignments to use for mad computation
      static double max_score_per_inner_chunk = std::numeric_limits<double>::infinity();
      static int num_permutation_trials = 0; // Number of trials for permutation test.
      static bool query_rescaling = true;
      static bool verbose = false;
      static bool reference_is_circular = false;
      static double min_mad = 1.0; // Minimum mad to use when computing mad scores.

  }
}

static const char* shortopts = "q:r:hv";
enum {
  OPT_QUERY_MAX_MISSES = 1,
  OPT_REF_MAX_MISSES,
  OPT_SD_RATE,
  OPT_MIN_SD,
  OPT_MAX_CHUNK_SIZING_ERROR,
  OPT_ALIGNMENTS_PER_REFERENCE,
  OPT_MAX_ALIGNMENTS,
  OPT_NUM_PERMUTATION_TRIALS,
  OPT_REF_MAX_MISS_RATE,
  OPT_QUERY_MAX_MISS_RATE,
  OPT_MAX_SCORE_PER_INNER_CHUNK,
  OPT_VERBOSE,
  OPT_NO_QUERY_RESCALING,
  OPT_REFERENCE_IS_CIRCULAR
};

static const struct option longopts[] = {

    { "query-miss-penalty", required_argument, NULL, 'q'},
    { "ref-miss-penalty", required_argument, NULL, 'r'},
    { "query-max-misses", required_argument, NULL, OPT_QUERY_MAX_MISSES},
    { "ref-max-misses", required_argument, NULL, OPT_REF_MAX_MISSES},
    { "query-max-miss-rate", required_argument, NULL, OPT_QUERY_MAX_MISS_RATE},
    { "ref-max-miss-rate", required_argument, NULL, OPT_REF_MAX_MISS_RATE},
    { "sd-rate", required_argument, NULL, OPT_SD_RATE},
    { "min-sd", required_argument, NULL, OPT_MIN_SD},
    { "max-chunk-sizing-error", required_argument, NULL, OPT_MAX_CHUNK_SIZING_ERROR},
    { "max-alignments-per-reference", required_argument, NULL, OPT_ALIGNMENTS_PER_REFERENCE},
    { "max-alignments", required_argument, NULL, OPT_MAX_ALIGNMENTS},
    { "max-score-per-inner-chunk", required_argument, NULL, OPT_MAX_SCORE_PER_INNER_CHUNK},
    { "num-permutation-trials", required_argument, NULL, OPT_NUM_PERMUTATION_TRIALS},
    { "no-query-rescaling", no_argument, NULL, OPT_NO_QUERY_RESCALING},
    { "reference-is-circular", no_argument, NULL, OPT_REFERENCE_IS_CIRCULAR},
    { "verbose", no_argument, NULL, OPT_VERBOSE},
    { "help",     no_argument,       NULL, 'h' },
    { "version",  no_argument,       NULL, 'v'},
    { NULL, 0, NULL, 0 }

};

void parse_args(int argc, char** argv)
{
    using namespace maligner_dp;

    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) 
        {
            case 'q': arg >> opt::query_miss_penalty; break;
            case 'r': arg >> opt::ref_miss_penalty; break;
            case OPT_QUERY_MAX_MISSES: arg >> opt::query_max_misses; break;
            case OPT_REF_MAX_MISSES: arg >> opt::ref_max_misses; break;
            case OPT_QUERY_MAX_MISS_RATE: arg >> opt::query_max_miss_rate; break;
            case OPT_REF_MAX_MISS_RATE: arg >> opt::ref_max_miss_rate; break;
            case OPT_SD_RATE: arg >> opt::sd_rate; break;
            case OPT_MIN_SD: arg >> opt::min_sd; break;
            case OPT_MAX_CHUNK_SIZING_ERROR: arg >> opt::max_chunk_sizing_error; break;
            case OPT_MAX_SCORE_PER_INNER_CHUNK: arg >> opt::max_score_per_inner_chunk; break;
            case OPT_ALIGNMENTS_PER_REFERENCE: arg >> opt::alignments_per_reference; break;
            case OPT_MAX_ALIGNMENTS: arg >> opt::max_alignments; break;
            case OPT_NUM_PERMUTATION_TRIALS: arg >> opt::num_permutation_trials; break;
            case OPT_VERBOSE: opt::verbose = true; break;
            case OPT_NO_QUERY_RESCALING: opt::query_rescaling = false; break;
            case OPT_REFERENCE_IS_CIRCULAR: 
              opt::reference_is_circular = true;
              opt::ref_is_bounded = true;
              break;
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

        if(arg.fail()) {
          std::cerr << "Trouble parsing input: " << optarg << "\n";
          die = true;
          break;
        }
    }

    if (argc - optind < NUM_POSITION_ARGS) 
    {
        std::cerr << ": missing arguments\n";
        die = true;
    } 
    else if (argc - optind > NUM_POSITION_ARGS) 
    {
        std::cerr << ": too many arguments\n";
        die = true;
    }

    if(opt::query_miss_penalty < 0.0) {
      std::cerr << "Query miss penalty must be non-negative.\n";
      die = true;
    }

    if(opt::ref_miss_penalty < 0.0) {
      std::cerr << "Reference miss penalty must be non-negative.\n";
      die = true;
    }

    if(opt::query_max_miss_rate < 0.0) {
      std::cerr << "Query max miss rate must be non-negative.\n";
      die = true;
    }

    if(opt::ref_max_miss_rate < 0.0) {
      std::cerr << "Reference max miss rate must be non-negative.\n";
      die = true;
    }

    if(opt::query_max_misses < 0) {
      std::cerr << "Query max misses must be positive\n";
      die = true;
    }

    if(opt::ref_max_misses < 0) {
      std::cerr << "Ref max misses must be positive\n";
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


std::ostream& print_args(std::ostream& os) {
  using namespace maligner_dp::opt;

  os << VERSION_MESSAGE << "\n"
     << "Settings:\n"
     << "\tquery_maps_file: " << query_maps_file << "\n"
     << "\tref_maps_file: " << ref_maps_file << "\n"
     << "\tverbose: " << verbose << "\n"
     << "\tquery_miss_penalty: " << query_miss_penalty << "\n"
     << "\tref_miss_penalty: " << ref_miss_penalty << "\n"
     << "\tquery_max_misses: " << query_max_misses << "\n"
     << "\tref_max_misses: " << ref_max_misses << "\n"
     << "\tquery_max_miss_rate: " << query_max_miss_rate << "\n"
     << "\tref_max_miss_rate: " << ref_max_miss_rate << "\n"
     << "\tsd_rate: " << sd_rate << "\n"
     << "\tmin_sd: " << min_sd << "\n"
     << "\tmax_chunk_sizing_error: " << max_chunk_sizing_error << "\n"
     << "\tmax_score_per_inner_chunk: " << max_score_per_inner_chunk << "\n"
     << "\talignments_per_reference: " << alignments_per_reference << "\n"
     << "\tmax_alignments: " << max_alignments << "\n"
     << "\tmin_alignment_spacing: " << min_alignment_spacing << "\n"
     << "\tnum_permutation_trials: " << num_permutation_trials << "\n"
     << "\tneighbor_delta: " << neighbor_delta << "\n"
     << "\treference_is_circular: " << reference_is_circular << "\n"
     << "\tquery_rescaling: " << query_rescaling << "\n"
     << "\tquery_is_bounded: " << query_is_bounded << "\n"
     << "\tref_is_bounded: " << ref_is_bounded << "\n"
     << "\treference_is_circular: " << reference_is_circular << "\n";

  return os;
}