//
// Getopt
//
#define PACKAGE_NAME "maligner dp"

static const char *VERSION_MESSAGE = "Version " PACKAGE_VERSION "\n"
"Written by " AUTHOR "(" AUTHOR_EMAIL ") \n"
"\n";

static const int NUM_POSITION_ARGS = 2;

static const char *USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " [OPTION] ... QUERY_MAPS_FILE REFERENCE_MAPS_FILE\n"
"Align the maps in the QUERY_MAPS_FILE to the maps in the REFERENCE_MAPS_FILE using dynamic programming.\n"
"\n"
"      -h, --help                       display this help and exit\n"
"      -v, --version                    display the version and exit\n"
"      -q,--query-miss-penalty          Query unmatched site penalty\n"
"      -r,--ref-miss-penalty            Reference unmatched site penalty\n"
"      --query-max-misses               Query maximum consecutive unmatched\n"
"      --ref-max-misses                 Reference maximum consecutive unmatched\n"
"      --sd-rate                        Standard deviation rate.\n"
"      --min-sd                         Minimum standard deviation (bp)\n"
"      --max-chunk-sizing-error         \n"
"      --alignments-per-reference       \n"
"      --max-alignments                 Max. number of alignments to output\n";


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
      static double max_chunk_sizing_error = 16.0;
      static int alignments_per_reference = 100;
      static int min_alignment_spacing = 10;
      static int neighbor_delta = 0;
      static bool query_is_bounded = false;
      static bool ref_is_bounded = false;
      static int max_alignment_seeds = 1000; // Max number of cells to seed alignments from for each reference
      static int max_alignments = 100;

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
  OPT_MAX_ALIGNMENT_SEEDS,
  OPT_MAX_ALIGNMENTS
};

static const struct option longopts[] = {

    { "query-miss-penalty", required_argument, NULL, 'q'},
    { "ref-miss-penalty", required_argument, NULL, 'r'},
    { "query-max-misses", required_argument, NULL, OPT_QUERY_MAX_MISSES},
    { "ref-max-misses", required_argument, NULL, OPT_REF_MAX_MISSES},
    { "sd-rate", required_argument, NULL, OPT_SD_RATE},
    { "min-sd", required_argument, NULL, OPT_MIN_SD},
    { "max-chunk-sizing-error", required_argument, NULL, OPT_MAX_CHUNK_SIZING_ERROR},
    { "alignments-per-reference", required_argument, NULL, OPT_ALIGNMENTS_PER_REFERENCE},
    { "max-alignments", required_argument, NULL, OPT_MAX_ALIGNMENTS},
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
            case OPT_SD_RATE: arg >> opt::sd_rate; break;
            case OPT_MIN_SD: arg >> opt::min_sd; break;
            case OPT_MAX_CHUNK_SIZING_ERROR: arg >> opt::max_chunk_sizing_error; break;
            case OPT_ALIGNMENTS_PER_REFERENCE: arg >> opt::alignments_per_reference; break;
            case OPT_MAX_ALIGNMENTS: arg >> opt::max_alignments; break;
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
     << "\tquery_miss_penalty: " << query_miss_penalty << "\n"
     << "\tref_miss_penalty: " << ref_miss_penalty << "\n"
     << "\tquery_max_misses: " << query_max_misses << "\n"
     << "\tref_max_misses: " << ref_max_misses << "\n"
     << "\tsd_rate: " << sd_rate << "\n"
     << "\tmin_sd: " << min_sd << "\n"
     << "\tmax_chunk_sizing_error: " << max_chunk_sizing_error << "\n"
     << "\talignments_per_reference: " << alignments_per_reference << "\n"
     << "\tmax_alignments: " << max_alignments << "\n"
     << "\tmin_alignment_spacing: " << min_alignment_spacing << "\n"
     << "\tneighbor_delta: " << neighbor_delta << "\n"
     << "\tquery_is_bounded: " << query_is_bounded << "\n"
     << "\tref_is_bounded: " << ref_is_bounded << "\n";

  return os;
}