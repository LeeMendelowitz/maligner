#include <vector>
#include <cassert>
#include <iostream>
#include <utility>
#include <algorithm>
#include <set>
#include <queue>

using std::cerr;

#include "utils.h"
#include "types.h"
#include "align.h"
#include "globals.h"
#include "alignment.h"
#include "ScoreMatrix.h"
#include "ScoreCell.h"

#define DEBUG 0
#define GET_BEST_DEBUG 0
#define BUILD_TRAIL_DEBUG 0
#define BREAKS_DEBUG 0
#define RESCALE_DEBUG 0
#define NEIGHBORHOOD_DEBUG 0
#define QUEUE_DEBUG 0

namespace maligner_dp {

  using Constants::INF;

  int MapData::num_copies = 0;
  int MapData::num_constructs = 0;

  class ColRowPair {
  public:
    ColRowPair(int col, int row) :
      col_(col),
      row_(row) {};
    int col_;
    int row_;

    bool operator<(const ColRowPair& rhs) const {
      return (col_ < rhs.col_) || (!(rhs.col_ < col_) && (row_ < rhs.row_));
    }

  };

  double sizing_penalty(int query_size, int ref_size, const AlignOpts& align_opts) {

    /* TODO: This can be baked into the dynamic programming routine */
    /* For each reference fragment, compute the standard deviation ahead of time */
    
    double delta = query_size - ref_size;
    double sd = align_opts.sd_rate * ref_size;

    if (sd < align_opts.min_sd) {
      sd = align_opts.min_sd;
    }

    double sd_1 = 1.0 / sd;
    double penalty = delta * sd_1;
    penalty = penalty * penalty;
    return penalty;

  }

  void build_trail(ScoreCell* pCell, ScoreCellPVec& trail) {

    #if BUILD_TRAIL_DEBUG > 0
    std::cerr << "Building trail from " << pCell << std::endl;
    #endif

    ScoreCell* pCur = pCell;
    while (pCur != nullptr) {
      #if BUILD_TRAIL_DEBUG > 0
      std::cerr << "pushing back: " << pCur << std::endl;
      #endif
      trail.push_back(pCur);
      pCur = pCur->backPointer_;
    }
  }




  void print_chunk_trail(const ChunkVec& query_chunks, const ChunkVec& ref_chunks) {

    assert(query_chunks.size() == ref_chunks.size());
    int cs = query_chunks.size();
    cerr << "\n\ntrail: ";
    for(int i = 0; i < cs; i++) {

      const Chunk& qc = query_chunks[i];
      const Chunk& rc = ref_chunks[i];

      #if DEBUG > 0
      cerr << "q: " << qc << " r: " << rc << "\n";
      #endif

    }

  }

  std::ostream& operator<<(std::ostream& os, const AlignOpts& ao) {
      
      os << "AlignOpts:"
         <<  "\n\t query_miss_penalty: " << ao.query_miss_penalty
         << "\n\t ref_miss_penalty: " << ao.ref_miss_penalty
         << "\n\t query_max_misses: " << ao.query_max_misses
         << "\n\t ref_max_misses: " << ao.ref_max_misses
         << "\n\t query_max_miss_rate: " << ao.query_max_miss_rate
         << "\n\t ref_max_miss_rate: " << ao.ref_max_miss_rate
         << "\n\t sd_rate: " << ao.sd_rate
         << "\n\t min_sd: " << ao.min_sd
         << "\n\t max_chunk_sizing_error: " << ao.max_chunk_sizing_error
         << "\n\t alignments_per_reference: " << ao.alignments_per_reference 
         << "\n\t min_alignment_spacing: " << ao.min_alignment_spacing 
         << "\n\t neighbor_delta: " << ao.neighbor_delta 
         << "\n\t rescale_query: " << ao.rescale_query
         << "\n\t query_is_bounded: " << ao.query_is_bounded 
         << "\n\t ref_is_bounded: " << ao.ref_is_bounded
         << "\n";
      return os;
  }



}
