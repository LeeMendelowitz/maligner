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
#define FILL_DEBUG 0
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

  
  // Make partial sums of the preceeding fragment sizes, up to (missed_sites + 1) fragments.
  PartialSums make_partial_sums(const IntVec& frags, const int missed_sites) {
  /*

    Consider fragments with indices i and fragment sizes f
    ...|---------|----------|--------|...
         i-2        i-1        i
         f_(i-2)    f_(i-1)   f_(i)

    The partial sums fragment i for the case missed_sites = 2 will be:
     [ f_i, f_i + f_(i-1), f_i + f_(i-1) + f(i-2)] 

   */
    const IntVec zero_sums(missed_sites + 1, 0);
    const int num_frags = frags.size();

    PartialSums partial_sums(num_frags, zero_sums);

    for (int i = 0; i < num_frags; i++) {
      IntVec& ps = partial_sums[i];
      const int lower_index = i - missed_sites; // inclusive
      int ind = 0;
      int cur_sum = 0;
      for (int j = i, ind = 0; j >= lower_index; j--, ind++) {
        if (j < 0) break;
        cur_sum += frags[j];
        ps[ind] = cur_sum;
      }
    }

    return partial_sums;
  }

  // Make partial sums of the preceeding fragment sizes, up to (missed_sites + 1) fragments.
  PartialSumsPtr make_partial_sums_new(const IntVec& frags, const int missed_sites) {
  /*

    Consider fragments with indices i and fragment sizes f
    ...|---------|----------|--------|...
         i-2        i-1        i
         f_(i-2)    f_(i-1)   f_(i)

    The partial sums fragment i for the case missed_sites = 2 will be:
     [ f_i, f_i + f_(i-1), f_i + f_(i-1) + f(i-2)] 

   */
    const IntVec zero_sums(missed_sites + 1, 0);
    const int num_frags = frags.size();

    PartialSumsPtr p_partial_sums = PartialSumsPtr( new PartialSums(num_frags, zero_sums) );

    for (int i = 0; i < num_frags; i++) {
      IntVec& ps = (*p_partial_sums)[i];
      const int lower_index = i - missed_sites; // inclusive
      int ind = 0;
      int cur_sum = 0;
      for (int j = i, ind = 0; j >= lower_index; j--, ind++) {
        if (j < 0) break;
        cur_sum += frags[j];
        ps[ind] = cur_sum;
      }
    }

    return p_partial_sums;
  }

}
