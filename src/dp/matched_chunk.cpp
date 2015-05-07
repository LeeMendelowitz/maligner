#include <vector>
#include <ostream>

#include "matched_chunk.h"


namespace maligner_dp {

  std::ostream& operator<<(std::ostream& os, const MatchedChunk& chunk) {

    const Chunk& q = chunk.query_chunk;
    const Chunk& r = chunk.ref_chunk;

    os  << q.start << "," << q.end << "," << q.size << ","
        << r.start << "," << r.end << "," << r.size;

    return os;
  }

  std::ostream& operator<<(std::ostream& os, const Score& score) {

    os << "(" << score.query_miss_score
       << ", " << score.ref_miss_score
       << ", " << score.sizing_score
       << ")";

    return os;
  }

  Score sum_scores(const MatchedChunkVec& mcs) {
    Score score;
    const size_t N = mcs.size();
    for(size_t i = 0; i < N; i++) {
      score += mcs[i].score;
    }
    return score;
  }
  
}

