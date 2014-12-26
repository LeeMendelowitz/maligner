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
  
}

