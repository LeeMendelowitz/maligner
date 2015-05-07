#ifndef MATCHEDCHUNK_H
#define MATCHEDCHUNK_H
#include <vector>

#include "chunk.h"


namespace maligner_dp {

  // Score for a Matched Chunk
  class Score {
  public:
    
    Score() :
      query_miss_score(0.0),
      ref_miss_score(0.0),
      sizing_score(0.0) {};

    Score(double qms, double rms, double ss) :
      query_miss_score(qms),
      ref_miss_score(rms),
      sizing_score(ss) {};

    double total() const {
      return query_miss_score + ref_miss_score + sizing_score;
    }

    Score operator+(const Score& other) const {
      return Score(query_miss_score + other.query_miss_score,
                   ref_miss_score + other.ref_miss_score,
                   sizing_score + other.sizing_score);
    }

    Score& operator+=(const Score& other) {
      query_miss_score += other.query_miss_score;
      ref_miss_score += other.ref_miss_score;
      sizing_score += other.sizing_score;
      return *this;
    }

    bool operator==(const Score& other ) const {

      return (query_miss_score == other.query_miss_score &&
              ref_miss_score == other.ref_miss_score &&
              sizing_score == other.sizing_score );
    }

    double query_miss_score;
    double ref_miss_score;
    double sizing_score;

  };


  class MatchedChunk {
  public:

    MatchedChunk() {};

    MatchedChunk(const Chunk& qc, const Chunk& rc, const Score& s) :
      query_chunk(qc), ref_chunk(rc), score(s) {};

    Chunk query_chunk;
    Chunk ref_chunk;
    Score score;

    bool is_boundary() const {
      return (query_chunk.is_boundary || ref_chunk.is_boundary);
    }

    bool operator==(const MatchedChunk& other) const {
      return (query_chunk == other.query_chunk &&
              ref_chunk == other.ref_chunk &&
              score == other.score);
    }

    int ref_start() const {
      return ref_chunk.start;
    }

    int ref_end() const {
      return ref_chunk.end;
    }

    int ref_size() const {
      return ref_chunk.size;
    }

    int query_start() const {
      return query_chunk.start;
    }

    int query_end() const {
      return query_chunk.end;
    }

    int query_size() const {
      return query_chunk.size;
    }


  };

  typedef std::vector<MatchedChunk> MatchedChunkVec;

  Score sum_scores(const MatchedChunkVec&);

  std::ostream& operator<<(std::ostream& os, const MatchedChunk& chunk);
  std::ostream& operator<<(std::ostream& os, const Score& score);


  

  
}

#endif

