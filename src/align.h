
#include "ScoreCell.h"
#include "types.h"

#include <iostream>
#include <utility>

// Forward Declarations
class ScoreMatrix;


class AlignOpts {

public:

  AlignOpts(double p1, double p2, int p3, int p4,
            double p5) : 
    query_miss_penalty(p1), // penalty for having a site in query unaligned to reference
    ref_miss_penalty(p2), // penalty for having a site in reference unaligned to query
    query_max_misses(p3),
    ref_max_misses(p4),
    max_chunk_sizing_error(p5)
  {};

  double query_miss_penalty;
  double ref_miss_penalty;
  int query_max_misses;
  int ref_max_misses;
  double max_chunk_sizing_error;

};

// Bundle the alignment options into a single object that
// can be passed around.
class AlignTask {
  
public:
  AlignTask(IntVec& q, IntVec& r, ScoreMatrix& m, AlignOpts& ao) :
    query(q), ref(r), mat(m), align_opts(ao) {};

  IntVec& query;
  IntVec& ref;
  ScoreMatrix& mat;
  AlignOpts& align_opts;
};

class Chunk {
public:
  Chunk(int s, int e, int sz, bool isb) :
    start(s), end(e), size(sz), is_boundary(isb) {};
  int start; // start coordinate, inclusive
  int end; // end coordinate, exclusive
  int size; // total size of chunk (bp)
  bool is_boundary; // is boundary chunk

  bool operator==(const Chunk& other) const {
    return (start == other.start &&
            end == other.end &&
            size == other.size &&
            is_boundary == other.is_boundary);
  }

  int num_misses() const {
    return end - start - 1;
  }

};

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

  double query_miss_score;
  double ref_miss_score;
  double sizing_score;

};

std::ostream& operator<<(std::ostream& os, const Score& score);

class MatchedChunk {
public:

  MatchedChunk(Chunk& qc, Chunk& rc, Score& s) :
    query_chunk(qc), ref_chunk(rc), score(s) {};

  Chunk query_chunk;
  Chunk ref_chunk;
  Score score;

  bool is_boundary() const {
    return (query_chunk.is_boundary || ref_chunk.is_boundary);
  }

};


typedef std::vector<Chunk> ChunkVec;
typedef std::vector<MatchedChunk> MatchedChunkVec;
std::ostream& operator<<(std::ostream& os, const Chunk& chunk);
std::ostream& operator<<(std::ostream& os, const MatchedChunk& chunk);

class Alignment {
public:

  //Alignment(MatchedChunkVec& mc) : matched_chunks(mc) {};
  Alignment() {
    std::cerr <<"Default constructor\n";
    reset_stats();
  }

  Alignment(MatchedChunkVec&& mc, Score& s) : matched_chunks(mc), score(s) {
    std::cerr << "Alignment constructor from matched_chunks\n";
    summarize();
  }

  Alignment(const Alignment& a) : matched_chunks(a.matched_chunks), score(a.score) {
    std::cerr << "Alignment copy constructor\n";
    summarize();
  }

  Alignment(Alignment&& a) : matched_chunks(a.matched_chunks), score(a.score) {
    std::cerr << "Alignment move constructor!\n";
    summarize();
  }

  // Compute summary statistics from matched chunks.
  void summarize() {

    query_misses = 0;
    ref_misses = 0;
    query_interior_size = 0;
    ref_interior_size = 0;
    num_matched_sites = 0;

    const size_t l = matched_chunks.size();
    for (size_t i = 0; i < l; i++) {
      const MatchedChunk& mc = matched_chunks[i];
      query_misses += mc.query_chunk.num_misses();
      ref_misses += mc.ref_chunk.num_misses();
      if (!mc.query_chunk.is_boundary && !mc.ref_chunk.is_boundary) {
        query_interior_size += mc.query_chunk.size;
        ref_interior_size += mc.ref_chunk.size;
      }
    }

    // Count the number of matched sites.
    //  - each non-boundary chunk begins/ends with a matched site.
    //  - The first chunk can potentially begin with a matched site, in local alignment.
    //  - The last chunk can potentially edn with a matched site, in local alignment.

    num_matched_sites = matched_chunks.size()-1;
    if (l > 0 && !matched_chunks[0].is_boundary()) {
      num_matched_sites++; // The first matched_chunk is not a boundary chunk
    }
    if (l > 1 && !matched_chunks[l-1].is_boundary()) {
      num_matched_sites++;
    }

    ref_miss_rate = ((double) ref_misses)/((double) num_matched_sites + ref_misses);
    query_miss_rate = ((double) query_misses)/((double) num_matched_sites + query_misses);
    total_miss_rate = ((double) ref_misses + query_misses) / ((double) ref_misses + query_misses + num_matched_sites);
    interior_size_ratio = ((double) query_interior_size) / ref_interior_size;

  }

  void reset_stats() {
    query_misses = 0;
    ref_misses = 0;
    query_miss_rate = 0;
    ref_miss_rate = 0;
    total_miss_rate = 0;
    query_interior_size = 0;
    ref_interior_size = 0;
    num_matched_sites = 0;
    interior_size_ratio = 0;

  }

  // Attributes

  MatchedChunkVec matched_chunks;
  Score score;

  // summary statistics of an alignment.
  // These are computable from the matched_chunks
  int num_matched_sites;
  int query_misses;
  int ref_misses;
  double query_miss_rate;
  double ref_miss_rate;
  double total_miss_rate;
  int query_interior_size; // total size of non-boundary fragments
  int ref_interior_size; // total size of non-boundary fragments
  double interior_size_ratio;

};

std::ostream& operator<<(std::ostream& os, const Alignment& aln);



double sizing_penalty(int query_size, int ref_size, const AlignOpts& align_opts);

/*
  Populate a score matrix using dynamic programming for ungapped alignment.

  The score matrix should have the same number of columns as the reference.

  The ScoreMatrix should already have the same nubmer of columns as the reference,
  and should have enough rows to accomodate the query.
*/
void fill_score_matrix(AlignTask& task);

// Build the trail which starts at pCell by following its backpointers.
void build_trail(ScoreCell* pCell, ScoreCellPVec& trail);

// Create a vector of query chunks and reference chunks for the given trail.
void build_chunk_trail(AlignTask& task, ScoreCellPVec& trail, ChunkVec& query_chunks, ChunkVec& ref_chunks);

// Print the vector of query chunks and reference chunks.
void print_chunk_trail(const ChunkVec& query_chunks, const ChunkVec& ref_chunks);

// Build the trail for the best alignment.
bool get_best_alignment(AlignTask& task, ScoreCellPVec& trail);

// Make and return an alignment from the trail through the
// score matrix.
Alignment alignment_from_trail(AlignTask& task, ScoreCellPVec& trail);



