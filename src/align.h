
#include "ScoreCell.h"

#include <iostream>
#include <utility>

// Forward Declarations
class ScoreMatrix;


class AlignOpts {

public:

  AlignOpts(double p1, double p2, int p3, int p4,
            double p5) : 
    query_miss_penalty(p1),
    ref_miss_penalty(p2),
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
  Chunk(int s, int e, int sz) :
    start(s), end(e), size(sz) {};
  int start; // start coordinate, inclusive
  int end; // end coordinate, exclusive
  int size; // total size of chunk (bp)

  bool operator==(const Chunk& other) const {
    return (start == other.start &&
            end == other.end &&
            size == other.size);
  }

};

// Score for a Matched Chunk
class Score {
public:
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
  }

  Alignment(MatchedChunkVec&& mc) : matched_chunks(mc) {
    std::cerr << "Alignment constructor from matched_chunks\n";
  }

  Alignment(const Alignment& a) : matched_chunks(a.matched_chunks) {
    std::cerr << "Alignment copy constructor\n";
  }

  Alignment(Alignment&& a) : matched_chunks(a.matched_chunks) {
    std::cerr << "Alignment move constructor!\n";

  }

  Alignment& operator=(const Alignment& a) {
    this->matched_chunks = a.matched_chunks;
    std::cerr << "Alignment copy assignment!\n";
    return *this;
  }

  Alignment& operator=(Alignment&& a) {
    this->matched_chunks = std::move(a.matched_chunks);
    std::cerr << "Alignment move assignment!\n";
    return *this;
  }

  MatchedChunkVec matched_chunks;
  
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



