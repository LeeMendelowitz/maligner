#include "types.h"
#include "ScoreCell.h"

// Forward Declarations
class ScoreMatrix;


class AlignOpts {

public:

  AlignOpts(double qmp, double rmp, int qmm, int rmm) : 
    query_miss_penalty(qmp),
    ref_miss_penalty(rmp),
    query_max_misses(qmm),
    ref_max_misses(rmm)
  {};

  double query_miss_penalty;
  double ref_miss_penalty;
  int query_max_misses;
  int ref_max_misses;

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

typedef std::vector<Chunk> ChunkVec;
std::ostream& operator<<(std::ostream& os, const Chunk& chunk);

double sizing_penalty(int query_size, int ref_size, const AlignOpts& align_opts);


/*
  Populate a score matrix using dynamic programming for ungapped alignment.

  The score matrix should have the same number of columns as the reference.

  The ScoreMatrix should already have the same nubmer of columns as the reference,
  and should have enough rows to accomodate the query.
*/
void fill_score_matrix(AlignTask& task);
void build_trail(ScoreCell* pCell, ScoreCellPVec& trail);
void build_chunk_trail(AlignTask& task, ScoreCellPVec& trail, ChunkVec& query_chunks, ChunkVec& ref_chunks);
void print_chunk_trail(const ChunkVec& query_chunks, const ChunkVec& ref_chunks);
bool get_best_alignment(AlignTask& task, ScoreCellPVec& trail);


