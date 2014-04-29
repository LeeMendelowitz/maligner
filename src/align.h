
#include "ScoreCell.h"
#include "types.h"

#include <iostream>
#include <vector>
#include <utility>

// Forward Declarations
class ScoreMatrix;
class Alignment;

typedef std::vector<IntVec> PartialSums;
typedef std::vector<Alignment *> AlignmentPVec;

class AlignOpts {

public:

  AlignOpts(double p1, double p2, int p3, int p4,
            double p5, double p6, int p7 = 1, int p8 = 1,
            bool p9 = true) : 
    query_miss_penalty(p1), // penalty for having a site in query unaligned to reference
    ref_miss_penalty(p2), // penalty for having a site in reference unaligned to query
    query_max_misses(p3),
    ref_max_misses(p4),
    sd_rate(0.1), // Fraction of reference fragment to use as standard deviation
    min_sd(p5), // minimum standard deviation imposed in sizing error model, bp
    max_chunk_sizing_error(p6),
    alignments_per_reference(p7),
    min_alignment_spacing(p8),
    query_is_bounded(p9)
  {};

  double query_miss_penalty;
  double ref_miss_penalty;
  int query_max_misses;
  int ref_max_misses;
  double sd_rate;
  double min_sd;
  double max_chunk_sizing_error;

  int alignments_per_reference; // max number of alignments per reference
  int min_alignment_spacing; // minimum amount of spacing between multiple accepted
                            // alignments to the same reference.
  bool query_is_bounded; // true if the first/last fragments in query map are bounded by restriction sites.
};

double sizing_penalty(int query_size, int ref_size, const AlignOpts& align_opts);

//////////////// ///////////////////////////////////////////////////////////////////////
// Bundle the alignment options and data structures
// into a single object that
// can be passed around.
// These structures are for a single alignment task (i.e. one query to one reference)
// Note this just collects non-const pointers
// to external objects into a single object
class AlignTask {
  
public:

  AlignTask(IntVec& q, IntVec& r, PartialSums& qps,
            PartialSums& rps, ScoreMatrix * m,
            AlignmentPVec& alns,
            AlignOpts& ao) :
    query(&q),
    ref(&r),
    query_partial_sums(&qps),
    ref_partial_sums(&rps),
    mat(m),
    alignments(&alns),
    align_opts(&ao)
  {
  }

  // Do we need this copy constructor??
  AlignTask(const AlignTask& other) :
    query(other.query),
    ref(other.ref),
    query_partial_sums(other.query_partial_sums),
    ref_partial_sums(other.ref_partial_sums),
    mat(other.mat),
    alignments(other.alignments),
    align_opts(other.align_opts)
  {
  }

  ~AlignTask() {
  }

  IntVec * query;
  IntVec * ref;
  PartialSums * query_partial_sums;
  PartialSums * ref_partial_sums;
  ScoreMatrix * mat;
  AlignmentPVec * alignments;
  AlignOpts * align_opts;
};

std::ostream& print_align_task(std::ostream& os, const AlignTask& task);


class Chunk {
public:
  
  Chunk() : start(0), end(0), size(0), is_boundary(false) {};

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

  bool operator==(const Score& other ) const {

    return (query_miss_score == other.query_miss_score &&
            ref_miss_score == other.ref_miss_score &&
            sizing_score == other.sizing_score );
  }

  double query_miss_score;
  double ref_miss_score;
  double sizing_score;

};

std::ostream& operator<<(std::ostream& os, const Score& score);

class MatchedChunk {
public:

  MatchedChunk() {};

  MatchedChunk(Chunk& qc, Chunk& rc, Score& s) :
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

};


typedef std::vector<Chunk> ChunkVec;
typedef std::vector<MatchedChunk> MatchedChunkVec;
std::ostream& operator<<(std::ostream& os, const Chunk& chunk);
std::ostream& operator<<(std::ostream& os, const MatchedChunk& chunk);

class Alignment {
public:

  //Alignment(MatchedChunkVec& mc) : matched_chunks(mc) {};
  Alignment() {
    //std::cerr <<"Default constructor\n";
    reset_stats();
  }

  Alignment(MatchedChunkVec&& mc, Score& s) :
    matched_chunks(std::move(mc)),
    rescaled_matched_chunks(matched_chunks),
    score(s),
    query_scaling_factor(1.0)
  {
    //std::cerr << "Alignment constructor from matched_chunks\n";
    summarize();
  }

  Alignment(const Alignment& a) :
    matched_chunks(a.matched_chunks),
    rescaled_matched_chunks(matched_chunks),
    score(a.score),
    query_scaling_factor(1.0)
  {
    //std::cerr << "Alignment copy constructor\n";
    summarize();
  }

  Alignment(Alignment&& a) :
    matched_chunks(std::move(a.matched_chunks)),
    rescaled_matched_chunks(matched_chunks),
    score(a.score),
    query_scaling_factor(1.0) {
    //std::cerr << "Alignment move constructor!\n";
    summarize();
  }


  // rescale the query chunks using the query_scaling_factor, and
  // recompute the sizing error for those chunks.
  void rescale_matched_chunks(AlignOpts& align_opts);

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

    query_scaling_factor = ((double) ref_interior_size) / query_interior_size;

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
    total_miss_rate = ((double) ref_misses + query_misses) / ((double) ref_misses + query_misses + 2.0*num_matched_sites);
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
    query_scaling_factor = 0;

  }

  // Attributes

  MatchedChunkVec matched_chunks;
  MatchedChunkVec rescaled_matched_chunks;
  Score score;
  Score rescaled_score;

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
  double query_scaling_factor;
};



std::ostream& operator<<(std::ostream& os, const Alignment& aln);




////////////////////////////////
// Lower Level functions
////////////////////////////////

/*
  Populate a score matrix using dynamic programming for ungapped alignment.

  The score matrix should have the same number of columns as the reference.

  The ScoreMatrix should already have the same nubmer of columns as the reference,
  and should have enough rows to accomodate the query.
*/
void fill_score_matrix(const AlignTask& task);

// Build the trail which starts at pCell by following its backpointers.
void build_trail(ScoreCell* pCell, ScoreCellPVec& trail);

// Create a vector of query chunks and reference chunks for the given trail.
void build_chunk_trail(const AlignTask& task, ScoreCellPVec& trail, ChunkVec& query_chunks, ChunkVec& ref_chunks);

// Print the vector of query chunks and reference chunks.
void print_chunk_trail(const ChunkVec& query_chunks, const ChunkVec& ref_chunks);

// Build the trail for the best alignment.
bool get_best_alignment(const AlignTask& task, ScoreCellPVec& trail);

// Get the best alignments and store them in the task.
int get_best_alignments(const AlignTask& task);

// Make and return an alignment from the trail through the
// score matrix.
Alignment * alignment_from_trail(const AlignTask& task, ScoreCellPVec& trail);

// Build an alignment by tracing back from ScoreCell.
Alignment * alignment_from_cell(const AlignTask& task, ScoreCell* p_cell);

PartialSums make_partial_sums(const IntVec& frags, const int missed_sites);
PartialSums* make_partial_sums_new(const IntVec& frags, const int missed_sites);

void fill_score_matrix_using_partials(const AlignTask& align_task);


/////////////////////////////////////////
// Friendly functions
/////////////////////////////////////////


// Fill score matrix, find best alignment, and return it.
Alignment * make_best_alignment(const AlignTask& task);
Alignment * make_best_alignment_using_partials(const AlignTask& task);
int make_best_alignments_using_partials(const AlignTask& task);


