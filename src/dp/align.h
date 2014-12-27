
#include "ScoreMatrix.h"
#include "ScoreCell.h"
#include "types.h"
#include "alignment.h"
#include "map_data.h"

#include <iostream>
#include <vector>
#include <utility>
#include <memory>
#include <string>

namespace maligner_dp {


  // Forward Declarations
  class ScoreMatrix;

  typedef std::vector<IntVec> PartialSums;
  using std::shared_ptr;

  // Using shared pointers:
  // typedef std::shared_ptr< PartialSums > PartialSumsPtr;
  // typedef std::shared_ptr< MapData > MapDataPtr;
  // typedef std::shared_ptr< Alignment > AlignmentPtr;
  // typedef std::vector< AlignmentPtr > AlignmentPVec;
  // typedef std::shared_ptr< AlignmentPVec > AlignmentPVecPtr;

  // Using raw pointers:
  typedef  PartialSums* PartialSumsPtr;
  typedef  MapData* MapDataPtr;


  class AlignOpts {

  public:

    AlignOpts(double query_miss_penality_in,
              double ref_miss_penalty_in,
              int query_max_misses_in,
              int ref_max_misses_in,
              double sd_rate_in,
              double min_sd_in,
              double max_chunk_sizing_error_in,
              int alignments_per_reference_in = 1,
              int min_alignment_spacing_in = 1,
              int neighborhood_delta_in = 0,
              bool query_is_bounded_in = true,
              bool ref_is_bounded_in = false) : 
      query_miss_penalty(query_miss_penality_in), // penalty for having a site in query unaligned to reference
      ref_miss_penalty(ref_miss_penalty_in), // penalty for having a site in reference unaligned to query
      query_max_misses(query_max_misses_in),
      ref_max_misses(ref_max_misses_in),
      sd_rate(sd_rate_in), // Fraction of reference fragment to use as standard deviation
      min_sd(min_sd_in), // minimum standard deviation imposed in sizing error model, bp
      max_chunk_sizing_error(max_chunk_sizing_error_in),
      alignments_per_reference(alignments_per_reference_in),
      min_alignment_spacing(min_alignment_spacing_in),
      neighbor_delta(neighborhood_delta_in), // Return alignments within +/- neighbor_delta ref. fragments
                              // of each selected alignment.
      rescale_query(true), // Rescale the query chunks post-alignment & adjust sizing score.
      query_is_bounded(query_is_bounded_in),
      ref_is_bounded(ref_is_bounded_in),
      query_miss_penalties(query_max_misses + 1, 0.0),
      ref_miss_penalties(ref_max_misses + 1, 0.0)
    {

      // Initialize query_miss_penalties array
      for(int i = 0, val = 0;
          i < query_max_misses + 1;
          val += query_miss_penalty, ++i)
      {
        query_miss_penalties[i] = val;
      }

      // Initialize ref_miss_penalties array
      for(int i = 0, val = 0;
          i < ref_max_misses + 1;
          val += ref_miss_penalty, ++i)
      {
        ref_miss_penalties[i] = val;
      }

    };

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
    int neighbor_delta;
    bool rescale_query;
    bool query_is_bounded; // true if the first/last fragments in query map are bounded by restriction sites.
    bool ref_is_bounded; // true if the first/last fragments in reference map are bounded by restriction sites.
    
    // Precomputed multiples of query_miss_penalties and ref_miss_penalties
    DoubleVec query_miss_penalties;
    DoubleVec ref_miss_penalties;

  };

  double sizing_penalty(int query_size, int ref_size, const AlignOpts& align_opts);


  //////////////// ///////////////////////////////////////////////////////////////////////
  // Bundle the alignment options and data structures
  // into a single object that
  // can be passed around.
  // These structures are for a single alignment task (i.e. one query to one reference)
  // Note this just bundles non-const pointers
  // to external objects into a single object

  class AlignTask {
    
  public:

    AlignTask(MapDataPtr qmd, MapDataPtr rmd, const std::vector<int>* q,
              const std::vector<int>* r, PartialSumsPtr qps,
              PartialSumsPtr rps,
              ScoreMatrixPtr m,
              AlignmentVec * alns,
              AlignOpts& ao) :
      query_map_data(qmd),
      ref_map_data(rmd),
      query(q),
      ref(r),
      query_partial_sums(qps),
      ref_partial_sums(rps),
      ref_offset(0),
      mat(m),
      alignments(alns),
      align_opts(&ao)
    {
    }

    AlignTask(MapDataPtr qmd, MapDataPtr rmd,
              const std::vector<int>* q,
              const std::vector<int>* r,
              PartialSumsPtr qps,
              PartialSumsPtr rps,
              int ro,
              ScoreMatrixPtr m,
              AlignmentVec * alns,
              AlignOpts& ao) :
      query_map_data(qmd),
      ref_map_data(rmd),
      query(q),
      ref(r),
      query_partial_sums(qps),
      ref_partial_sums(rps),
      ref_offset(ro),
      mat(m),
      alignments(alns),
      align_opts(&ao)
    {
    }

    ~AlignTask() {
    }

    MapDataPtr query_map_data;
    MapDataPtr ref_map_data;
    const std::vector<int>* query; // query fragments
    const std::vector<int>* ref; // reference fragments
    PartialSumsPtr query_partial_sums;
    PartialSumsPtr ref_partial_sums;
    int ref_offset; // index of the first fragment in ref. This will be nonzero if aligning to slice of reference.
    ScoreMatrixPtr mat;
    AlignmentVec* alignments;
    AlignOpts * align_opts;
  };

  std::ostream& print_align_task(std::ostream& os, const AlignTask& task);




  ////////////////////////////////
  // Lower Level functions
  ////////////////////////////////
  /*
    Populate a score matrix using dynamic programming for ungapped alignment.

    The score matrix should have the same number of columns as the reference.

    The ScoreMatrix should already have the same number of columns as the reference,
    and should have enough rows to accomodate the query.
  */
  void fill_score_matrix(const AlignTask& task);
  void fill_score_matrix_using_partials(const AlignTask& align_task);
  void fill_score_matrix_using_partials_with_cell_queue(const AlignTask& align_task);
  void fill_score_matrix_using_partials_with_cell_mark(const AlignTask& align_task);
  void fill_score_matrix_using_partials_no_size_penalty(const AlignTask& align_task);
  void fill_score_matrix_using_partials_with_size_penalty_class(const AlignTask& align_task);

  // Build the trail which starts at pCell by following its backpointers.
  void build_trail(ScoreCell* pCell, ScoreCellPVec& trail);

  // Create a vector of query chunks and reference chunks for the given trail.
  void build_chunk_trail(const AlignTask& task, ScoreCellPVec& trail, ChunkVec& query_chunks, ChunkVec& ref_chunks);

  // Print the vector of query chunks and reference chunks.
  void print_chunk_trail(const ChunkVec& query_chunks, const ChunkVec& ref_chunks);

  // Build the trail for the best alignment.
  bool get_best_alignment_trail(const AlignTask& task, ScoreCellPVec& trail);

  // Get the best alignments and store them in the task.
  int get_best_alignments(const AlignTask& task);

  // Make and return an alignment from the trail through the
  // score matrix.
  Alignment alignment_from_trail(const AlignTask& task, ScoreCellPVec& trail);

  // Build an alignment by tracing back from ScoreCell.
  Alignment alignment_from_cell(const AlignTask& task, ScoreCell* p_cell);

  PartialSums make_partial_sums(const IntVec& frags, const int missed_sites);
  PartialSumsPtr make_partial_sums_new(const IntVec& frags, const int missed_sites);


  /////////////////////////////////////////
  // Friendly functions
  /////////////////////////////////////////

  // Fill score matrix, find best alignment, and return it.
  Alignment make_best_alignment(const AlignTask& task);
  Alignment make_best_alignment_using_partials(const AlignTask& task);
  int make_best_alignments_using_partials(const AlignTask& task);

}


