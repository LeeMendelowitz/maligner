#include <iostream>
#include <vector>
#include <utility>
#include <memory>
#include <string>
#include <cassert>
#include <queue>
#include <ostream>
#include <limits>

#include "ScoreMatrix.h"
#include "ScoreCell.h"
#include "types.h"
#include "alignment.h"
#include "map_data.h"
#include "map.h"
#include "map_wrappers.h"
#include "partialsums.h"
#include "bitcover.h"

namespace maligner_dp {

  using std::shared_ptr;
  using namespace maligner_maps;

  typedef  MapData* MapDataPtr;



  //////////////////////////////////////////////////////////////
    

  class AlignOpts {

  public:

    AlignOpts(double query_miss_penality_in,
              double ref_miss_penalty_in,
              int query_max_misses_in,
              int ref_max_misses_in,
              double sd_rate_in,
              double min_sd_in,
              double max_chunk_sizing_error_in,
              double ref_max_miss_rate_in,
              double query_max_miss_rate_in,
              int alignments_per_reference_in,
              int min_alignment_spacing_in,
              int neighborhood_delta_in,
              bool query_is_bounded_in,
              bool ref_is_bounded_in,
              bool rescale_query_in) : 
      query_miss_penalty(query_miss_penality_in), // penalty for having a site in query unaligned to reference
      ref_miss_penalty(ref_miss_penalty_in), // penalty for having a site in reference unaligned to query
      query_max_misses(query_max_misses_in),
      ref_max_misses(ref_max_misses_in),
      sd_rate(sd_rate_in), // Fraction of reference fragment to use as standard deviation
      min_sd(min_sd_in), // minimum standard deviation imposed in sizing error model, bp
      max_chunk_sizing_error(max_chunk_sizing_error_in),
      ref_max_miss_rate(ref_max_miss_rate_in),
      query_max_miss_rate(query_max_miss_rate_in),
      alignments_per_reference(alignments_per_reference_in),
      min_alignment_spacing(min_alignment_spacing_in),
      neighbor_delta(neighborhood_delta_in), // Return alignments within +/- neighbor_delta ref. fragments
                              // of each selected alignment.
      rescale_query(rescale_query_in), // Rescale the query chunks post-alignment & adjust sizing score.
      query_is_bounded(query_is_bounded_in),
      ref_is_bounded(ref_is_bounded_in),
      query_miss_penalties(query_max_misses + 1, 0.0),
      ref_miss_penalties(ref_max_misses + 1, 0.0)
    {

      // Initialize query_miss_penalties array
      for(int i = 0;i < query_max_misses + 1; ++i)
      {
        query_miss_penalties[i] = i*query_miss_penalty;
      }

      // Initialize ref_miss_penalties array
      for(int i = 0; i < ref_max_misses + 1; ++i)
      {
        ref_miss_penalties[i] = i*ref_miss_penalty;
      }

    };

    double query_miss_penalty;
    double ref_miss_penalty;
    int query_max_misses;
    int ref_max_misses;
    double sd_rate;
    double min_sd;
    double max_chunk_sizing_error;
    double ref_max_miss_rate; // Maximum unmatched site rate in the reference
    double query_max_miss_rate; // Maximum unmatched site rate in the query

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

  std::ostream& operator<<(std::ostream& os, const AlignOpts& ao);


  //////////////////////////////////////////////////////
  // Sizing Penalty Definitions

  double sizing_penalty(int query_size, int ref_size, const AlignOpts& align_opts);

  class Chi2SizingPenalty {
  public:

    double operator()(int query_size, int ref_size, const AlignOpts& align_opts) const {

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

  };

  class NoSizingPenalty {
  public:
    double operator()(int, int, const AlignOpts&) const {
      return 0.0;
    }
  };
  /////////////////////////////////////////////////////


  //////////////// ///////////////////////////////////////////////////////////////////////
  // Bundle the alignment options and data structures
  // into a single object that
  // can be passed around.
  // These structures are for a single alignment task (i.e. one query to one reference)
  // Note this just bundles non-const pointers
  // to external objects into a single object

  template<class ScoreMatrixType, class SizingPenaltyType>
  class AlignTask {

    typedef ScoreMatrixType* ScoreMatrixPtr;

  public:

    AlignTask(const MapDataPtr qmd,
              const MapDataPtr rmd,
              const std::vector<int>* q,
              const std::vector<int>* r,
              const PartialSums* qps,
              const PartialSums* rps,
              const SDInv* rsd,
              const IntVec* query_ix_to_locs_in,
              const IntVec* ref_ix_to_locs_in,
              ScoreMatrixPtr m,
              AlignmentVec * alns,
              bool query_is_forward_in,
              bool ref_is_forward_in,
              bool is_circular_in,
              AlignOpts& ao) :
      query_map_data(qmd),
      ref_map_data(rmd),
      query(q),
      ref(r),
      query_partial_sums(qps),
      ref_partial_sums(rps),
      ref_sd_inv(rsd),
      query_ix_to_locs(query_ix_to_locs_in),
      ref_ix_to_locs(ref_ix_to_locs_in),      
      ref_offset(0),
      mat(m),
      alignments(alns),
      query_is_forward(query_is_forward_in),
      ref_is_forward(ref_is_forward_in),
      align_opts(&ao) { compute_max_misses(); }

    // Takes an additional parameter: ref_offset_in
    AlignTask(const MapDataPtr qmd,
              const MapDataPtr rmd,
              const std::vector<int>* q,
              const std::vector<int>* r,
              const PartialSums* qps,
              const PartialSums* rps,
              const SDInv* rsd,
              const IntVec* query_ix_to_locs_in,
              const IntVec* ref_ix_to_locs_in,
              int ref_offset_in,
              ScoreMatrixPtr m,
              AlignmentVec * alns,
              bool query_is_forward_in,
              bool ref_is_forward_in,
              AlignOpts& ao) :
      query_map_data(qmd),
      ref_map_data(rmd),
      query(q),
      ref(r),
      query_partial_sums(qps),
      ref_partial_sums(rps),
      ref_sd_inv(rsd),
      query_ix_to_locs(query_ix_to_locs_in),
      ref_ix_to_locs(ref_ix_to_locs_in),
      ref_offset(ref_offset_in),
      mat(m),
      alignments(alns),
      query_is_forward(query_is_forward_in),
      ref_is_forward(ref_is_forward_in),
      align_opts(&ao)
    { compute_max_misses(); }

    const MapDataPtr query_map_data;
    const MapDataPtr ref_map_data;
    const std::vector<int>* query; // query fragments
    const std::vector<int>* ref; // reference fragments
    const PartialSums* query_partial_sums;
    const PartialSums* ref_partial_sums;
    const SDInv* ref_sd_inv; // 1.0/(sd^2) precomputed for each reference chunk.
    const IntVec* query_ix_to_locs; // Fragment index to bp location (always oriented forward)
    const IntVec* ref_ix_to_locs; // Fragment index to bp location (always oriented forward)
    int query_max_total_misses;
    int ref_max_total_misses;

    int ref_offset; // index of the first fragment in ref. This will be nonzero if aligning to slice of reference.
    ScoreMatrixPtr mat;
    AlignmentVec* alignments; // Alignment vector to append all found alignments to.

    bool query_is_forward; // true if the query data is given as forward (i.e. query, query_partial_sums)
    bool ref_is_forward; // true if reference data is given as forward(i.e. ref, ref_partial_sums)
    SizingPenaltyType sizing_penalty;
    AlignOpts * align_opts;

    //////////////////////////////////////////
    // Compute the maximum number of allowed misses, given the 
    // settings for query_max_miss_rate and ref_max_miss_rate.
    void compute_max_misses() {

      int num_query_sites = query_map_data->is_bounded_ ? query->size() + 1: query->size()-1;
      query_max_total_misses = int(align_opts->query_max_miss_rate * num_query_sites);

      // Compute an upper bound on the number of allowed reference misses. It depends
      // on the number of sites in the query.
      const double ref_max_miss_rate = align_opts->ref_max_miss_rate;

      ref_max_total_misses = (ref_max_miss_rate < 1.0 ) ?
        (ref_max_miss_rate / (1 - ref_max_miss_rate) * num_query_sites) :
        std::numeric_limits<int>::max();

    };

    bool ref_is_circular() const {
      return ref_map_data->is_circular_;
    };

  };


  template<class ScoreMatrixType, class SizingPenaltyType>
  std::ostream& print_align_task(std::ostream& os, const AlignTask<ScoreMatrixType, SizingPenaltyType>& task);


  ////////////////////////////////
  // Lower Level functions
  ////////////////////////////////
  /*
    Populate a score matrix using dynamic programming for ungapped alignment.

    The score matrix should have the same number of columns as the reference.

    The ScoreMatrix should already have the same number of columns as the reference,
    and should have enough rows to accomodate the query.
  */
  template<class ScoreMatrixType, class SizingPenaltyType>
  void fill_score_matrix(const AlignTask<ScoreMatrixType, SizingPenaltyType>& task);

  template<class ScoreMatrixType, class SizingPenaltyType>
  void fill_score_matrix_using_partials(const AlignTask<ScoreMatrixType, SizingPenaltyType>& align_task);

  template<class ScoreMatrixType, class SizingPenaltyType>
  void fill_score_matrix_using_partials(const AlignTask<ScoreMatrixType, SizingPenaltyType>& align_task, row_order_tag);

  template<class ScoreMatrixType, class SizingPenaltyType>
  void fill_score_matrix_using_partials(const AlignTask<ScoreMatrixType, SizingPenaltyType>& align_task, column_order_tag);

  template<class ScoreMatrixType, class SizingPenaltyType>
  void fill_score_matrix_using_partials_with_cell_queue(const AlignTask<ScoreMatrixType, SizingPenaltyType>& align_task);

  template<class ScoreMatrixType, class SizingPenaltyType>
  void fill_score_matrix_using_partials_with_cell_mark(const AlignTask<ScoreMatrixType, SizingPenaltyType>& align_task);

  template<class ScoreMatrixType, class SizingPenaltyType>
  void fill_score_matrix_using_partials_no_size_penalty(const AlignTask<ScoreMatrixType, SizingPenaltyType>& align_task);

  template<class ScoreMatrixType, class SizingPenaltyType>
  void fill_score_matrix_using_partials_with_size_penalty_class(const AlignTask<ScoreMatrixType, SizingPenaltyType>& align_task);





  // Build the trail which starts at pCell by following its backpointers.
  void build_trail(ScoreCell* pCell, ScoreCellPVec& trail);

  // Create a vector of query chunks and reference chunks for the given trail.
  template<class ScoreMatrixType, class SizingPenaltyType>
  void build_chunk_trail(const AlignTask<ScoreMatrixType, SizingPenaltyType>& task, ScoreCellPVec& trail, ChunkVec& query_chunks, ChunkVec& ref_chunks);

  // Print the vector of query chunks and reference chunks.
  void print_chunk_trail(const ChunkVec& query_chunks, const ChunkVec& ref_chunks);

  // Build the trail for the best alignment.
  template<class ScoreMatrixType, class SizingPenaltyType>
  bool get_best_alignment_trail(const AlignTask<ScoreMatrixType, SizingPenaltyType>& task, ScoreCellPVec& trail);

  // Get the best alignments and store them in the task.
  template<class ScoreMatrixType, class SizingPenaltyType>
  int get_best_alignments(const AlignTask<ScoreMatrixType, SizingPenaltyType>& task);

  template<class ScoreMatrixType, class SizingPenaltyType>
  int get_best_alignments_try_all(const AlignTask<ScoreMatrixType, SizingPenaltyType>& task);

  // Make and return an alignment from the trail through the
  // score matrix.
  template<class ScoreMatrixType, class SizingPenaltyType>
  Alignment alignment_from_trail(const AlignTask<ScoreMatrixType, SizingPenaltyType>& task, ScoreCellPVec& trail);

  // Build an alignment by tracing back from ScoreCell.
  template<class ScoreMatrixType, class SizingPenaltyType>
  Alignment alignment_from_cell(const AlignTask<ScoreMatrixType, SizingPenaltyType>& task, ScoreCell* p_cell);

  PartialSums make_partial_sums(const IntVec& frags, const int missed_sites);
  PartialSums* make_partial_sums_new(const IntVec& frags, const int missed_sites);


  /////////////////////////////////////////
  // Friendly functions
  /////////////////////////////////////////

  // Fill score matrix, find best alignment, and return it.
  template<class ScoreMatrixType, class SizingPenaltyType>
  Alignment make_best_alignment(const AlignTask<ScoreMatrixType, SizingPenaltyType>& task);

  template<class ScoreMatrixType, class SizingPenaltyType>
  Alignment make_best_alignment_using_partials(const AlignTask<ScoreMatrixType, SizingPenaltyType>& task);

  template<class ScoreMatrixType, class SizingPenaltyType>
  int make_best_alignments_using_partials(const AlignTask<ScoreMatrixType, SizingPenaltyType>& task);

  // Assigns mad scores to the alignments, and returns the mad.
  double compute_mad_scores(AlignmentVec& alignments, const int max_alignments_mad, double min_mad);


  class SizingPenalty {
  public:

    SizingPenalty(int ref_size, const AlignOpts& align_opts) {

      ref_size_ = ref_size;
      sd_ = align_opts.sd_rate * ref_size_;
      if( sd_ < align_opts.min_sd ) {
        sd_ = align_opts.min_sd;
      }
      sd_1_ = 1.0 / sd_;

    }

    double operator()(int query_size) {
      double delta = query_size - ref_size_;
      double penalty = delta*sd_1_;
      return penalty*penalty;
    }

  // private:
    int ref_size_;
    double sd_;
    double sd_1_;

  };

}

// Bring in the definitions of templated functions
#include "templated_align_functions.h"


