#include <iostream>
#include <vector>
#include <utility>
#include <memory>
#include <string>
#include <cassert>
#include <queue>
#include <ostream>

#include "ScoreMatrix.h"
#include "ScoreCell.h"
#include "types.h"
#include "alignment.h"
#include "map_data.h"
#include "map.h"

namespace maligner_dp {

  using std::shared_ptr;
  using namespace maligner_maps;

  typedef  MapData* MapDataPtr;


  class PartialSums {
    // This precomputes chunk lengths to assist in speeding up dynamic programming.
    // It behaves as a 2D array with chunk lengths, indexed the following way:
    // ps(i, j): The length of the chunk of j fragments ending with fragment i.
    //
    // Example:
    //     -----------|-----------|-----------|-----|-----
    //      f[i-2]      f[i-1]     f[i]
    //
    // ps(i, 0) = f[i]
    // ps(i, 1) = f[i] + f[i-1]
    // ps(i, 2) = f[i] + f[i-1] + f[i-2]
    //
    // This indexing may be counterintuitive, but it works more seamlessly with the dynamic programming

    friend class SDInv2;
    IntVec d_;
    const size_t m_;

  public:

    struct forward_tag {};
    struct reverse_tag {};

    PartialSums(const IntVec& d, int max_miss) : m_(max_miss + 1) {
      fill_forward(d, max_miss); 
    }

    PartialSums(const IntVec& d, int max_miss, forward_tag) : m_(max_miss + 1) {
      fill_forward(d, max_miss);
    }

    PartialSums(const IntVec& d, int max_miss, reverse_tag) : m_(max_miss + 1) {
      fill_reverse(d, max_miss);
    }

    int operator()(size_t i, size_t num_miss) const {
      return d_[i*m_ + num_miss];
    }

  private:

    void fill_forward(const IntVec& d, int max_miss) {
      const size_t n = d.size();
      d_ = IntVec(n*m_);

      for(size_t s = 0; s < n; s++) {

        size_t offset = s*m_;
        
        // Constraints: s - i >= 0 -> i <= s -> i < s + 1
        // i < m_
        const size_t L = min(m_, s + 1);
        int a = 0;
        for(size_t i = 0; i < L; i++) {
          a += d[s - i];
          d_[offset+i] = a;
        }

      }

    }

    void fill_reverse(const IntVec& d, int max_miss) {
      const size_t n = d.size();
      d_ = IntVec(n*m_);

      for(size_t s = 0; s < n; s++) {

        const size_t sr = n - s - 1; // Extract right to left from d
        size_t offset = s*m_; // but fill d_ left to right
      
        // Constraints: sr - i >= 0 -> i <= sr -> i < sr + 1
        // i < m_
        const size_t L = min(m_, sr + 1);
        int a = 0;
        for(size_t i = 0; i < L; i++) {
          a += d[sr - i];
          d_[offset+i] = a;
        }

      }

    }

  };

  class SDInv2 {

    // Precompute 1/(sd^2) for each reference chunk

    public:
      SDInv2(const PartialSums& ps, double sd_rate, double min_sd) :
        d_(ps.d_.size()),
        m_(ps.m_) {

          const size_t n = d_.size();
          for(size_t i = 0; i < n; i++) {
            int ref_size = ps.d_[i];
            double sd = sd_rate * ref_size;
            if (sd < min_sd) { sd = min_sd; }
            double sd_1 = 1.0/sd;
            d_[i] = sd_1*sd_1;
          }

     }     

    DoubleVec d_;
    const size_t m_;
  };
  

  ////////////////////////////////////////////////////////////
  // These wrappers around a maligner_maps::Map
  // compute the MapData (i.e. metadata)
  // and potentially PartialSums and SDInv2, 
  // depending on which wrapper is used.
  // This is a convenience to bundle all map data together

  struct MapWrapper {

    MapWrapper(const Map& m, int num_missed_sites) :
      m_(m),
      md_(m_.name_, m_.frags_.size())
    {

    }

    Map m_;
    MapData md_;

  };

  struct QueryMapWrapper {

    QueryMapWrapper(const Map& m, int num_missed_sites) :
        m_(m),
        md_(m_.name_, m_.frags_.size()),
        ps_forward_(m_.frags_, num_missed_sites, PartialSums::forward_tag()),
        ps_reverse_(m_.frags_, num_missed_sites, PartialSums::reverse_tag())
      {

      }

      Map m_;
      MapData md_;
      PartialSums ps_forward_;
      PartialSums ps_reverse_;

  };

  struct RefMapWrapper {

    RefMapWrapper(const Map& m, int num_missed_sites, double sd_rate, double min_sd) :
      m_(m),
      md_(m_.name_, m_.frags_.size()),
      ps_(m_.frags_, num_missed_sites),
      sd_inv_2_(ps_, sd_rate, min_sd )
    {

    }

    Map m_;
    MapData md_;
    PartialSums ps_;
    SDInv2 sd_inv_2_;

  };
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
              int max_alignment_seeds_in = 1000,
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
      max_alignment_seeds(max_alignment_seeds_in),
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
    int max_alignment_seeds;

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
              const SDInv2* rsd,
              ScoreMatrixPtr m,
              AlignmentVec * alns,
              bool is_forward_in,
              AlignOpts& ao) :
      query_map_data(qmd),
      ref_map_data(rmd),
      query(q),
      ref(r),
      query_partial_sums(qps),
      ref_partial_sums(rps),
      ref_sd_inv_2(rsd),
      ref_offset(0),
      mat(m),
      alignments(alns),
      is_forward(is_forward_in),
      align_opts(&ao) { }

    AlignTask(const MapDataPtr qmd,
              const MapDataPtr rmd,
              const std::vector<int>* q,
              const std::vector<int>* r,
              const PartialSums* qps,
              const PartialSums* rps,
              const SDInv2* rsd,
              int ref_offset_in,
              ScoreMatrixPtr m,
              AlignmentVec * alns,
              bool is_forward_in,
              AlignOpts& ao) :
      query_map_data(qmd),
      ref_map_data(rmd),
      query(q),
      ref(r),
      query_partial_sums(qps),
      ref_partial_sums(rps),
      ref_sd_inv_2(rsd),
      ref_offset(ref_offset_in),
      mat(m),
      alignments(alns),
      is_forward(is_forward_in),
      align_opts(&ao)
    { }

    const MapDataPtr query_map_data;
    const MapDataPtr ref_map_data;
    const std::vector<int>* query; // query fragments
    const std::vector<int>* ref; // reference fragments
    const PartialSums* query_partial_sums;
    const PartialSums* ref_partial_sums;
    const SDInv2* ref_sd_inv_2; // 1.0/(sd^2) precomputed for each reference chunk.
    int ref_offset; // index of the first fragment in ref. This will be nonzero if aligning to slice of reference.
    ScoreMatrixPtr mat;
    AlignmentVec* alignments; // Alignment vector to append all found alignments to.
    bool is_forward; // true if the alignment is forward in the reference, false otherwise.
    SizingPenaltyType sizing_penalty;
    AlignOpts * align_opts;

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


  ///////////////////////////////////////////////////////////////////////////
  // A BitVector for marking locations that have been covered in a reference.
  class BitCover {

  public:
    
    BitCover(size_t n) : cover_(n, false) {};

    // Mark locations from begin to end as covered.
    void cover(size_t begin, size_t end) {
      const size_t e = std::max(std::min(end, cover_.size()), begin);
      for(size_t i = begin; i < e; i++) {
        cover_[i] = true; 
      }
    }

    // Return true if any part of the interval is covered
    bool is_covered(size_t begin, size_t end) {
      for(size_t i = begin; i < end; i++) {
        if (cover_[i]) return true;
      }
      return false;
    }

    bool is_covered(size_t ind) {
      return cover_[ind];
    }

    bool operator[](size_t ind) {
      return cover_[ind];
    }

    // A safer version of cover that checks bounds
    void cover_safe(size_t start, size_t end) {

      size_t e = (start < end) ? end : start;

      if (e > cover_.size()) {
        cover_.resize(e, false);
      }

      for(size_t i = start; i < e; i++) {
        cover_[i] = true;
      }

    }

    void reset() {
      const size_t N = cover_.size();
      reset(N);
    }

    void reset(size_t n) {
      cover_.clear();
      cover_.resize(n, false);
    }


  private:
    std::vector<bool> cover_;

  };

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


