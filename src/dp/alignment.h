#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#include <vector>
#include <ostream>

#include "map_data.h"
#include "matched_chunk.h"

namespace maligner_dp {

  using std::size_t;

  // Forward Declarations
  class AlignOpts;

  class AlignmentHeader {};
  std::ostream& operator<<(std::ostream& os, const AlignmentHeader& a);

  class Alignment {
  public:

    //Alignment(MatchedChunkVec& mc) : matched_chunks(mc) {};

    // We must provide a default constructor because
    // we have provided a specialized constructor.
    Alignment() : is_valid(false) {
    }

    Alignment(const MatchedChunkVec& mc, const Score& s,
      const MapData& query_md,
      const MapData& ref_md,
      bool is_forward__) :
      query_map_data(query_md),
      ref_map_data(ref_md),
      matched_chunks(mc),
      rescaled_matched_chunks(matched_chunks),
      score(s),
      rescaled_score(s),
      query_scaling_factor(1.0),
      is_forward(is_forward__)
    {
      summarize();
    }

    
    // Use the default copy constructor/assignment
    Alignment(const Alignment&) = default;
    Alignment& operator=(const Alignment&) = default;

    // Use the default move constructor/assignment
    Alignment(Alignment&&) = default;
    Alignment& operator=(Alignment&&) = default;

    // rescale the query chunks using the query_scaling_factor, and
    // recompute the sizing error for those chunks.
    void rescale_matched_chunks(const AlignOpts& align_opts);

    // Compute summary statistics from matched chunks.
    void summarize(); 

    // Reset the summary
    void reset_stats();

    int ref_start() const {
      if (is_valid) {
        const MatchedChunk& first = matched_chunks.front();
        return first.ref_start();
      }
      return 0;
    }

    int ref_end() const {
      if (is_valid) {
        const MatchedChunk& last = matched_chunks.back();
        return last.ref_end();
      }
      return 0;
    }

    void flip_query_coords();

    // Attributes
    MapData query_map_data;
    MapData ref_map_data;
    MatchedChunkVec matched_chunks;
    MatchedChunkVec rescaled_matched_chunks;
    Score score;
    Score rescaled_score;

    double total_score;
    double total_rescaled_score;

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
    bool is_forward;
    bool is_valid;

  };

  typedef std::vector<Alignment> AlignmentVec;

    // Sort an AlignmentVec in ascending order of score.
  class AlignmentRescaledScoreComp {
  public:
      bool operator()(const Alignment& a1, const Alignment& a2) {
        return a1.total_rescaled_score < a2.total_rescaled_score;
      }
  };

  // Sort an AlignmentVec in ascending order of score.
  class AlignmentScoreComp {
  public:
      bool operator()(const Alignment& a1, const Alignment& a2) {
        return a1.total_score < a2.total_score;
      }
  };

  extern const Alignment INVALID_ALIGNMENT;


  inline void Alignment::summarize() {

      query_misses = 0;
      ref_misses = 0;
      query_interior_size = 0;
      ref_interior_size = 0;
      num_matched_sites = 0;

      total_score = score.total();
      total_rescaled_score = rescaled_score.total();

      is_valid = !matched_chunks.empty();

      // Fix query coordinates to match forward coordinates of query
      if(!is_forward) {
        flip_query_coords(); 
      }

      std::size_t l(matched_chunks.size());
      for (std::size_t i = 0; i < l; i++) {
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

  inline void Alignment::reset_stats() {
      is_valid = false;
      is_forward = false;
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
      total_score = 0;
      total_rescaled_score = 0;
  }

  // Flip the query coordinates in every matched chunk
  inline void Alignment::flip_query_coords() {
    const size_t num_chunks {matched_chunks.size()};
    const size_t num_query_frags = query_map_data.num_frags_;
    for (size_t i = 0; i < num_chunks; i++) {

      MatchedChunk& mc1 = matched_chunks[i];
      MatchedChunk& mc2 = rescaled_matched_chunks[i];

      Chunk& qc = mc1.query_chunk;
      std::swap(qc.start, qc.end);
      qc.start = num_query_frags - qc.start;
      qc.end = num_query_frags - qc.end;

      qc = mc2.query_chunk;
      std::swap(qc.start, qc.end);
      qc.start = num_query_frags - qc.start;
      qc.end = num_query_frags - qc.end;

    }
  }


  void print_alignment(std::ostream& os, const Alignment& aln);

  std::ostream& operator<<(std::ostream& os, const Alignment& aln);

}

#endif



 