#include "alignment.h"
#include "align.h"
#include "matched_chunk.h"

#include <ostream>

namespace maligner_dp {

  // rescale the query chunks using the query_scaling_factor, and
  // recompute the sizing error for those chunks.
  void Alignment::rescale_matched_chunks(const AlignOpts& align_opts) {

    //summarize();

    rescaled_score = score;
    rescaled_score.sizing_score = 0.0;

    #if RESCALE_DEBUG > 0
    std::cerr << "Rescaling with query factor: " << query_scaling_factor << "\n";
    #endif

    const size_t l = rescaled_matched_chunks.size();
    for (size_t i = 0; i < l; i++) {
      const MatchedChunk& mc = matched_chunks[i];
      MatchedChunk& rmc =  rescaled_matched_chunks[i];
      int old_query_size = mc.query_chunk.size;
      int new_query_size = query_scaling_factor*old_query_size;
      double old_sizing_score = mc.score.sizing_score;
      double new_sizing_score = sizing_penalty(new_query_size, mc.ref_chunk.size, align_opts);
      rmc.query_chunk.size = new_query_size;
      if (!mc.query_chunk.is_boundary && !mc.ref_chunk.is_boundary) {
        rmc.score.sizing_score = new_sizing_score;
        rescaled_score.sizing_score += new_sizing_score;
        #if RESCALE_DEBUG > 10
        std::cerr << "old_q: " << old_query_size << " new_q: " << new_query_size
                << " old_sizing_score: " << old_sizing_score << " new: " << new_sizing_score << "\n";
        #endif
      }
    }

    total_rescaled_score = rescaled_score.total();

  }


  std::ostream& operator<<(std::ostream& os, const AlignmentHeader&) {
    os << "query_map" << "\t"
       << "ref_map" << "\t"
       << "is_forward" << "\t"
       << "num_matched_chunks" << "\t"
       << "query_misses" << "\t"
       << "ref_misses" << "\t"
       << "query_miss_rate" << "\t"
       << "ref_miss_rate" << "\t"
       << "total_score" << "\t"
       << "total_rescaled_score" << "\t"
       << "m_score" << "\t"
       << "sizing_score" << "\t"
       << "sizing_score_rescaled" << "\t"
       << "query_scaling_factor" << "\t"
       << "chunk_string" << "\n";
    return os;
  }


  void print_alignment(std::ostream& os, const Alignment& aln) {

    os << aln.query_map_data.map_name_ << "\t"
       << aln.ref_map_data.map_name_ << "\t"
       << (aln.is_forward ? "F" : "R") << "\t"
       << aln << "\n";
  }


  
  std::ostream& operator<<(std::ostream& os, const Alignment& aln) {

    os << aln.matched_chunks.size() << "\t"
       << aln.query_misses << "\t"
       << aln.ref_misses << "\t"
       << aln.query_miss_rate << "\t"
       << aln.ref_miss_rate << "\t"
       << aln.total_score << "\t"
       << aln.total_rescaled_score << "\t"
       << aln.m_score << "\t"
       << aln.score.sizing_score << "\t"
       << aln.rescaled_score.sizing_score << "\t"
       << aln.query_scaling_factor << "\t";

    for(auto& mc : aln.rescaled_matched_chunks) {
      os << mc << ";";
    }

    return os;

  }


  const Alignment INVALID_ALIGNMENT;

}
