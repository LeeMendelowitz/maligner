#ifndef MALIGNER_SCORER_H
#define MALIGNER_SCORER_H

#include "matched_chunk.h"

namespace maligner_common {

  using maligner_dp::Score;
  using maligner_dp::MatchedChunk;
  using maligner_dp::Chunk;

  class Scorer {



    public:

      Scorer(double query_miss_penalty, double ref_miss_penalty, double min_sd, double sd_rate) :
        q_(query_miss_penalty),
        r_(ref_miss_penalty),
        min_sd_(min_sd),
        sd_rate_(sd_rate) { }


      Score compute_score(const Chunk& query_chunk, const Chunk& ref_chunk) const {

        const int num_query_misses = query_chunk.num_misses();
        const int num_ref_misses =  ref_chunk.num_misses();
        const int delta =  query_chunk.size - ref_chunk.size;

        double sd = sd_rate_ * ref_chunk.size;
        if (sd < min_sd_) { sd = min_sd_; }
        double sizing_error = delta/sd;
        sizing_error = sizing_error * sizing_error;

        double query_miss_score = num_query_misses*q_;
        double ref_miss_score = num_ref_misses*r_;


        return Score(query_miss_score, ref_miss_score, sizing_error);

      }


      Score compute_score(const MatchedChunk& mc) const {
        return compute_score(mc.query_chunk, mc.ref_chunk);
      }


      double q_;
      double r_;
      double min_sd_;
      double sd_rate_;


  };
} 

#endif