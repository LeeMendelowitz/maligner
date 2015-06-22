#include <stdexcept>

#include "score_matrix_profile.h"

namespace maligner_vd {


  std::ostream& operator<<(std::ostream& os, const AlignmentOrientation& o) {

    switch(o) {

      case AlignmentOrientation::RF_QF :
        os << "RF_QF";
        break;

      case AlignmentOrientation::RF_QR :
        os << "RF_QR";
        break;

      case AlignmentOrientation::RR_QF :
        os << "RR_QF";
        break;

      case AlignmentOrientation::RR_QR :
        os << "RR_QR";
        break;
    }

    return os;

  }


  std::ostream& operator<<(std::ostream& os, const ScoreMatrixRecord& smr) {

    os << smr.query_ << "\t"
       << smr.ref_ << "\t"
       << smr.orientation_ << "\t"
       << smr.row_ << "\t"
       << smr.col_ << "\t"
       << smr.m_score_;

    return os;

  }

  std::ostream& operator<<(std::ostream& os,  const ScoreMatrixProfile& p) {
    for(auto& rec : p) {
      os << rec << "\n";
    }
    return os;
  }

  ScoreMatrixProfile merge_profiles(const ScoreMatrixProfileVec& profiles) {
    
    if(profiles.empty()) {
      return ScoreMatrixProfile();
    }

    // Check that each profile in profiles has the same length.
    const size_t N = profiles.front().size();
    for(auto& p : profiles) {
      if (p.size() != N) {
        throw std::runtime_error("ScoreMatrixProfiles have inconsistent sizes.");
      }
    }

    ScoreMatrixProfile ret; // A std::vector of ScoreMatrixRecord's
    ret.reserve(N);

    for (size_t i = 0; i < N; i++) {

      const ScoreMatrixRecord * p_best = &profiles.front()[i];

      const ScoreMatrixProfileVec::const_iterator E = profiles.end();
      for(ScoreMatrixProfileVec::const_iterator p_iter = profiles.begin() + 1;
          p_iter < E;
          p_iter++) {
        
        const ScoreMatrixProfile& profile = *p_iter;
        const ScoreMatrixRecord * p_rec = &profile[i];
        if (p_rec->m_score_ < p_best->m_score_) {
          p_best = p_rec;
        }
        
      }

      ret.push_back(*p_best);

    }

    return ret;

  }

}