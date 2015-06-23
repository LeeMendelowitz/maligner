#ifndef SCORE_MATRIX_PROFILE_H
#define SCORE_MATRIX_PROFILE_H

#include <string>
#include <ostream>
#include <vector>
#include <limits>

namespace maligner_vd {

  using std::string;

  enum class AlignmentOrientation {RF_QF, RF_QR, RR_QF, RR_QR};

  /////////////////////////////////////////////////////////
  // A record for summarizing a score matrix.
  // Used for collecting diagnostic information about the
  // ScoreMatrix contents.
  struct ScoreMatrixRecord {

    ScoreMatrixRecord() = default;

    ScoreMatrixRecord(const string& query,
                      const string& ref,
                      AlignmentOrientation orientation) :
      query_(query),
      ref_(ref),
      orientation_(orientation),
      row_(-1),
      col_(-1),
      col_start_(-1),
      score_(-std::numeric_limits<double>::infinity()),
      m_score_(std::numeric_limits<double>::infinity()) { }

    ScoreMatrixRecord(const string& query,
                      const string& ref,
                      AlignmentOrientation orientation,
                      int row,
                      int col,
                      int col_start,
                      double score,
                      double m_score ) :
      query_(query),
      ref_(ref),
      orientation_(orientation),
      row_(row),
      col_(col),
      col_start_(col_start),
      score_(score),
      m_score_(m_score) { }

    ScoreMatrixRecord(const ScoreMatrixRecord& o) = default;
    ScoreMatrixRecord& operator=(const ScoreMatrixRecord& o) = default;

    bool operator<(const ScoreMatrixRecord& other) const { 
      return m_score_ < other.m_score_;
    }

    bool ref_is_forward() const {
      return (orientation_ == AlignmentOrientation::RF_QF || orientation_ == AlignmentOrientation::RF_QR);
    }

    // Get the starting reference coordinate, with respect
    // to the forward orientation of the reference
    int get_ref_start(int num_ref_frags) const {
      
      if (ref_is_forward()) {
        return col_start_;
      }
      
      //Reverse
      int ret = num_ref_frags - col_;
      
      // Ensure it's positive (might be negative due to circularization)
      if(ret < 0) {
        return ret + num_ref_frags;
      }

      return ret;
    }

    // Get the ending reference coordinate, with respect
    // to the forward orientation of the reference
    int get_ref_end(int num_ref_frags) const {

      if (ref_is_forward()) {
        return col_;
      }
      
      //Reverse
      int ret = num_ref_frags - col_start_;
      
      // Ensure it's positive (might be negative due to circularization)
      if(ret < 0) {
        return ret + num_ref_frags;
      }

      return ret;
    }

    string query_;
    string ref_;
    AlignmentOrientation orientation_;
    int row_;
    int col_; // The column of the ScoreMatrix where the corresponding alignment ends
    int col_start_; // The column of the ScoreMatrix where the corresponding alignment starts
    double score_;
    double m_score_;

  };

  struct ScoreMatrixRecordScoreCmp {
    bool operator()(const ScoreMatrixRecord r1, const ScoreMatrixRecord r2) {
      return r1.score_ > r2.score_;
    }
  };

  struct ScoreMatrixRecordMScoreCmp {
    bool operator()(const ScoreMatrixRecord r1, const ScoreMatrixRecord r2) {
      return r1.m_score_ > r2.m_score_;
    }
  };


  typedef std::vector<ScoreMatrixRecord> ScoreMatrixProfile;
  typedef std::vector<ScoreMatrixProfile> ScoreMatrixProfileVec;

  ScoreMatrixProfile merge_profiles(const ScoreMatrixProfileVec& profiles);

  // Remove overlapping records in a profile.
  // The profile must include records that are aligned to a *single* reference.
  void remove_overlaps(ScoreMatrixProfile& profile, int num_ref_frags);

  struct ScoreMatrixRecordHeader {};
  std::ostream& operator<<(std::ostream& os, const ScoreMatrixRecord& smr);
  std::ostream& operator<<(std::ostream& os, ScoreMatrixRecordHeader);
  std::ostream& operator<<(std::ostream& os, const ScoreMatrixProfile& p);
  std::ostream& operator<<(std::ostream& os, const AlignmentOrientation& o);


}

#endif