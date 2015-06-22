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
  struct ScoreMatrixRecord {

    ScoreMatrixRecord() = default;

    ScoreMatrixRecord(const string& query,
                      const string& ref,
                      AlignmentOrientation orientation) :
      query_(query),
      ref_(ref),
      orientation_(orientation),
      row_(-1), col_(-1),
      m_score_(std::numeric_limits<double>::infinity()) { }

    ScoreMatrixRecord(const string& query,
                      const string& ref,
                      AlignmentOrientation orientation,
                      int row,
                      int col,
                      double m_score ) :
      query_(query),
      ref_(ref),
      orientation_(orientation),
      row_(row),
      col_(col),
      m_score_(m_score) { }

    ScoreMatrixRecord(const ScoreMatrixRecord& o) = default;
    ScoreMatrixRecord& operator=(const ScoreMatrixRecord& o) = default;

    bool operator<(const ScoreMatrixRecord& other) const { 
      return m_score_ < other.m_score_;
    }

    string query_;
    string ref_;
    AlignmentOrientation orientation_;
    int row_;
    int col_;
    double m_score_;

  };

  typedef std::vector<ScoreMatrixRecord> ScoreMatrixProfile;
  typedef std::vector<ScoreMatrixProfile> ScoreMatrixProfileVec;

  ScoreMatrixProfile merge_profiles(const ScoreMatrixProfileVec& profiles);

  std::ostream& operator<<(std::ostream& os, const ScoreMatrixRecord& smr);
  std::ostream& operator<<(std::ostream& os, const ScoreMatrixProfile& p);
  std::ostream& operator<<(std::ostream& os, const AlignmentOrientation& o);


}

#endif