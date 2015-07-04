#ifndef SCORE_MATRIX_VD_DB
#define SCORE_MATRIX_VD_DB

#include <vector>
#include <iostream>
#include <string>
#include <algorithm>

#include "score_matrix_vd.h"
#include "score_matrix_profile.h"

namespace maligner_vd {

  using maligner_dp::RefMapWrapper;
  using maligner_dp::QueryMapWrapper;
  using maligner_dp::AlignmentVec;
  using maligner_dp::AlignOpts;
  using maligner_dp::Chi2SizingPenalty;
  using maligner_dp::MapData;
  using std::vector;
  using std::string;

  // Define a class for storing ScoreMatrices
  template<typename ScoreMatrixVDType >
  class RefScoreMatrixVDDB {


    typedef std::vector<ScoreMatrixVDType> ScoreMatrixVDVec;
    typedef std::vector<double> DoubleVec;
    typedef std::vector<DoubleVec> DoubleVecVec;

  public:

    void add_ref_map(RefMapWrapper& rmw) { sm_vec_.push_back(rmw); _compute_max_scores(); }
    void add_ref_map(RefMapWrapper&& rmw) { sm_vec_.push_back(rmw); _compute_max_scores(); }

    void aln_to_forward_refs(const QueryMapWrapper& query, const AlignOpts& align_opts);
    void aln_to_reverse_refs(const QueryMapWrapper& q, const AlignOpts& ao);

    // Compute mscores for prefix alignment and assign the m_score to the ScoreCell's in the ScoreMatrixes.
    void compute_query_prefix_mscores(size_t max_samples, double min_mad, const std::string& query);

    // Compute mscores for suffix alignment and assign the m_score to the ScoreCell's in the ScoreMatrixes.
    void compute_query_suffix_mscores(size_t max_samples, double min_mad, const std::string& query);
    
    // Compute and assign MScores for prefix and suffix alignment.
    void compute_mscores(size_t max_samples, double min_mad, const std::string& query);
   
    size_t num_maps() const { return sm_vec_.size(); }
    size_t size() const { return sm_vec_.size(); }

    const ScoreMatrixVDVec& get_score_matrix_vec() const { return sm_vec_; }

    ScoreMatrixProfile get_score_matrix_profile_rf_qf(const string& query);
    ScoreMatrixProfile get_score_matrix_profile_rf_qr(const string& query);
    ScoreMatrixProfile get_score_matrix_profile_rr_qf(const string& query);
    ScoreMatrixProfile get_score_matrix_profile_rr_qr(const string& query);

    AlignmentVec get_best_alignments_rf_qf(size_t max_alignments, int min_aln_chunks) const;
    AlignmentVec get_best_alignments_rf_qr(size_t max_alignments, int min_aln_chunks) const;
    AlignmentVec get_best_alignments_rr_qf(size_t max_alignments, int min_aln_chunks) const;
    AlignmentVec get_best_alignments_rr_qr(size_t max_alignments, int min_aln_chunks) const;

    template<typename MatrixGetter>
    ScoreMatrixProfile get_score_matrix_profile_helper(MatrixGetter g);



    template<typename AlignmentGetter>
    AlignmentVec get_alignments_helper(AlignmentGetter g,
      size_t max_alignments) const;

  private:


    // Compute the number of scores we may need to store in the forward + reverse reference.
    // For computing m-scores.
    void _compute_max_scores() {

      _max_num_scores = 0;
      for(const auto& sm : sm_vec_) _max_num_scores += sm.get_ref_map().num_frags_total() + 1;

      // Need to handle scores for forward and reverse alignments, hence double.
      _max_num_scores = 2*_max_num_scores; 
    }



    // Check that all score matrices for query suffix the same number of rows.
    bool check_query_suffix_sane() const;

    // Check that all score matrices for query prefix alignment have the same number of rows.
    bool check_query_prefix_sane() const;

    size_t num_rows_query_prefix() const;
    size_t num_rows_query_suffix() const;

    ScoreMatrixVDVec sm_vec_;
    size_t  _max_num_scores;

  };

  template<typename ScoreMatrixVDType>
  void RefScoreMatrixVDDB<ScoreMatrixVDType>::aln_to_forward_refs(const QueryMapWrapper& query, const AlignOpts& align_opts) {
    
    for(auto& m : sm_vec_) {
      
      std::cerr << "aligning to " << m.get_ref_map().get_name() << " forward ref\n";

      m.aln_to_forward_ref(query, align_opts);
      
      std::cerr << "num rows: " << m.num_rows_ref_forward() << " "
          << "usage (bytes): " << m.get_memory_usage() << " "
          << "capacity (bytes): " << m.get_memory_capacity() << "\n";
    }

  }

  template<typename ScoreMatrixVDType>
  void RefScoreMatrixVDDB<ScoreMatrixVDType>::aln_to_reverse_refs(const QueryMapWrapper& query, const AlignOpts& align_opts) {
   
    for(auto& m : sm_vec_) {

      std::cerr << "aligning to " << m.get_ref_map().get_name() << " reverse ref\n";

      m.aln_to_reverse_ref(query, align_opts);

      std::cerr << "num rows: " << m.num_rows_ref_reverse() << " "
          << "usage (bytes): " << m.get_memory_usage() << " "
          << "capacity (bytes): " << m.get_memory_capacity() << "\n";      
    }

  }

  /////////////////////////////////////////////////////////////////////////////
  // Compute m_scores for prefix alignment of the query.
  // This corresponds to sm_rf_qf_ and sm_rr_qf_ score matrixes, as
  // these align query in forward direction, so DP extending prefix alignments
  // of query.
  template<typename ScoreMatrixVDType>
  void RefScoreMatrixVDDB<ScoreMatrixVDType>::compute_query_prefix_mscores(size_t max_samples, double min_mad, const std::string& query) {
   
    // We must gather scores all all alignments by row.
    if (sm_vec_.size() == 0) return;

    std::vector<ScoreMatrixRecord> recs;
    recs.reserve(_max_num_scores);

    const size_t num_rows = num_rows_query_prefix();

    const bool allow_overlaps = false;
    for(size_t row_num = 0; row_num < num_rows; row_num++) {

      recs.clear();
      recs.reserve(_max_num_scores);

      // Compute mad & median of scores for this row.

      for(const auto& sm: sm_vec_) {
        sm.get_prefix_score_matrix_row_profile(row_num, query, allow_overlaps, max_samples, recs);
      }

      std::sort(recs.begin(), recs.end(), ScoreMatrixRecordScoreCmp());

      // Extract scores
      std::vector<double> scores(recs.size());
      std::transform(recs.begin(), recs.end(), scores.begin(), 
          [](const ScoreMatrixRecord& rec) {return rec.score_; });

      // For information purposes summarize all of the row scores.
      // double med = median_sorted(scores);
      // double md = mad(scores, med);
      // double min_value = scores.front();
      // double max_value = scores.back();

      if(recs.size() > max_samples) {
        recs.resize(max_samples);
        scores.resize(max_samples);
      }

      double med_top = median_sorted(scores);
      double md_top = std::max(mad(scores, med_top), min_mad);

      // std::cerr << "prefix row: " << row_num << " "
      //           << "n: " << recs.size() << " "
      //           << "min: " << min_value << " "
      //           << "max: " << max_value << " "
      //           << "median: " << med << " "
      //           << "mad: " << md << " "
      //           << "median (top): " << med_top << " "
      //           << "mad (top): " << md_top << "\n";

      // Print out the records
      // for(auto& rec: recs) {
      //   std::cerr << rec << "\n";
      // }

      // Assign m-scores for this row.
      for(auto& sm: sm_vec_) {
        sm.assign_prefix_mscores(row_num, med_top, md_top);
      }

    } 
  }

  /////////////////////////////////////////////////////////////////////////////
  // Compute m_scores for suffix alignment of the query.
  // This corresponds to sm_rf_qr_ and sm_rr_qr_ score matrixes, as
  // these align query in reverse direction, so DP extending suffix alignments
  // of query.
  template<typename ScoreMatrixVDType>
  void RefScoreMatrixVDDB<ScoreMatrixVDType>::compute_query_suffix_mscores(size_t max_samples,
    double min_mad, const std::string& query) {
   
    // Compute m_scores for suffix alignment of the query.
    // This corresponds to sm_rf_qr_ and sm_rr_qr_ score matrixes, as
    // these align query in reverse direction, so DP extending suffix alignments
    // of query.

    // We must gather scores all all alignments by row.
    if (sm_vec_.size() == 0) return;

    std::vector<ScoreMatrixRecord> recs;
    recs.reserve(_max_num_scores);

    const size_t num_rows = num_rows_query_suffix();
    const bool allow_overlaps = false;

    for(size_t row_num = 0; row_num < num_rows; row_num++) {

      // Compute mad & median of scores for this row.
      recs.clear();
      recs.reserve(_max_num_scores);
      
      for(const auto& sm: sm_vec_) {
        sm.get_suffix_score_matrix_row_profile(row_num, query, allow_overlaps, max_samples, recs);
      }

      std::sort(recs.begin(), recs.end(), ScoreMatrixRecordScoreCmp());

      // Extract scores
      std::vector<double> scores(recs.size());
      std::transform(recs.begin(), recs.end(), scores.begin(), 
          [](const ScoreMatrixRecord& rec) {return rec.score_; });

      // For information purposes summarize all of the row scores.
      // double med = median_sorted(scores);
      // double md = mad(scores, med);
      // double min_value = scores.front();
      // double max_value = scores.back();

      if(recs.size() > max_samples) {
        recs.resize(max_samples);
        scores.resize(max_samples);
      }

      double med_top = median_sorted(scores);
      double md_top = std::max(mad(scores, med_top), min_mad);

      // std::cerr << "prefix row: " << row_num << " "
      //           << "n: " << recs.size() << " "
      //           << "min: " << min_value << " "
      //           << "max: " << max_value << " "
      //           << "median: " << med << " "
      //           << "mad: " << md << " "
      //           << "median (top): " << med_top << " "
      //           << "mad (top): " << md_top << "\n";

      // Print out the records
      // for(auto& rec: recs) {
      //   std::cerr << rec << "\n";
      // }

      // Assign m-scores for this row.
      for(auto& sm: sm_vec_) {
        sm.assign_suffix_mscores(row_num, med_top, md_top);
      }

    } 
  }

  template<typename ScoreMatrixVDType>
  void RefScoreMatrixVDDB<ScoreMatrixVDType>::compute_mscores(size_t max_samples,
    double min_mad, const std::string& query) {

    compute_query_suffix_mscores(max_samples, min_mad, query);
    compute_query_suffix_mscores(max_samples, min_mad, query);
   
  }

  template<typename ScoreMatrixVDType>
  bool RefScoreMatrixVDDB<ScoreMatrixVDType>::check_query_suffix_sane() const {

    if (sm_vec_.empty()) return true;

    const size_t num_rows = sm_vec_.front().num_rows_query_suffix();
    for(const auto& sm : sm_vec_) {
      if(num_rows != sm.num_rows_query_suffix()) {
        throw std::runtime_error("Num rows do not match for suffix alignment.");
      }
    }

    return true;

  }

  template<typename ScoreMatrixVDType>
  bool RefScoreMatrixVDDB<ScoreMatrixVDType>::check_query_prefix_sane() const {

    if (sm_vec_.empty()) return true;

    const size_t num_rows = sm_vec_.front().num_rows_query_prefix();
    for(const auto& sm : sm_vec_) {
      if(num_rows != sm.num_rows_query_prefix()) {
        throw std::runtime_error("Num rows do not match for prefix alignment.");
      }
    }

    return true;   

  }

  template<typename ScoreMatrixVDType>
  size_t RefScoreMatrixVDDB<ScoreMatrixVDType>::num_rows_query_prefix() const {

    if (sm_vec_.empty()) return 0;
    check_query_prefix_sane();
    return sm_vec_.front().num_rows_query_prefix();

  }

  template<typename ScoreMatrixVDType>
  size_t RefScoreMatrixVDDB<ScoreMatrixVDType>::num_rows_query_suffix() const {

    if (sm_vec_.empty()) return 0;
    check_query_suffix_sane();
    return sm_vec_.front().num_rows_query_suffix();
    
  }

  template<typename ScoreMatrixVDType>
  template<typename MatrixGetter>
  ScoreMatrixProfile RefScoreMatrixVDDB<ScoreMatrixVDType>::get_score_matrix_profile_helper(MatrixGetter g) {

    // Get the profile for each underlying matix
    ScoreMatrixProfileVec profiles(sm_vec_.size());
    std::transform(sm_vec_.begin(), sm_vec_.end(), profiles.begin(), g);

    // Merge the profiles
    return merge_profiles(profiles);

  }

  template<typename ScoreMatrixVDType>
  template<typename AlignmentGetter>
  AlignmentVec RefScoreMatrixVDDB<ScoreMatrixVDType>::get_alignments_helper(AlignmentGetter g,
    size_t max_alignments) const {

    using maligner_dp::AlignmentMScoreComp;

    // Get the profile for each underlying matix
    AlignmentVec alns;

    for(const auto& sm: sm_vec_) {

      // Get alignments for this score matrix.
      AlignmentVec alns_m = g(sm);
      alns.insert(alns.end(), make_move_iterator(alns_m.begin()), make_move_iterator(alns_m.end()));

    }

    // Sort alignments by m_score
    std::sort(alns.begin(), alns.end(), AlignmentMScoreComp() );

    if (alns.size() > max_alignments) {
      alns.resize(max_alignments);
    }
    
    return alns;

  }

  template<typename ScoreMatrixVDType>
  ScoreMatrixProfile RefScoreMatrixVDDB<ScoreMatrixVDType>::get_score_matrix_profile_rf_qf(const string& query) {

    return get_score_matrix_profile_helper([&query](const ScoreMatrixVDType& sm) {
      return sm.get_score_matrix_profile_rf_qf(query);
    });

  }

  template<typename ScoreMatrixVDType>
  ScoreMatrixProfile RefScoreMatrixVDDB<ScoreMatrixVDType>::get_score_matrix_profile_rf_qr(const string& query) {

    return get_score_matrix_profile_helper([&query](const ScoreMatrixVDType& sm) {
      return sm.get_score_matrix_profile_rf_qr(query);
    });
  }

  
  template<typename ScoreMatrixVDType>
  ScoreMatrixProfile RefScoreMatrixVDDB<ScoreMatrixVDType>::get_score_matrix_profile_rr_qf(const string& query) {

    return get_score_matrix_profile_helper([&query](const ScoreMatrixVDType& sm) {
      return sm.get_score_matrix_profile_rr_qf(query);
    });

  }

  template<typename ScoreMatrixVDType>
  ScoreMatrixProfile RefScoreMatrixVDDB<ScoreMatrixVDType>::get_score_matrix_profile_rr_qr(const string& query) {

    return get_score_matrix_profile_helper([&query](const ScoreMatrixVDType& sm) {
      return sm.get_score_matrix_profile_rr_qr(query);
    });

  }

  template<typename ScoreMatrixVDType>
  AlignmentVec RefScoreMatrixVDDB<ScoreMatrixVDType>::get_best_alignments_rf_qf(size_t max_alignments, int min_aln_chunks) const {

      auto aln_getter = [max_alignments, min_aln_chunks](const ScoreMatrixVDType& sm) {
        return sm.get_best_alignments_rf_qf(max_alignments, min_aln_chunks);
      };

      return get_alignments_helper(aln_getter, max_alignments);
  }

  template<typename ScoreMatrixVDType>
  AlignmentVec RefScoreMatrixVDDB<ScoreMatrixVDType>::get_best_alignments_rf_qr(size_t max_alignments, int min_aln_chunks) const {

      auto aln_getter = [max_alignments, min_aln_chunks](const ScoreMatrixVDType& sm) {
        return sm.get_best_alignments_rf_qr(max_alignments, min_aln_chunks);
      };

      return get_alignments_helper(aln_getter, max_alignments);
  }

  template<typename ScoreMatrixVDType>
  AlignmentVec RefScoreMatrixVDDB<ScoreMatrixVDType>::get_best_alignments_rr_qf(size_t max_alignments, int min_aln_chunks) const {

      auto aln_getter = [max_alignments, min_aln_chunks](const ScoreMatrixVDType& sm) {
        return sm.get_best_alignments_rr_qf(max_alignments, min_aln_chunks);
      };

      return get_alignments_helper(aln_getter, max_alignments);
  }

  template<typename ScoreMatrixVDType>
  AlignmentVec RefScoreMatrixVDDB<ScoreMatrixVDType>::get_best_alignments_rr_qr(size_t max_alignments, int min_aln_chunks) const {

      auto aln_getter = [max_alignments, min_aln_chunks](const ScoreMatrixVDType& sm) {
        return sm.get_best_alignments_rr_qr(max_alignments, min_aln_chunks);
      };

      return get_alignments_helper(aln_getter, max_alignments);
  }

}

#endif