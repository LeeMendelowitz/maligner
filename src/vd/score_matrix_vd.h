#ifndef SCORE_MATRIX_DB
#define SCORE_MATRIX_DB

#include <vector>
#include <iostream>
#include <sstream>

#include "ScoreMatrix.h"
#include "alignment.h"
#include "map_wrappers.h"
#include "align.h"
#include "common_math.h"
#include "globals.h"

namespace maligner_vd {

  using maligner_dp::RefMapWrapper;
  using maligner_dp::QueryMapWrapper;
  using maligner_dp::AlignmentVec;
  using maligner_dp::AlignOpts;
  using maligner_dp::Chi2SizingPenalty;
  using maligner_dp::MapData;
  using maligner_dp::Alignment;
  using maligner_dp::ScoreCell;
  using maligner_dp::Constants::INF;

  /////////////////////////////////////////////////////////
  // Define a class for storing ScoreMatrices
  // We store 4 ScoreMatrices for the for different types of alignment:
  //
  // sm_rf_qf_: Reference forward, query forward
  // sm_rf_qr_: Reference forward, query reverse 
  // sm_rf_qf_: Reference reverse, query forward
  // sm_rf_qr_: Reference reverse, query reverse
  //
  template<typename ScoreMatrixType >
  class RefScoreMatrixVD {

    typedef maligner_dp::AlignTask<ScoreMatrixType, Chi2SizingPenalty> AlignTaskType;

  public:

    RefScoreMatrixVD(const RefMapWrapper& ref) :
      ref_(ref) {};

    RefScoreMatrixVD(RefMapWrapper&& ref) :
      ref_(ref) {};

    const RefMapWrapper& get_ref_map() const { return ref_; }

    void aln_to_forward_ref(const QueryMapWrapper& q, const AlignOpts& ao);
    void aln_to_reverse_ref(const QueryMapWrapper& q, const AlignOpts& ao);
    void get_prefix_scores(size_t row_number, std::vector<double>& scores);
    void get_suffix_scores(size_t row_number, std::vector<double>& scores);
    void assign_prefix_mscores(size_t row_number, double median, double mad);
    void assign_suffix_mscores(size_t row_number, double median, double mad);

    Alignment best_aln_rf_qf() const { return get_best_alignment_from_filled_scorematrix(aln_task_rf_qf_);}
    Alignment best_aln_rf_qr() const { return get_best_alignment_from_filled_scorematrix(aln_task_rf_qr_);}
    Alignment best_aln_rr_qr() const { return get_best_alignment_from_filled_scorematrix(aln_task_rr_qr_);}
    Alignment best_aln_rr_qf() const { return get_best_alignment_from_filled_scorematrix(aln_task_rr_qf_);}

    size_t num_rows() const;
    size_t num_rows_ref_forward() const;
    size_t num_rows_ref_reverse() const;
    size_t num_rows_query_suffix() const;
    size_t num_rows_query_prefix() const;

    size_t num_rows_rf_qf() const { return sm_rf_qf_.getNumRows(); }
    size_t num_rows_rf_qr() const { return sm_rf_qr_.getNumRows(); }
    size_t num_rows_rr_qf() const { return sm_rr_qf_.getNumRows(); }
    size_t num_rows_rr_qr() const { return sm_rr_qr_.getNumRows(); }

    const ScoreMatrixType * get_score_matrix_rf_qf() const { return &sm_rf_qf_; }
    const ScoreMatrixType * get_score_matrix_rf_qr() const { return &sm_rf_qr_; }
    const ScoreMatrixType * get_score_matrix_rr_qf() const { return &sm_rr_qf_; }
    const ScoreMatrixType * get_score_matrix_rr_qr() const { return &sm_rr_qr_; }

    void reset();

    bool check_sane() const;
    bool check_sane_query_prefix() const;
    bool check_sane_query_suffix() const;

    uint64_t get_memory_usage() const;
    uint64_t get_memory_capacity() const;


  private:

    RefMapWrapper ref_;

    ScoreMatrixType sm_rf_qf_; // ref forward, query forward
    ScoreMatrixType sm_rf_qr_; // ref forward, query reverse
    ScoreMatrixType sm_rr_qf_; // ref reverse, query forward
    ScoreMatrixType sm_rr_qr_; // ref reverse, query reverse

    AlignTaskType aln_task_rf_qf_;
    AlignTaskType aln_task_rf_qr_;
    AlignTaskType aln_task_rr_qf_;
    AlignTaskType aln_task_rr_qr_;

    AlignmentVec aln_rf_qf_;
    AlignmentVec aln_rf_qr_;
    AlignmentVec aln_rr_qf_;
    AlignmentVec aln_rr_qr_;

  };

  template<typename ScoreMatrixType>
  void RefScoreMatrixVD<ScoreMatrixType>::aln_to_forward_ref(const QueryMapWrapper& query, const AlignOpts& align_opts) {
    

    aln_task_rf_qf_ = AlignTaskType(
      const_cast<MapData*>(&query.map_data_),
      const_cast<MapData*>(&ref_.map_data_),
      &query.get_frags(),
      &ref_.get_frags(), 
      &query.get_partial_sums_forward(),
      &ref_.get_partial_sums(),
      &ref_.sd_inv_,
      &query.ix_to_locs_,
      &ref_.ix_to_locs_,
      0, // ref_offset
      &sm_rf_qf_,
      &aln_rf_qf_,
      true, // query_is_forward
      true, // ref_is_forward
      align_opts
    );

    aln_task_rf_qr_ = AlignTaskType(
      const_cast<MapData*>(&query.map_data_),
      const_cast<MapData*>(&ref_.map_data_),
      &query.get_frags_reverse(),
      &ref_.get_frags(), 
      &query.get_partial_sums_reverse(),
      &ref_.get_partial_sums(),
      &ref_.sd_inv_,
      &query.ix_to_locs_,
      &ref_.ix_to_locs_,
      0, // ref_offset
      &sm_rf_qr_,
      &aln_rf_qr_,
      false, // query_is_forward
      true, // ref_is_forward
      align_opts
    );


    fill_score_matrix_using_partials_with_breaks_hardcode_penalty_max_miss(aln_task_rf_qf_);
    fill_score_matrix_using_partials_with_breaks_hardcode_penalty_max_miss(aln_task_rf_qr_);


  }

  template<typename ScoreMatrixType>
  void RefScoreMatrixVD<ScoreMatrixType>::aln_to_reverse_ref(const QueryMapWrapper& query, const AlignOpts& align_opts) {


    aln_task_rr_qf_ = AlignTaskType(
      const_cast<MapData*>(&query.map_data_),
      const_cast<MapData*>(&ref_.map_data_),
      &query.get_frags(),
      &ref_.get_frags_reverse(), 
      &query.get_partial_sums_forward(),
      &ref_.get_partial_sums_reverse(),
      &ref_.sd_inv_reverse_,
      &query.ix_to_locs_,
      &ref_.ix_to_locs_,
      0, // ref_offset
      &sm_rr_qf_,
      &aln_rr_qf_,
      true, // query_is_forward
      false, // ref_is_forward
      align_opts
    );

    aln_task_rr_qr_ = AlignTaskType(
      const_cast<MapData*>(&query.map_data_),
      const_cast<MapData*>(&ref_.map_data_),
      &query.get_frags_reverse(),
      &ref_.get_frags_reverse(), 
      &query.get_partial_sums_reverse(),
      &ref_.get_partial_sums_reverse(),
      &ref_.sd_inv_reverse_,
      &query.ix_to_locs_,
      &ref_.ix_to_locs_,
      0, // ref_offset
      &sm_rr_qr_,
      &aln_rr_qr_,
      false, // query_is_forward
      false, // ref_is_forward
      align_opts
    );


    fill_score_matrix_using_partials_with_breaks_hardcode_penalty_max_miss(aln_task_rr_qf_);
    fill_score_matrix_using_partials_with_breaks_hardcode_penalty_max_miss(aln_task_rr_qr_);

  }

  template<typename ScoreMatrixType>
  void RefScoreMatrixVD<ScoreMatrixType>::reset() {
    aln_rf_qf_.clear();
    aln_rf_qr_.clear();
    aln_rr_qf_.clear();
    aln_rr_qr_.clear();
  }


  template<typename ScoreMatrixType>
  void RefScoreMatrixVD<ScoreMatrixType>::get_prefix_scores(size_t row_number,
    std::vector<double>& scores) {

    check_sane_query_prefix();

    {
      const size_t num_rows = sm_rf_qf_.getNumRows();
      const size_t num_cols = sm_rf_qf_.getNumCols();
      for(size_t j = 0; j < num_cols; j++) {
        const ScoreCell * pCell = sm_rf_qf_.getCell(row_number, j);
        if (pCell->score_ != -INF) {
          scores.push_back(pCell->score_);
        }
      }
    }

    {
      const size_t num_rows = sm_rr_qf_.getNumRows();
      const size_t num_cols = sm_rr_qf_.getNumCols();
      for(size_t j = 0; j < num_cols; j++) {
        const ScoreCell * pCell = sm_rr_qf_.getCell(row_number, j);
        if (pCell->score_ != -INF) {
          scores.push_back(pCell->score_);
        }
      }
    }

  }


  template<typename ScoreMatrixType>
  void RefScoreMatrixVD<ScoreMatrixType>::get_suffix_scores(size_t row_number,
    std::vector<double>& scores) {

    check_sane_query_suffix();

    {
      const size_t num_rows = sm_rf_qr_.getNumRows();
      const size_t num_cols = sm_rf_qr_.getNumCols();
      for(size_t j = 0; j < num_cols; j++) {
        const ScoreCell * pCell = sm_rf_qr_.getCell(row_number, j);
        if (pCell->score_ != -INF) {
          scores.push_back(pCell->score_);
        }
      }
    }   

    {
      const size_t num_rows = sm_rr_qr_.getNumRows();
      const size_t num_cols = sm_rr_qr_.getNumCols();
      for(size_t j = 0; j < num_cols; j++) {
        const ScoreCell * pCell = sm_rr_qr_.getCell(row_number, j);
        if (pCell->score_ != -INF) {
          scores.push_back(pCell->score_);
        }
      }
    }

  }


  template<typename ScoreMatrixType>  
  bool RefScoreMatrixVD<ScoreMatrixType>::check_sane() const {

    const size_t M = sm_rf_qf_.getNumRows();
    const size_t N = sm_rf_qf_.getNumCols();

    bool rows_okay = ( M == sm_rf_qr_.getNumRows() &&
                  M == sm_rf_qf_.getNumRows() &&
                  M == sm_rr_qf_.getNumRows() &&
                  M == sm_rr_qr_.getNumRows());


    if (!rows_okay) {
      std::ostringstream oss;
      oss << "Num rows do not match across all score matrices: " 
         << sm_rf_qr_.getNumRows() << " "
         << sm_rf_qf_.getNumRows() << " "
         << sm_rr_qf_.getNumRows() << " "
         << sm_rr_qr_.getNumRows();

      throw std::runtime_error(oss.str());
    }

    bool cols_okay = ( N == sm_rf_qr_.getNumCols() &&
                  N == sm_rf_qf_.getNumCols() &&
                  N == sm_rr_qf_.getNumCols() &&
                  N == sm_rr_qr_.getNumCols());

    if (!cols_okay) {
      throw std::runtime_error("Num cols do not match");
    }

    return true;

  }

  template<typename ScoreMatrixType>  
  bool RefScoreMatrixVD<ScoreMatrixType>::check_sane_query_prefix() const {

    const size_t M = sm_rf_qf_.getNumRows();
    const size_t N = sm_rf_qf_.getNumCols();

    bool rows_okay = ( M == sm_rf_qf_.getNumRows() &&
                  M == sm_rf_qf_.getNumRows());


    if (!rows_okay) {
      std::ostringstream oss;
      oss << "Num rows do not match: for query prefix: " 
         << sm_rf_qf_.getNumRows() << " "
         << sm_rr_qf_.getNumRows();

      throw std::runtime_error(oss.str());
    }

    bool cols_okay = ( N == sm_rf_qf_.getNumCols() &&
                  N == sm_rr_qf_.getNumCols());

    if (!cols_okay) {
      throw std::runtime_error("Num cols do not match");
    }

    return true;

  }

  template<typename ScoreMatrixType>  
  bool RefScoreMatrixVD<ScoreMatrixType>::check_sane_query_suffix() const {

    const size_t M = sm_rf_qr_.getNumRows();
    const size_t N = sm_rf_qr_.getNumCols();

    bool rows_okay = ( M == sm_rf_qr_.getNumRows() &&
                  M == sm_rf_qr_.getNumRows());


    if (!rows_okay) {
      std::ostringstream oss;
      oss << "Num rows do not match for query suffix: " 
         << sm_rf_qr_.getNumRows() << " "
         << sm_rr_qr_.getNumRows();

      throw std::runtime_error(oss.str());
    }

    bool cols_okay = ( N == sm_rf_qr_.getNumCols() &&
                  N == sm_rr_qr_.getNumCols());

    if (!cols_okay) {
      throw std::runtime_error("Num cols do not match");
    }

    return true;

  }
  
  template<typename ScoreMatrixType>  
  uint64_t RefScoreMatrixVD<ScoreMatrixType>::get_memory_capacity() const {
      return sm_rf_qf_.getMemoryCapacity() + 
             sm_rf_qr_.getMemoryCapacity() +
             sm_rr_qf_.getMemoryCapacity() +
             sm_rr_qr_.getMemoryCapacity();
  }

  template<typename ScoreMatrixType>
  uint64_t RefScoreMatrixVD<ScoreMatrixType>::get_memory_usage() const {
    return sm_rf_qf_.getMemoryUsage() + 
           sm_rf_qr_.getMemoryUsage() +
           sm_rr_qf_.getMemoryUsage() +
           sm_rr_qr_.getMemoryUsage();
  }


  template<typename ScoreMatrixType>
  size_t RefScoreMatrixVD<ScoreMatrixType>::num_rows() const {
    check_sane();
    return sm_rf_qf_.getNumRows();
  }

  template<typename ScoreMatrixType>
  size_t RefScoreMatrixVD<ScoreMatrixType>::num_rows_query_prefix() const {
    check_sane_query_prefix();
    return sm_rf_qf_.getNumRows();
  }

  template<typename ScoreMatrixType>
  size_t RefScoreMatrixVD<ScoreMatrixType>::num_rows_query_suffix() const {
    check_sane_query_suffix();
    return sm_rf_qr_.getNumRows();
  }

  template<typename ScoreMatrixType>
  size_t RefScoreMatrixVD<ScoreMatrixType>::num_rows_ref_forward() const {

    const size_t num_rows = sm_rf_qf_.getNumRows();
    if(num_rows != sm_rf_qr_.getNumRows()) {
      throw std::runtime_error("Number of rows ref forward do not match!");
    }
    return num_rows;

  }

  template<typename ScoreMatrixType>
  size_t RefScoreMatrixVD<ScoreMatrixType>::num_rows_ref_reverse() const {

    const size_t num_rows = sm_rr_qf_.getNumRows();

    if(num_rows != sm_rr_qr_.getNumRows()) {
      throw std::runtime_error("Number of rows ref forward do not match!");
    }

    return num_rows;

  }


  template<typename ScoreMatrixType>
  void RefScoreMatrixVD<ScoreMatrixType>::assign_prefix_mscores(size_t row_number, double median, double mad) {

    const size_t num_rows = num_rows_query_prefix();
    const size_t num_cols = sm_rf_qf_.getNumCols();

    if (row_number >= num_rows) {
      throw std::runtime_error("Invalid row_number in assign_prefix_mscores. Exceeds number of rows.");
    }

    ScoreCell * p_cell;
    for(size_t j = 1; j < num_cols; j++) {

      p_cell = sm_rf_qf_.getCell(row_number, j);
      if(p_cell->is_valid()) {
        p_cell->m_score_ = (p_cell->score_ - median)/mad;
      }

      p_cell = sm_rr_qf_.getCell(row_number, j);
      if(p_cell->is_valid()) {
        p_cell->m_score_ = (p_cell->score_ - median)/mad;
      }

    }


  }

  
  template<typename ScoreMatrixType>
  void RefScoreMatrixVD<ScoreMatrixType>::assign_suffix_mscores(size_t row_number, double median, double mad) {

    const size_t num_rows = num_rows_query_suffix();
    const size_t num_cols = sm_rf_qr_.getNumCols();

    if (row_number >= num_rows) {
      throw std::runtime_error("Invalid row_number in assign_suffix_mscores. Exceeds number of rows.");
    }

    ScoreCell * p_cell;
    for(size_t j = 1; j < num_cols; j++) {

      p_cell = sm_rf_qr_.getCell(row_number, j);
      if(p_cell->is_valid()) {
        p_cell->m_score_ = (p_cell->score_ - median)/mad;
      }

      p_cell = sm_rr_qr_.getCell(row_number, j);
      if(p_cell->is_valid()) {
        p_cell->m_score_ = (p_cell->score_ - median)/mad;
      }

    }

  }


}

#endif