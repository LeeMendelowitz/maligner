#ifndef SCORE_MATRIX_DB
#define SCORE_MATRIX_DB

#include <vector>
#include <iostream>
#include <sstream>
#include <string>

#include "ScoreMatrix.h"
#include "alignment.h"
#include "map_wrappers.h"
#include "align.h"
#include "common_math.h"
#include "globals.h"
#include "score_matrix_profile.h"

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

    ScoreMatrixProfile get_score_matrix_profile_rf_qf(const string& query) const;
    ScoreMatrixProfile get_score_matrix_profile_rf_qr(const string& query) const;
    ScoreMatrixProfile get_score_matrix_profile_rr_qf(const string& query) const;
    ScoreMatrixProfile get_score_matrix_profile_rr_qr(const string& query) const;

    void get_prefix_score_matrix_row_profile(size_t row_number, 
      const string& query,
      bool allow_overlaps,
      size_t max_records,
      ScoreMatrixProfile& vec) const;

    void get_suffix_score_matrix_row_profile(size_t row_number, 
      const string& query,
      bool allow_overlaps,
      size_t max_records,      
      ScoreMatrixProfile& vec) const;

    void get_score_matrix_row_profile(
      const ScoreMatrixType& sm,
      const size_t row_number,
      const string& query,
      AlignmentOrientation orientation,
      bool allow_overlaps,
      size_t max_records,
      ScoreMatrixProfile& recs) const;

    void reset();

    bool check_sane() const;
    bool check_sane_query_prefix() const;
    bool check_sane_query_suffix() const;

    uint64_t get_memory_usage() const;
    uint64_t get_memory_capacity() const;

    void print_filled_by_row() const;


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

    ScoreMatrixProfile get_score_matrix_profile(
      const ScoreMatrixType& sm,
      const string& query,
      AlignmentOrientation orientation) const;

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
        if (pCell->is_valid()) {
          scores.push_back(pCell->score_);
        }
      }
    }

    {
      const size_t num_rows = sm_rr_qf_.getNumRows();
      const size_t num_cols = sm_rr_qf_.getNumCols();
      for(size_t j = 0; j < num_cols; j++) {
        const ScoreCell * pCell = sm_rr_qf_.getCell(row_number, j);
        if (pCell->is_valid()) {
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
        if (pCell->is_valid()) {
          scores.push_back(pCell->score_);
        }
      }
    }   

    {
      const size_t num_rows = sm_rr_qr_.getNumRows();
      const size_t num_cols = sm_rr_qr_.getNumCols();
      for(size_t j = 0; j < num_cols; j++) {
        const ScoreCell * pCell = sm_rr_qr_.getCell(row_number, j);
        if (pCell->is_valid()) {
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

  template<typename ScoreMatrixType>
  void RefScoreMatrixVD<ScoreMatrixType>::print_filled_by_row() const {

    std::vector<double> filled_by_row;

    filled_by_row = sm_rf_qf_.percentFilledByRow();
    std::cerr << ref_.get_name() << " " << std::fixed << std::setprecision(3);
    for(auto v: filled_by_row) {
      std::cerr << v << " ";
    }
    std::cerr << "\n";

    // ScoreMatrixType sm_rf_qf_; // ref forward, query forward
    // ScoreMatrixType sm_rf_qr_; // ref forward, query reverse
    // ScoreMatrixType sm_rr_qf_; // ref reverse, query forward
    // ScoreMatrixType sm_rr_qr_; // ref reverse, query reverse

  }


  template<typename ScoreMatrixType>
  ScoreMatrixProfile RefScoreMatrixVD<ScoreMatrixType>::get_score_matrix_profile_rf_qf(const string& query) const {
    return get_score_matrix_profile(sm_rf_qf_, query, AlignmentOrientation::RF_QF);
  }

  template<typename ScoreMatrixType>
  ScoreMatrixProfile RefScoreMatrixVD<ScoreMatrixType>::get_score_matrix_profile_rf_qr(const string& query) const {
    return get_score_matrix_profile(sm_rf_qr_, query, AlignmentOrientation::RF_QR);
  }


  template<typename ScoreMatrixType>
  ScoreMatrixProfile RefScoreMatrixVD<ScoreMatrixType>::get_score_matrix_profile_rr_qf(const string& query) const {
    return get_score_matrix_profile(sm_rr_qf_, query, AlignmentOrientation::RR_QF);
  }

  template<typename ScoreMatrixType>
  ScoreMatrixProfile RefScoreMatrixVD<ScoreMatrixType>::get_score_matrix_profile_rr_qr(const string& query) const {
    return get_score_matrix_profile(sm_rr_qr_, query, AlignmentOrientation::RR_QR);
  }   


  template<typename ScoreMatrixType>
  ScoreMatrixProfile RefScoreMatrixVD<ScoreMatrixType>::get_score_matrix_profile(
    const ScoreMatrixType& sm,
    const string& query,
    AlignmentOrientation orientation ) const {
    
    

    const size_t num_rows = sm.getNumRows();
    const size_t num_cols = sm.getNumCols();
    const string ref = ref_.get_name();

    ScoreMatrixProfile ret;
    ret.reserve(num_rows);

    for(size_t i = 0; i < num_rows; i++) {

      ScoreMatrixRecord rec(query, ref, orientation);
      rec.m_score_ = std::numeric_limits<double>::infinity();

      for(size_t j = 0; j < num_cols; j++) {

        const ScoreCell* p_cell = sm.getCell(i, j);

        if (p_cell->m_score_ < rec.m_score_) {

          rec.score_ = p_cell->score_;
          rec.m_score_ = p_cell->m_score_;
          rec.row_ = i;
          rec.col_ = j;
          rec.col_start_ = p_cell->ref_start_;

        }

      }

      ret.push_back(std::move(rec));

    }

    return ret;

  }


  template<typename ScoreMatrixType>
  void RefScoreMatrixVD<ScoreMatrixType>::get_score_matrix_row_profile(
    const ScoreMatrixType& sm,
    const size_t row_number,
    const string& query,
    AlignmentOrientation orientation,
    bool allow_overlaps,
    size_t max_records,
    ScoreMatrixProfile& recs) const {
    
    using std::copy;
    using std::make_move_iterator;
    using std::back_inserter;


    const size_t num_rows = sm.getNumRows();
    const size_t num_cols = sm.getNumCols();

    if(row_number >= num_rows) {
      throw std::runtime_error("Invalid row number.");
    }

    const string ref = ref_.get_name();
    const size_t num_ref_frags = ref_.num_frags();

    ScoreMatrixProfile my_recs;
    my_recs.reserve(max_records);

    BitCover bit_cover(num_cols);

    // Extract the score cells from the row of interest
    std::vector<const ScoreCell*> score_cells;
    score_cells.reserve(num_cols);
    for(size_t j = 0; j < num_cols; j++) {

      const ScoreCell* p_cell = sm.getCell(row_number, j);

      if(!p_cell->is_valid()) {
        continue;
      }

      score_cells.push_back(p_cell);

    }

    // Sort the score cells by score
    std::sort(score_cells.begin(), score_cells.end(), maligner_dp::ScoreCellPointerCmp());

    // Select the best scoring, non-overlapping records

    const size_t N = score_cells.size();
    for(size_t i = 0; i < N; i++) {

      const ScoreCell * p_cell = score_cells[i];

      ScoreMatrixRecord rec(query, ref, orientation,
        row_number, p_cell->r_, p_cell->ref_start_, p_cell->score_, p_cell->m_score_);

      if( !allow_overlaps ) {
        int start = rec.get_ref_start(num_ref_frags);
        int end = rec.get_ref_end(num_ref_frags);
        if(bit_cover.is_covered(start, end))
          continue;
        bit_cover.cover_safe(start, end);
      }

      my_recs.push_back(std::move(rec));

      if(my_recs.size() == max_records) break;

    }

    copy(make_move_iterator(my_recs.begin()),
      make_move_iterator(my_recs.end()),
      back_inserter(recs));

    my_recs.clear();

  }

  template<typename ScoreMatrixType>
  void RefScoreMatrixVD<ScoreMatrixType>::get_prefix_score_matrix_row_profile(size_t row_number, 
      const string& query,
      bool allow_overlaps,
      size_t max_records,
      ScoreMatrixProfile& vec) const {

    get_score_matrix_row_profile(sm_rf_qf_, row_number, query, AlignmentOrientation::RF_QF, allow_overlaps, max_records, vec);
    get_score_matrix_row_profile(sm_rr_qf_, row_number, query, AlignmentOrientation::RR_QF, allow_overlaps, max_records, vec);

  }

  template<typename ScoreMatrixType>
  void RefScoreMatrixVD<ScoreMatrixType>::get_suffix_score_matrix_row_profile(size_t row_number, 
      const string& query,
      bool allow_overlaps,
      size_t max_records,
      ScoreMatrixProfile& vec) const {

    get_score_matrix_row_profile(sm_rf_qf_, row_number, query, AlignmentOrientation::RF_QF, allow_overlaps, max_records, vec);
    get_score_matrix_row_profile(sm_rr_qf_, row_number, query, AlignmentOrientation::RR_QF, allow_overlaps, max_records, vec);

  }
}


#endif