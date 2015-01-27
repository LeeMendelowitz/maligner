#ifndef SCOREMATRIX_H
#define SCOREMATRIX_H

/**********************************************************

This ScoreMatrix will maintain a vector to hold ScoreCell
objects. Ideally a single ScoreMatrix should be used per thread
for the lifetime of the thread.
    
Before the first alignment, a ScoreMatrix should be initialized with ScoreCells, and
ideally given a capacity large enough to handle any alignment problem
it will see.

**************************************************************/
#include <vector>
#include <memory>
#include <iostream>

#include "ScoreCell.h"
#include "types.h"


namespace maligner_dp {


    using std::shared_ptr;

    struct row_order_tag{};
    struct column_order_tag{};

    template<class OrderTag>
    class ScoreMatrix
    {
    public:

        typedef OrderTag order_tag;

        ScoreMatrix(size_t numRows = 0, size_t numCols = 0);
        void resize(size_t numRows, size_t numCols);

        size_t getSize() const { return size_; }
        size_t getCapacity() const { return data_.size(); }
        size_t getNumRows() const { return numRows_; }
        size_t getNumCols() const { return numCols_; }

        ScoreCell * getCell(size_t row, size_t col);
        ScoreCell * getCell(const IntPair& coord);
        ScoreCell * getCell(size_t n);
        ScoreCellVec& getData() { return data_; }

        const ScoreCell * getCell(size_t row, size_t col) const;
        const ScoreCell * getCell(const IntPair& coord) const;
        const ScoreCell * getCell(size_t n) const;
        const ScoreCellVec& getData() const { return data_; }

        // Indicators for whether cell in last row is in play.
        bool cell_in_play(size_t col) const { return last_row_in_play_[col];}
        void mark_cell_in_play(size_t col, bool in_play) { last_row_in_play_[col] = in_play; }
        void reset();

        // Summary Functions
        size_t countFilledCells() const;
        size_t countFilledByRow(size_t row) const;
        double percentFilled() const;
        double getMaxScore() const;
        double getMaxScoreByRow(size_t row) const;
        size_t countColor(ScoreCellColor col) const;

        template<class OtherOrderTag>
        bool operator==(const ScoreMatrix<OtherOrderTag>& sm_o);

    private:
        size_t numRows_;
        size_t numCols_;
        size_t size_;
        ScoreCellVec data_;
        BoolVec last_row_in_play_; // Indicate which cells in last row are in play.

        // Helper functions using tag dispatch
        ScoreCell * _get_cell(size_t row, size_t col, row_order_tag);
        ScoreCell * _get_cell(size_t row, size_t col, column_order_tag);
        const ScoreCell * _get_cell(size_t row, size_t col, row_order_tag) const;
        const ScoreCell * _get_cell(size_t row, size_t col, column_order_tag) const;
        void _reset(row_order_tag);
        void _reset(column_order_tag);

        
    };

    // typedef std::shared_ptr< ScoreMatrix> ScoreMatrixPtr;
    typedef ScoreMatrix<column_order_tag>* ScoreMatrixPtr;
    typedef ScoreMatrix<column_order_tag> ColumnScoreMatrix;
    typedef ScoreMatrix<row_order_tag> RowScoreMatrix;

    ////////////////////////////////////////////////////////////
    template<class OrderTag>
    inline ScoreMatrix<OrderTag>::ScoreMatrix(size_t numRows, size_t numCols) :
        numRows_(numRows),
        numCols_(numCols),
        size_(numRows_*numCols_)
    {
        resize(numRows, numCols);
    }

    template<class OrderTag>
    inline void ScoreMatrix<OrderTag>::reset() {
        _reset(OrderTag());
        for (size_t j = 0; j < numCols_; j++) {
            last_row_in_play_[j] = true;
        }
    }

    template<class OrderTag>
    inline void ScoreMatrix<OrderTag>::_reset(row_order_tag) {
        // Set the coordinates for those cells that are in bounds
        for (size_t i = 0; i < numRows_; i++) {
            for (size_t j = 0; j < numCols_; j++) {
                ScoreCell * pCell = getCell(i, j);
                pCell->reset();
                pCell->q_ = i;
                pCell->r_ = j;
                pCell->backPointer_ = nullptr;
                pCell->score_ = -INF;
                pCell->color_ = ScoreCellColor::WHITE;
            }
        }
    }

    template<class OrderTag>
    inline void ScoreMatrix<OrderTag>::_reset(column_order_tag) {
        // Set the coordinates for those cells that are in bounds
        for (size_t j = 0; j < numCols_; j++) {
            for (size_t i = 0; i < numRows_; i++) {
                ScoreCell * pCell = getCell(i, j);
                pCell->reset();
                pCell->q_ = i;
                pCell->r_ = j;
                pCell->backPointer_ = nullptr;
                pCell->score_ = -INF;
                pCell->color_ = ScoreCellColor::WHITE;
            }
        }
    }


    template<class OrderTag>
    inline void ScoreMatrix<OrderTag>::resize(size_t numRows, size_t numCols) {
        
        numRows_ = numRows;
        numCols_ = numCols;
        size_ = numRows*numCols;
        if (size_ > data_.size())
        {
            data_.resize(size_);
        }

        last_row_in_play_.resize(numCols_, true);

        reset();
    }

    ///////////////////////////////////////////////////////////////////////////////////////
    template<class OrderTag>
    inline ScoreCell * ScoreMatrix<OrderTag>::getCell(size_t row, size_t col) {
        return _get_cell(row, col, OrderTag());
    }

    template<class OrderTag>
    inline ScoreCell * ScoreMatrix<OrderTag>::getCell(const IntPair& coords)
    {
        return _get_cell(coords.first, coords.second, OrderTag());
    }

    template<class OrderTag>
    inline ScoreCell * ScoreMatrix<OrderTag>::getCell(size_t n) {
        return &data_[n];
    }

    template<class OrderTag>
    inline const ScoreCell * ScoreMatrix<OrderTag>::getCell(size_t row, size_t col) const {
        return _get_cell(row, col, OrderTag());
    }

    template<class OrderTag>
    inline const ScoreCell * ScoreMatrix<OrderTag>::getCell(const IntPair& coords) const
    {
        return _get_cell(coords.first, coords.second, OrderTag());
    }

    template<class OrderTag>
    inline const ScoreCell * ScoreMatrix<OrderTag>::getCell(size_t n) const {
        return &data_[n];
    }
    //////////////////////////////////////////////////////////////////////////////////////////


    //////////////////////////////////////////////////////////
    // These functions actually retrieve the cell. They are called using tag dispatch based on the orientation
    // of the matrix (column major ordered or row major ordered)
    template<class OrderTag>
    inline const ScoreCell* ScoreMatrix<OrderTag>::_get_cell(size_t row, size_t col, row_order_tag) const {
        return &data_[row*numCols_ + col];
    }

    template<class OrderTag>
    inline ScoreCell* ScoreMatrix<OrderTag>::_get_cell(size_t row, size_t col, row_order_tag) {
        return &data_[row*numCols_ + col];
    }

    template<class OrderTag>
    inline const ScoreCell* ScoreMatrix<OrderTag>::_get_cell(size_t row, size_t col, column_order_tag) const {
        return &data_[col*numRows_ + row]; 
    }

    template<class OrderTag>
    inline ScoreCell* ScoreMatrix<OrderTag>::_get_cell(size_t row, size_t col, column_order_tag) {
        return &data_[col*numRows_ + row]; 
    }
    /////////////////////////////////////////////////////////

   template<class OrderTag>
   size_t ScoreMatrix<OrderTag>::countFilledCells() const {
        size_t filled = 0;
        for (int j = 0; j < getNumCols(); j++) {
            for (int i = 0; i < getNumRows(); i++) {
                const ScoreCell* pCell = getCell(i, j);
                if (pCell->score_ > -INF) {
                    filled++;
                }
            }
        }
        return filled;
    }

    template<class OrderTag>
    size_t ScoreMatrix<OrderTag>::countFilledByRow(size_t row) const {
        
        if (row >= numRows_) {
            return 0;
        }

        size_t count = 0;

        for (size_t j = 0; j < numCols_; j++) {
            const ScoreCell * pCell = getCell(row, j);
            if ( pCell->score_ > -INF) {
                count++;
            }
        }
        return count;
    }

    template<class OrderTag>
    double ScoreMatrix<OrderTag>::percentFilled() const {
        double size = (double) getSize();
        int filled = 0;
        for (int j = 0; j < getNumCols(); j++) {
            for (int i = 0; i < getNumRows(); i++) {
                const ScoreCell* pCell = getCell(i, j);
                if (pCell->score_ > -INF) {
                    filled++;
                }
            }
        }
        return ((double) filled)/size*100.0;
    }

    template<class OrderTag>
    double ScoreMatrix<OrderTag>::getMaxScore() const {
        double maxScore = -INF;
        for(size_t i = 0; i < numCols_; i++)
        {
            const ScoreCell * pCell = getCell(numRows_-1, i);
            if (pCell->score_ > maxScore) 
                maxScore = pCell->score_;
        }
        //cout << "Max Score: " << maxScore << "\n";
        return maxScore;
    }

    template<class OrderTag>
    double ScoreMatrix<OrderTag>::getMaxScoreByRow(size_t row) const {
        double maxScore = -INF;
        if (row >= numRows_) {
            return maxScore;
        }
        for (size_t i = 0; i < numCols_; i++){
            const ScoreCell * pCell = getCell(row-1, i);
                maxScore = pCell->score_;
        }
        return maxScore;
    }

    template<class OrderTag>
    size_t ScoreMatrix<OrderTag>::countColor(ScoreCellColor col) const {
        size_t count = 0;
        for (int j = 0; j < getNumCols(); j++) {
            for (int i = 0; i < getNumRows(); i++) {
                const ScoreCell* pCell = getCell(i, j);
                if (pCell->color_ == col) {
                    count++;
                }
            }
        }
        return count;
    }

    template<class OrderTag>
    template<class OtherOrderTag>
    bool ScoreMatrix<OrderTag>::operator==(const ScoreMatrix<OtherOrderTag>& sm_o) {

        bool dim_match = (getNumCols() == sm_o.getNumCols()) && (getNumRows() == sm_o.getNumRows());
        if (!dim_match) { std::cout << "Dimensions don't match!\n"; return false; }
        size_t m = getNumRows();
        size_t n = getNumCols();
        for(size_t i = 0; i < m; i++) {
            for(size_t j = 0; j < n; j++) {
                const ScoreCell* p1 = getCell(i, j);
                const ScoreCell* p2 = sm_o.getCell(i, j);
                if ( !(*p1 == *p2) ) { std::cout << "Cells don't match! 1: " << *p1 << ", 2: " << *p2 << "\n"; return false; }
            }
        }
        return true;
    }

    template<class OrderTag>
    std::ostream& operator<<(std::ostream& os, const ScoreMatrix<OrderTag>& mat) {
        const size_t nrow = mat.getNumRows();
        const size_t ncol = mat.getNumCols();
        for (size_t i = 0; i < nrow; i++) {
            for (size_t j = 0; j < ncol; j++) {
                os << *mat.getCell(i, j) << " , ";
            }
            os << "\n";
        }
        return os;
    }

    template<class OrderTag>
    void print_filled_by_row(std::ostream& os, const ScoreMatrix<OrderTag>& sm) {
    
        const size_t num_rows = sm.getNumRows();
        const size_t num_cols = sm.getNumCols();
        for(size_t i = 0; i < num_rows; i++) {
            size_t num_filled = sm.countFilledByRow(i);
            os << i << "\t" << num_filled << "\t" << double(num_filled)/double(num_cols) << "\n";
        }

    }




}

#endif
