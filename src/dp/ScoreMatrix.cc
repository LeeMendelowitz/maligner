    #include <iostream>

#include "ScoreMatrix.h"
#include "globals.h"

namespace maligner_dp {

    using Constants::INF;
    using namespace maligner_dp;

    size_t ScoreMatrix::countFilledCells() const {
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

    size_t ScoreMatrix::countFilledByRow(size_t row) const {
        
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

    double ScoreMatrix::percentFilled() const {
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

    double ScoreMatrix::getMaxScore() const {
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

    double ScoreMatrix::getMaxScoreByRow(size_t row) const {
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

    size_t ScoreMatrix::countColor(ScoreCellColor col) const {
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

    std::ostream& operator<<(std::ostream& os, ScoreMatrix& mat) {
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

   void print_filled_by_row(std::ostream& os, const ScoreMatrix& sm) {
    
        const size_t num_rows = sm.getNumRows();
        const size_t num_cols = sm.getNumCols();
        for(size_t i = 0; i < num_rows; i++) {
            size_t num_filled = sm.countFilledByRow(i);
            os << i << "\t" << num_filled << "\t" << double(num_filled)/double(num_cols) << "\n";
        }

   }

}
