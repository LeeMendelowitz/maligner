#include <iostream>

#include "ScoreMatrix.h"
#include "globals.h"
//using Constants::INF;

size_t ScoreMatrix::countFilledCells() {
    size_t count = 0;
    for (size_t i = 0; i < data_.size(); i++)
    {
        if (data_[i].score_ > -INF)
            count++;
        /*
         cout << "Score: " << data_[i].score_ << "\tBackPtrs: " << data_[i].backPointers_.size()
                                            << "\tForwardPtrs: " << data_[i].forwardPointers_.size()
                                            << endl;*/
    }
    return count;
}

double ScoreMatrix::getMaxScore() {
    double maxScore = -INF;
    for(size_t i = 0; i < numCols_; i++)
    {
        ScoreCell * pCell = getCell(numRows_-1, i);
        if (pCell->score_ > maxScore) 
            maxScore = pCell->score_;
    }
    //cout << "Max Score: " << maxScore << "\n";
    return maxScore;
}

double ScoreMatrix::getMaxScoreByRow(size_t row) {
    double maxScore = -INF;
    if (row >= numRows_) {
        return maxScore;
    }
    for (size_t i = 0; i < numCols_; i++){
        ScoreCell * pCell = getCell(row-1, i);
            maxScore = pCell->score_;
    }
    return maxScore;
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
