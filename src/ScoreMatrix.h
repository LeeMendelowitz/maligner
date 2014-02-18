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

#include "ScoreCell.h"


class ScoreMatrix
{
public:

    ScoreMatrix(size_t numRows = 0, size_t numCols = 0);
    void resize(size_t numRows, size_t numCols);
    size_t getSize() const { return size_; }
    size_t getCapacity() const { return data_.size(); }
    size_t getNumRows() const { return numRows_; }
    size_t getNumCols() const { return numCols_; }
    ScoreCell * getCell(size_t row, size_t col);
    ScoreCell * getCell(const IntPair& coord);
    ScoreCell * getCell(size_t n);
    const ScoreCell * getCell(size_t row, size_t col) const;
    const ScoreCell * getCell(const IntPair& coord) const;
    const ScoreCell * getCell(size_t n) const;
    void resetCells();

    // Summary Functions
    size_t countFilledCells() const;
    double percentFilled() const;
    double getMaxScore() const;
    double getMaxScoreByRow(size_t row) const;

private:
    size_t numRows_;
    size_t numCols_;
    size_t size_;
    ScoreCellVec data_;
    
};


std::ostream& operator<<(std::ostream& os, ScoreMatrix& mat);



////////////////////////////////////////////////////////////
inline ScoreMatrix::ScoreMatrix(size_t numRows, size_t numCols) :
    numRows_(numRows),
    numCols_(numCols),
    size_(numRows_*numCols_)
{
    resize(numRows, numCols);
}


inline void ScoreMatrix::resetCells() {

    const size_t N = data_.size();

    // Wipe out all cells
    for (size_t i = 0; i < N; i++) {
        getCell(i)->reset();
    }

    // Set the coordinates for those cells that are in bounds
    for (size_t i = 0; i < numRows_; i++) {
        for (size_t j = 0; j < numCols_; j++) {
            ScoreCell * pCell = getCell(i, j);
            pCell->q_ = i;
            pCell->r_ = j;
        }
    }

}

inline void ScoreMatrix::resize(size_t numRows, size_t numCols) {
    numRows_ = numRows;
    numCols_ = numCols;
    size_ = numRows*numCols;
    if (size_ > data_.size())
    {
        data_.resize(size_);
    }

    resetCells();
}

inline ScoreCell * ScoreMatrix::getCell(size_t row, size_t col) {
    // Column-major order
    return &data_[col*numRows_ + row];
}

inline ScoreCell * ScoreMatrix::getCell(const IntPair& coords)
{
    // Check that coords within bounds?
    return getCell(coords.first, coords.second);
}

inline ScoreCell * ScoreMatrix::getCell(size_t n) {
    // Check that coord within bounds?
    return &data_[n];
}

inline const ScoreCell * ScoreMatrix::getCell(size_t row, size_t col) const {
    // Check that coords within bounds?
    return &data_[col*numRows_ + row];
}

inline const ScoreCell * ScoreMatrix::getCell(const IntPair& coords) const
{
    // Check that coords within bounds?
    return getCell(coords.first, coords.second);
}

inline const ScoreCell * ScoreMatrix::getCell(size_t n) const {
    // Check that coord within bounds?
    return &data_[n];
}



#endif
