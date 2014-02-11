#include "ScoreCell.h"
#include <algorithm>

std::ostream& operator<<(std::ostream& os, const ScoreCell& cell)
{
    os << "[(" << cell.q_ << "," << cell.r_ << "), " << cell.score_ << ",";
    const ScoreCell * pTgt = cell.backPointer_;
    if (pTgt)
        os << " tgt: (" << pTgt->q_ << "," << pTgt->r_ << ") ]";
    else
        os << " tgt: nullptr ]";
    return os;
}

