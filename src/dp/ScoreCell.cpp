#include "ScoreCell.h"
#include <algorithm>

using namespace maligner_dp;

std::ostream& maligner_dp::operator<<(std::ostream& os, const ScoreCell& cell)
{
    os << "[(" << cell.q_ << "," << cell.r_ << "), " << cell.score_ << ",";
    const ScoreCell * pTgt = cell.backPointer_;
    if (pTgt)
        os << " tgt: (" << pTgt->q_ << "," << pTgt->r_ << ") ]";
    else
        os << " tgt: nullptr ]";
    return os;
}

bool maligner_dp::ScoreCellPointerCmp(const ScoreCell* p1, const ScoreCell* p2) {
  return (p1->score_ > p2->score_);
}
