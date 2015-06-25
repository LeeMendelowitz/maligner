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

std::ostream& maligner_dp::operator<<(std::ostream& os, const ScoreCellFullOutput& cell)
{
    const ScoreCell* p_cell = cell.get_cell();

    if (p_cell == nullptr) {
      os << "[ nullptr ]";
      return os;
    }

    const ScoreCell& c = *p_cell;

    os << "[("
       << c.q_ << " "
       << c.r_ << " "
       << c.qm_ << " "
       << c.rm_ << " "
       << c.ref_start_ << " "
       << c.score_ << " "
       << c.m_score_ << " "
       << ") tgt: " 
       << ScoreCellFullOutput(c.backPointer_);

    return os;

}


// bool maligner_dp::ScoreCellPointerCmp(const ScoreCell* p1, const ScoreCell* p2) {
//   return (p1->score_ > p2->score_);
// }
