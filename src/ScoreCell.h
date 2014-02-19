#ifndef SCORECELL_H
#define SCORECELL_H

// ScoreCell is a cell in the dynamic programming table and
// identifies its successors/predecessors.
// It acts as a node in the DAG underlying the dynamic programming
// table.

#include<vector>
#include<set>
#include <unordered_map>
#include <unordered_set>
#include <ostream>

#include "globals.h"
using Constants::INF;

class ScoreCell;
typedef std::vector<ScoreCell *> ScoreCellPVec;
typedef std::unordered_set<ScoreCell *> ScoreCellSet;
typedef std::pair<int, int> IntPair;

class ScoreCell
{
    public:

    ScoreCell() :
        q_(-1),
        r_(-1),
        backPointer_(nullptr)
    { };

    ScoreCell(int q, int r, bool inPlay = false) :
        q_(q),
        r_(r),
        backPointer_(nullptr)
    { };

    ScoreCell(IntPair ip, bool inPlay = false) :
        q_(ip.first),
        r_(ip.second),
        backPointer_(nullptr)
    { };

    IntPair key() const { 
        return IntPair(q_, r_);
    }

    bool operator==(const ScoreCell& other) const {
        return (key() == other.key());
    } 

    void reset();

    int q_;
    int r_;
    double score_;
    ScoreCell * backPointer_; // back pointer for DP solution path
};


inline void ScoreCell::reset()
{
    q_ = -1;
    r_ = -1;
    backPointer_ = nullptr;
    score_ = -INF;
}


/*
namespace std {
    template <>
    struct std::hash< ScoreCell > {
        public:
        size_t operator()(const ScoreCell& x) const throw() {
            std::hash<IntPair> hasher;
            return hasher(x.key());
        }
    };
}
*/

typedef std::vector<ScoreCell> ScoreCellVec;
std::ostream& operator<<(std::ostream&, const ScoreCell& cell);

#endif
