#ifndef SCORECELL_H
#define SCORECELL_H

// ScoreCell is a cell in the dynamic programming table and
// identifies its successors/predecessors.
// It acts as a node in the DAG underlying the dynamic programming
// table.

#include <vector>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <ostream>

#include "globals.h"

namespace maligner_dp {
    
    using maligner_dp::Constants::INF;

    class ScoreCell;
    typedef std::vector<ScoreCell *> ScoreCellPVec;
    typedef std::unordered_set<ScoreCell *> ScoreCellSet;
    typedef std::pair<int, int> IntPair;

    enum class ScoreCellColor { WHITE, BLACK, RED, GREEN};
    class ScoreCell
    {
        public:

        ScoreCell() :
            q_(-1),
            r_(-1),
            backPointer_(nullptr),
            color_(ScoreCellColor::WHITE)
        { };

        ScoreCell(int q, int r) :
            q_(q),
            r_(r),
            backPointer_(nullptr),
            color_(ScoreCellColor::WHITE)
        { };

        ScoreCell(IntPair ip) :
            q_(ip.first),
            r_(ip.second),
            backPointer_(nullptr),
            color_(ScoreCellColor::WHITE)
        { };


        IntPair key() const { 
            return IntPair(q_, r_);
        }

        bool operator==(const ScoreCell& other) const {
            return (key() == other.key());
        } 

        void reset();

        void setColor(ScoreCellColor color) {
            color_ = color;
        }

        int q_;
        int r_;
        double score_;
        ScoreCell * backPointer_; // back pointer for DP solution path
        ScoreCellColor color_;
    };


    inline void ScoreCell::reset()
    {
        q_ = -1;
        r_ = -1;
        backPointer_ = nullptr;
        score_ = -INF;
        color_ = ScoreCellColor::WHITE;
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

    bool ScoreCellPointerCmp(const ScoreCell* p1, const ScoreCell* p2);

    typedef std::vector<ScoreCell> ScoreCellVec;
    std::ostream& operator<<(std::ostream&, const ScoreCell& cell);

}

#endif
