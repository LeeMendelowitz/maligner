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
#include <cmath>

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
            qm_(0),
            rm_(0),
            ref_start_(-1),
            backPointer_(nullptr),
            color_(ScoreCellColor::WHITE)
        { };

        ScoreCell(int q, int r) :
            q_(q),
            r_(r),
            qm_(0),
            rm_(0),
            ref_start_(-1),
            backPointer_(nullptr),
            color_(ScoreCellColor::WHITE)
        { };

        ScoreCell(IntPair ip) :
            q_(ip.first),
            r_(ip.second),
            qm_(0),
            rm_(0),
            ref_start_(-1),
            backPointer_(nullptr),
            color_(ScoreCellColor::WHITE)
        { };


        IntPair key() const { 
            return IntPair(q_, r_);
        }

        bool operator==(const ScoreCell& o) const {

            const ScoreCell* b1 = backPointer_;
            const ScoreCell* b2 = o.backPointer_;
            const double TOL = 1E-9;

            bool coord_score_match = (q_ == o.q_) && 
                                     (r_ == o.r_) &&
                                     (qm_ == o.qm_) &&
                                     (rm_ == o.rm_) &&
                                     (std::abs(score_ - o.score_) < TOL) &&
                                     (std::abs(m_score_ - o.m_score_) < TOL) &&
                                     (b1 == nullptr) == (b2 == nullptr);
            
            if (!coord_score_match) return false;

            // Check that the cells pointed to match in coordinates.
            if (b1 != nullptr) {
                if (b1->q_ != b2->q_ || b1->r_ != b2->r_) return false;
            }

            return true;

        }  

        bool is_valid() const {
            return score_ > -INF;
        }

        void reset();

        void setColor(ScoreCellColor color) {
            color_ = color;
        }

        int q_;
        int r_;
        int qm_; // cumulative query misses
        int rm_; // cumulative reference misses
        int ref_start_; // starting location in the reference of the trail that passes through.
        double score_;
        double m_score_;
        ScoreCell * backPointer_; // back pointer for DP solution path
        ScoreCellColor color_;

        friend std::ostream& operator<<(std::ostream&, const ScoreCell& cell);
    };


    inline void ScoreCell::reset()
    {
        q_ = -1;
        r_ = -1;
        qm_ = 0;
        rm_ = 0;
        ref_start_ = -1;
        backPointer_ = nullptr;
        score_ = -INF;
        m_score_ = INF;
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

    struct ScoreCellPointerCmp {
        bool operator()(const ScoreCell* p1, const ScoreCell* p2) {
            return p1->score_ > p2->score_;
        }
    };


    typedef std::vector<ScoreCell> ScoreCellVec;
    std::ostream& operator<<(std::ostream&, const ScoreCell& cell);

}

#endif
