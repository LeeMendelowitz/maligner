#ifndef GLOBALS_H
#define GLOBALS_H
#include <string>
#include <map>
#include <vector>

using namespace std;

namespace Constants
{
    extern const double INF;
    extern const int MAX_GAP_SIZE;
    extern const double MAX_LOST_FRAC; // maximum fraction of contig that is allowed to be lost
    extern const int FRAG_CUTOFF; // small fragment cutoff
    extern const double SIGMA2; // Optical Frag ~ N(C, SIGMA^2*C) where C is contig frag size. Units: (bp)
    extern const double SIGMA; // sqrt of SIGMA2
}

// Options set by ParseArgs
namespace opt
{
    extern string silicoMap;
    extern vector<string> opticalMapList;
    extern bool   circular;
//    extern bool allowFalseCuts;
//    extern bool   matchFragmentOnce;
    extern bool noReverse;
    extern bool oneToOneMatch;
    extern bool useBoundaries; // if true, treat first and last contig fragment as boundary fragment (bounded by only 1 restriction site, instead of 2).
    extern bool allowGaps;
    extern double pThreshold;
    extern string outputPrefix;
    extern double sdMin;
    extern double sdMax;
    extern double C_r_optical;
    extern double C_r_contig;
    //extern double C_s;

    // Match/Filter options
    extern double maxMissRateContig;
    extern double minLengthRatio;
    extern double avgChi2Threshold;
    extern int maxMatchesPerContig;
    extern int minContigHits;

//    extern double falseCutRate;
//    extern int maxGapSize;
    extern int numPermutationTrials;
    extern int numThreads;
    extern int maxChunkMissesQuery; // maximum number of unaligned sites in interior of alignment block
    extern int maxChunkMissesReference; // maximum number of unaligned sites in interior of alignment block
    extern int smallFrag; // Cutoff of small fragments (bp)
    extern double smallFragSlope;
    extern bool localAlignment;

    ///////////////////////////////////////////////////////
    // Parameters for local alignment

    // Local sizing scoring function shape parameters
    extern double H; // Height of parabola. Positive parameter specifying reward for zero sizing error.
    extern double T; // Root of parabola. Specifies the number of standard deviations of sizing error which gets zero cost.
}


#endif
