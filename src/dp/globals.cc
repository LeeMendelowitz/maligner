#include "globals.h"
#include <math.h>
#include <limits>

namespace maligner_dp {

    namespace Constants
    {
        const double INF = std::numeric_limits<double>::max();
        const double MAX_LOST_FRAC = 0.20;
        const int FRAG_CUTOFF = 0; // bp

        // Note: this value for SIGMA is taken from Valouev et al (10.1089/cmb.2006.13.442)
        // They use SIGMA^2 = 0.300 kbp for fragments measured in kbp.
        // Error Model: Optical Frag Length ~ N(C, SIGMA^2*C) where C is contig frag length
        const double SIGMA2 = 100; // Units: bp. 
        const double SIGMA = sqrt(SIGMA2); // Units: (bp)^(0.5)
    }


    // Options set by ParseArgs
    namespace opt
    {
        string silicoMap;
        vector<string> opticalMapList;
        bool   circular = false;
    //    bool   allowFalseCuts = false;
    //  bool   matchFragmentOnce = false;
        bool   noReverse = false;
        bool   oneToOneMatch = false;
        bool   useBoundaries = true;
        bool   allowGaps = false;
        double pThreshold = 0.05;
        string outputPrefix = "output";
        double sdMin = 4;
        double sdMax = 12;
        //int smallFragCutoff = 800; // 800 bp
        double C_r_contig = 3.0; // Cost for missed contig site
        double C_r_optical = 5.0; // Cost for missed optical site
    //    double falseCutRate = 1E-6; // 1 false cut per 1E6 bp (1 Mbp)

        // Match/Filter params
        double maxMissRateContig = 0.50; // percent of missed sites allowed
        double minLengthRatio = 0.90; // minimum length ratio between aligned inner contig and optical map fragments
        double avgChi2Threshold = 10.0; // Max. avg. chi2 score across aligned blocks for an alignment to be considered
        int maxMatchesPerContig = 5;
        int minContigHits = 5; 

    //   int maxGapSize = 0; // Maximum open gap allowed for missed fragments
        int numPermutationTrials = 0;
        int numThreads = 1;
        int maxChunkMissesQuery = 0;
        int maxChunkMissesReference = 3;
        int smallFrag = 2000; // Cutoff of small fragments (bp)
        double smallFragSlope = C_r_contig/smallFrag;
        bool localAlignment = false;

        ///////////////////////////////////////////////////////
        // Parameters for local alignment

        // Local sizing scoring function shape parameters
        double H = 5.0; // Height of parabola. Positive parameter specifying reward for zero sizing error.
        double T = 2.0; // Root of parabola. Specifies the number of standard deviations of sizing error which gets zero cost.
    }
    
}
