#ifndef MALIGNER_COMMON_TYPES
#define MALIGNER_COMMON_TYPES

#include <vector>

typedef std::vector<int> IntVec;  
typedef std::vector<double> DoubleVec;
typedef std::vector<bool> BoolVec;

// typedef std::shared_ptr< IntVec > IntVecPtr;
// typedef std::shared_ptr< BoolVec > BoolVecPtr;
typedef IntVec* IntVecPtr;
typedef BoolVec* BoolVecPtr;

#endif