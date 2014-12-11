#ifndef TYPES_H
#define TYPES_H

#include <vector>
#include <boost/shared_ptr.hpp>

typedef std::vector<int> IntVec;
typedef std::vector<bool> BoolVec;

typedef boost::shared_ptr< IntVec > IntVecPtr;
typedef boost::shared_ptr< BoolVec > BoolVecPtr;


#endif