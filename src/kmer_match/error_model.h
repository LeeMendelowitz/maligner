#ifndef ERROR_MODEL_H
#define ERROR_MODEL_H

#include <vector>
#include <utility>

typedef std::vector<int> IntVec;
typedef std::pair<int, int> IntPair;
typedef std::vector<IntPair> IntPairVec;


class ErrorModel {
public:
  
  ErrorModel(double rel_error, double min_error) :
    rel_error_(rel_error), min_error_(min_error) {};


  IntPairVec compute_bounds(const IntVec& frags);
  IntPairVec compute_bounds(IntVec::const_iterator s, IntVec::const_iterator e);

  double rel_error_;
  double min_error_;

};



#endif