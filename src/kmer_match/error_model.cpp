#include "error_model.h"
#include <algorithm>
#include <utility>

using namespace std;

IntPairVec ErrorModel::compute_bounds(const IntVec& frags) {

  IntPairVec ret(frags.size());

  IntPairVec::iterator ret_i = ret.begin();
  for(IntVec::const_iterator i = frags.begin(); i < frags.end(); i++, ret_i++) {
    const int size = *i;
    const double delta = std::max(rel_error_ * double(size), min_error_);
    int lb = int(std::max(0.0, size - delta)); // Floor
    int ub = int(size + delta) + 1; // Ceiling
    *ret_i = std::move(IntPair(lb, ub));
  }

  return ret;

}

IntPairVec ErrorModel::compute_bounds(const IntVec::const_iterator s, const IntVec::const_iterator e) {

  IntPairVec ret(e - s);
  IntPairVec::iterator ret_i = ret.begin();
  for(IntVec::const_iterator i = s; i < e; i++, ret_i++) {
    const int size = *i;
    const double delta = std::max(rel_error_ * double(size), min_error_);
    int lb = int(std::max(0.0, size - delta)); // Floor
    int ub = int(size + delta) + 1; // Ceiling
    *ret_i = std::move(IntPair(lb, ub));
  }

  return ret;

}