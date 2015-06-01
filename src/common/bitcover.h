#ifndef BITCOVER_H
#define BITCOVER_H

#include <vector>
#include <algorithm>

///////////////////////////////////////////////////////////////////////////
// A BitVector for marking locations that have been covered in a reference.
class BitCover {

public:
  
  BitCover() : cover_(0, false) {};
  BitCover(size_t n) : cover_(n, false) {};

  // Mark locations from begin to end as covered.
  void cover(size_t begin, size_t end) {
    const size_t e = std::max(std::min(end, cover_.size()), begin);
    for(size_t i = begin; i < e; i++) {
      cover_[i] = true; 
    }
  }

  // Return true if any part of the interval is covered
  bool is_covered(size_t begin, size_t end) {
    for(size_t i = begin; i < end; i++) {
      if (cover_[i]) return true;
    }
    return false;
  }

  bool is_covered(size_t ind) {
    return cover_[ind];
  }

  bool operator[](size_t ind) {
    return cover_[ind];
  }

  // A safer version of cover that checks bounds
  void cover_safe(size_t start, size_t end) {

    size_t e = (start < end) ? end : start;

    if (e > cover_.size()) {
      cover_.resize(e, false);
    }

    for(size_t i = start; i < e; i++) {
      cover_[i] = true;
    }

  }

  void reset() {
    const size_t N = cover_.size();
    reset(N);
  }

  void reset(size_t n) {
    cover_.clear();
    cover_.resize(n, false);
  }


private:
  std::vector<bool> cover_;

  };

#endif