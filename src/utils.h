#ifndef UTILS_H
#define UTILS_H

#include <vector>

// Sum the values in the vector v from start index to end (exclusive)
template< class T>
T sum(const std::vector<T>& v, int start, int end) {

  T s = 0;

  assert(start <= end);
  if (start >= v.size()) {
    return 0;
  }
  
  for (int i = start; i < end; i++) {
    s += v[i];
  }

  return s;
}

template< class T>
T sum(const std::vector<T>& v) {

  T s = 0;

  const int start = 0;
  const int end = v.size();
  for (int i = start; i < end; i++) {
    s += v[i];
  }

  return s;
}



#endif