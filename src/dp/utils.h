#ifndef UTILS_H
#define UTILS_H

#include <vector>

namespace maligner_dp {

  // Sum the values in the vector v from start index to end (exclusive)
  template< class T>
  T sum(const std::vector<T>& v, int start, int end) {

    T s = 0;

    assert(start <= end);
    if (start >= (int) v.size()) {
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

  template< typename T>
  std::ostream& operator<<(std::ostream& os, std::vector<T> vec) {
    size_t n = vec.size();
    os << "[ ";
    for (size_t i = 0; i < n; i++) {
      os << vec[i];
      if (i != n-1) {
        os << ", ";
      }
    }
    os << " ]";
    return os;
  }

  template< typename T>
  std::ostream& operator<<(std::ostream& os, std::vector<T*> vec) {
    size_t n = vec.size();
    os << "[ ";
    for (size_t i = 0; i < n; i++) {
      os << *vec[i];
      if (i != n-1) {
        os << ", ";
      }
    }
    os << " ]";
    return os;
  }

}
#endif