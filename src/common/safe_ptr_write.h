#ifndef LMM_SAFE_WRITE_H
#define LMM_SAFE_WRITE_H

#include <ostream>
#include <sstream>

namespace lmm_utils {

  // template<class T>
  // struct SAFE_WRITE {

  //   SAFE_WRITE(const T* p) : p_{p} {};

  //   T * p_;
  // };

  // template<class T>
  // std::ostream& operator<<(std::ostream& os, const SAFE_WRITE<T>& w) {

  //   if (w == nullptr) {
  //     os << "nullptr";
  //     return os;
  //   }

  //   os << *w;
  //   return os;
  // }

  template<class T>
  std::string SAFE_WRITE(const T* w) {

    if (w == nullptr) {

      return std::string{"nullptr"};

    } else {

      std::ostringstream oss;
      oss << *w;
      return oss.str();

    }

  }

}

#endif
