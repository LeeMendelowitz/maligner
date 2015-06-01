#include "partialsums.h"

namespace maligner_dp {
  std::ostream& operator<<(std::ostream& os, const maligner_dp::PartialSums& ps) {
    const size_t m = ps.m_;
    const size_t n = (ps.d_.size()/m);
    for(size_t i = 0; i < n; i++) {
        os << i << "| ";
        for(size_t j = 0; j < m; j++) {
          os << ps(i, j) << " ";
        }
        os << "\n";
    }
    return os;
  } 
}