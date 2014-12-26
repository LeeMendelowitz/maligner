#include <vector>
#include "chunk.h"

namespace maligner_dp {

  std::ostream& operator<<(std::ostream& os, const Chunk& chunk) {

    os << "([" << chunk.start << ", " << chunk.end << "], " 
       << chunk.size << ")";

    return os;

  }
  
}
