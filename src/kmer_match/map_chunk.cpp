#include "map_chunk.h"

namespace kmer_match {

  std::ostream& operator<<(std::ostream& os, const MapChunk& mc) {

    os << "(" << mc.start_ << "," << mc.end_ << "," << mc.size_ << ")";
    
    return os;

  }
  
}