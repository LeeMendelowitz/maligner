#include "map_chunk.h"

std::ostream& operator<<(std::ostream& os, const MapChunk& mc) {

  os << "(" << mc.start_ << "," << mc.end_ << "," << mc.size_ << ")";
  
  return os;

}