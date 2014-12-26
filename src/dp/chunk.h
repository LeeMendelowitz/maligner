#ifndef CHUNK_H
#define CHUNK_H

#include <vector>
#include <ostream>

namespace maligner_dp {

  class Chunk {
    public:
      
      Chunk() : start(0), end(0), size(0), is_boundary(false) {};

      Chunk(int s, int e, int sz, bool isb) :
        start(s), end(e), size(sz), is_boundary(isb) {};
      int start; // start coordinate, inclusive
      int end; // end coordinate, exclusive
      int size; // total size of chunk (bp)
      bool is_boundary; // is boundary chunk

      bool operator==(const Chunk& other) const {
        return (start == other.start &&
                end == other.end &&
                size == other.size &&
                is_boundary == other.is_boundary);
      }

      int num_misses() const {
        return end - start - 1;
      }

  };

  typedef std::vector<Chunk> ChunkVec;

  std::ostream& operator<<(std::ostream& os, const Chunk& chunk);
  
}

#endif

