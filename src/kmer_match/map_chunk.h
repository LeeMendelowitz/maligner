#ifndef MAP_CHUNK
#define MAP_CHUNK

#include <vector>
#include <numeric>
#include <iostream>

#include "map.h"

namespace kmer_match {
  
  class MapChunk {
  public:

    MapChunk(const Map * pMap, size_t s, size_t e) : 
      pMap_(pMap), start_(s), end_(e),
      size_(std::accumulate(pMap->frags_.begin() + s, pMap->frags_.begin() + e, 0))
    {};

    const Map * pMap_;
    size_t start_;
    size_t end_;
    int size_;

    size_t num_frags() const {
      return end_ - start_;
    }

    size_t num_unmatched() const {
      return end_ - start_ - 1;
    }

  };


  typedef std::vector<MapChunk> MapChunkVec;
  typedef std::vector<MapChunk*> MapChunkPVec;
  typedef std::vector<const MapChunk*> ConstMapChunkPVec;

  std::ostream& operator<<(std::ostream& os, const MapChunk& mc);


}
  #endif