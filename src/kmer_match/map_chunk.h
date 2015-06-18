#ifndef MAP_CHUNK
#define MAP_CHUNK

#include <vector>
#include <numeric>
#include <iostream>

#include "map_wrapper_base.h"

using namespace maligner_maps;

namespace kmer_match {
  
  class MapChunk {
  public:

    MapChunk(const MapWrapper * pMap, size_t s, size_t e) : 
      pMap_(pMap), start_(s), end_(e),
      size_(std::accumulate(pMap->map_.frags_.begin() + s, pMap->map_.frags_.begin() + e, 0))
    {};

    const MapWrapper * pMap_;
    size_t start_;
    size_t end_;
    int size_;

    const MapWrapper * get_map_wrapper() const {
      return pMap_;
    }

    const Map * get_map() const {
      return &pMap_->map_;
    }

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