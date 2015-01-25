#ifndef MAP_FRAG_H
#define MAP_FRAG_H

#include <string>
#include <vector>
#include <iostream>
#include <iterator>
#include "map.h"

namespace kmer_match {
  
  using namespace maligner_maps;

  enum MapFragStatus {VALID, INVALID};

  /////////////////////////////////////////////////////////////////
  /*
  MapFrag is an iterator like object which refers to a fragment
  from a map. We can increment it to determine the next fragment or
  decrement it to determine reverse map.
  */
  class MapFragForward {

  public:

    MapFragForward() = delete; // No default constructor

    MapFragForward(const Map* pMap, size_t index) :
      pMap_(pMap),
      index_(index),
      valid_(true) { };

    MapFragForward(const MapFragForward& o) = default;
    MapFragForward(MapFragForward&& o) = default;
    MapFragForward& operator=(const MapFragForward& o) = default;
    MapFragForward& operator=(MapFragForward&& o) = default;

    MapFragForward& operator++() {

      if (!valid_) { return *this; }

      index_++;
      valid_ = (index_ < pMap_->frags_.size()) && (index_ >= 0);

      return *this;

    }

    MapFragForward& operator--() {

      if (!valid_) { return *this; }

      index_--;
      valid_ = (index_ < pMap_->frags_.size()) && (index_ >= 0);

      return *this;

    }

    MapFragForward& operator+=(int o) {

      if (!valid_) { return *this; }

      index_ += o;
      valid_ = (index_ < pMap_->frags_.size()) && (index_ >= 0);

      return *this;

    }

    bool operator<(const MapFragForward& o) const {
      return this->get_size() < o.get_size();
    }

    bool operator<(int o) const {
      return this->get_size() < o;
    }

    std::ostream& write_context(std::ostream& os, size_t N) const;

    int get_size() const{
      if(!valid_) return -1;
      return pMap_->frags_[index_];
    }





    ////////////////////////////////////////////
    const Map* pMap_;

    // These can change during search
    size_t index_;
    bool valid_;

  };

  ///////////////////////////////////////////////////////////////////
  class MapFragReverse {
  public:

    MapFragReverse() = delete; // No default constructor

    MapFragReverse(const Map* pMap, size_t index) :
      pMap_(pMap),
      index_(index),
      valid_(true) { };

    MapFragReverse(const MapFragReverse& o) = default;
    MapFragReverse(MapFragReverse&& o) = default;
    MapFragReverse& operator=(const MapFragReverse& o) = default;
    MapFragReverse& operator=(MapFragReverse&& o) = default;

    MapFragReverse& operator++() {

      if (!valid_) {
        return *this;
      }

      index_--;

      valid_ = (index_ < pMap_->frags_.size()) && (index_ >= 0);


      return *this;

    }

    MapFragReverse& operator--() {

      if(!valid_) { return *this; }

      index_++;
      valid_ = (index_ < pMap_->frags_.size()) && (index_ >= 0);

      return *this;

    }

    MapFragReverse& operator+=(int o) {

      if (!valid_) { return *this; }

      index_ -= o; // Notice -=, because forward is in the reverse directino for MapFragReverse
      valid_ = (index_ < pMap_->frags_.size()) && (index_ >= 0);

      return *this;

    }

    bool operator<(const MapFragReverse& o) const {
      return this->get_size() < o.get_size();
    }

    bool operator<(int o) const {
      return this->get_size() < o;
    }

    std::ostream& write_context(std::ostream& os, size_t N) const;

    int get_size() const{
      return pMap_->frags_[index_];
    }

    ////////////////////////////////////////////
    const Map* pMap_;

    // These can change during search
    size_t index_;
    bool valid_;

  };

  template< class T >
  class FragIsInvalid {
      public:
      bool operator()(const T& frag) const {
          return !frag.valid_;
      }
  };

  template< class T>
  class FragIsValid {
      public:
      bool operator()(const T& frag) const {
          return frag.valid_;
      }
  };

  bool operator<(int val, const MapFragForward& frag);
  bool operator<(int val, const MapFragReverse& frag);

  typedef std::vector<MapFragForward> MapFragForwardVec;
  typedef std::vector<MapFragReverse> MapFragReverseVec;

  ////////////////////////////////////////////////////////
  // Get fragments from a map. This is templated so that
  // we can return either forward fragments or reverse fragments
  // depending on the context.
  template< class T>
  std::vector<T> get_frags(const Map& map) {
    
    std::vector<T> ret;

    size_t index;
    size_t L = map.frags_.size();

    for(index = 0; index < L; index++) {
      ret.emplace_back(&map, index);
    }

    return ret;
  }

  //////////////////////////////////////////////////////////////
  // Get fragments from a vector of maps.
  template< class T>
  std::vector<T> get_frags(const MapVec& maps) {
    
    std::vector<T> ret;
    auto ret_back_inserter = std::back_inserter(ret);
    
    size_t num_frags = 0;
    for(MapVec::const_iterator mi = maps.begin(); mi != maps.end();  mi++){
      num_frags += mi->frags_.size();
    }
    ret.reserve(num_frags);

    for (MapVec::const_iterator mi = maps.begin(); mi != maps.end(); mi++) {
      
      std::vector<T> map_frags = get_frags<T>(*mi);

      // Move all fragments into the database
      std::move(map_frags.begin(), map_frags.end(), ret_back_inserter);
    }

    return ret;
  }

  std::ostream& operator<<(std::ostream& os, const MapFragForward& m);
  std::ostream& operator<<(std::ostream& os, const MapFragReverse& m);

  template< class T>
  std::ostream& write_context(std::ostream& os, const Map& map, const std::vector<T>& frags, size_t N) {
    for(auto i = frags.begin(); i != frags.end(); i++) {
      os << map.name_ << " ";
      i->write_context(os, N);
      os << "\n";
    }
    return os;
  } 

  template< class T>
  std::ostream& operator<<(std::ostream& os, const std::vector<T> vec) {
    for (auto i = vec.begin(); i != vec.end(); i++) {
      os << *i << " ";
    }
    return os;
  }

}
#endif
