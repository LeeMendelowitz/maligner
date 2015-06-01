#ifndef MAP_WRAPPER_BASE
#define MAP_WRAPPER_BASE

#include <vector>
#include <numeric>
#include <iostream>
#include <string>

#include "map.h"
#include "map_data.h"
#include "bitcover.h"

namespace maligner_maps {


  //////////////////////////////////////////////////////////////////////////////////////
  // The MapWrapper class brings together structures for a single map.
  //
  // The Map class is meant for parsing maps from a .maps file into a
  // fragment vector.
  //
  // The MapData class maintains meta data about the map (its name, size, is circular,
  // is bounded by restriction sites).
  //
  // The MapWrapper contains these two, plus a structure ix_to_lcos_ for
  // converting an index to a bp location
  //
  // NOTE: If the map is circular (most prokaryotic genomes),
  // then we "double" the vector of fragments to allow use to 
  // align fragments past the origin of replication. We double the fragments in
  // map_.frags_.
  //
  // However, we maintain the number of fragments in the **original** map 
  // in the map_data_.
  ////////////////////////////////////////////////////////////////////////////////////
  class MapWrapper {

  public:

    MapWrapper() = default;

    MapWrapper(const Map& map, bool is_circular, bool is_bounded) :
      map_(map),
      map_data_(map_.name_, map_.frags_.size(), map_.size_, is_circular, is_bounded),
      frags_reverse_(map_.frags_.rbegin(), map_.frags_.rend()),
      bit_cover_(map_.frags_.size())
    {

      // If the map is circular, double the fragments (for the purposes of alignment)
      if(is_circular) {
        __circularize();
      }

      __compute_index();

    }

    MapWrapper(Map&& map, bool is_circular, bool is_bounded) :
      map_(map),
      map_data_(map_.name_, map_.frags_.size(), map_.size_, is_circular, is_bounded),
      frags_reverse_(map_.frags_.rbegin(), map_.frags_.rend()),
      bit_cover_(map_.frags_.size())
    {

      // If the map is circular, double the fragments (for the purposes of alignment)
      if(is_circular) {
        __circularize();
      }

      __compute_index();

    }

    MapWrapper(const Map& map, const MapData& map_data) :
      map_(map),
      map_data_(map_data),
      frags_reverse_(map_.frags_.rbegin(), map_.frags_.rend()) {

        // If the map is circular, double the fragments (for the purposes of alignment)
        if(map_data_.is_circular_) {
          __circularize();
        }

        __compute_index();

    }

    /////////////////////////////////////////////////////////////////////////////
    // Return the non-negative starting location in bp of the fragment at index.
    // -1 signals error
    int index_to_loc(size_t index) {
      if(index >= ix_to_locs_.size()) {
        return -1;
      }
      return ix_to_locs_[index];
    }

    const FragVec& get_frags() const {
      return map_.frags_;
    }

    FragVec get_frags_noncircularized() const {
      const size_t N = num_frags();
      const FragVec& frags = get_frags();
      FragVec ret(N);
      for(size_t i = 0; i < N; i++) {
        ret[i] = frags[i];
      }
      return ret;
    }

    const FragVec& get_frags_reverse() const {
      return frags_reverse_;
    }

    // This returns the number of fragments in the original map,
    // before doubling for circularization.
    size_t num_frags() const {
      return map_data_.num_frags_;
    }

    Map& get_map() {
      return map_;
    }

    std::string get_name() {
      return map_.name_;
    }


    ///////////////////////////////////////
    // Members
    Map map_;
    maligner_maps::MapData map_data_; // Metadata about the map
    FragVec frags_reverse_;
    IntVec ix_to_locs_;
    BitCover bit_cover_;

    ////////////////////////////////////////
    private:


    void __circularize() {

      // Repeat the fragment vector.
      FragVec& frags = map_.frags_;

      if(frags.size() == 2*num_frags()) {
        // Looks like it's already circularized!
        return;
      }

      frags.reserve(2*frags.size());
      frags.insert(frags.end(), frags.begin(), frags.end()); 
      frags_reverse_ = FragVec(frags.rbegin(), frags.rend());
      bit_cover_.reset(frags.size());

    }

    void __compute_index() {

      const FragVec& frags = map_.frags_;
      
      // Compute the index to start loc
      ix_to_locs_.resize(frags.size() + 1);

      int s = 0;

      for(size_t i = 0; i < frags.size(); i++) {
        ix_to_locs_[i] = s;
        s += frags[i];
      }

      ix_to_locs_[frags.size()] = s;

    } 
    

  };

  typedef std::vector<MapWrapper> MapWrapperVec;

}

#endif