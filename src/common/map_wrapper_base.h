#ifndef MAP_WRAPPER_BASE
#define MAP_WRAPPER_BASE

#include <vector>
#include <numeric>
#include <iostream>
#include <string>

#include "map.h"
#include "map_data.h"

namespace maligner_maps {

  class MapWrapper {

  public:

    MapWrapper() = default;

    MapWrapper(const Map& map, bool is_circular, bool is_bounded) :
      map_(map),
      map_data_(map_.name_, map_.frags_.size(), map_.size_, is_circular, is_bounded)
    {
      __compute_index(is_circular);

    }

    MapWrapper(Map&& map, bool is_circular, bool is_bounded) :
      map_(map),
      map_data_(map_.name_, map_.frags_.size(), map_.size_, is_circular, is_bounded)
    {

      __compute_index(is_circular);

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
    IntVec ix_to_locs_;



    ////////////////////////////////////////
    private:

    void __compute_index(bool is_circular) {

      FragVec& frags = map_.frags_;

      // If the map is circular, double the fragments (for the purposes of alignment)
      if(is_circular) {
        frags.reserve(2*frags.size());
        frags.insert(frags.end(), frags.begin(), frags.end()); 
      }

      // Compute the index to start loc
      ix_to_locs_.resize(frags.size());
      int s = 0;
      for(size_t i = 0; i < frags.size(); i++) {
        ix_to_locs_[i] = s;
        s += frags[i];
      }
    } 
    

  };

  typedef std::vector<MapWrapper> MapWrapperVec;

}

#endif