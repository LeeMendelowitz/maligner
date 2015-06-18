#ifndef MAPDATA_H
#define MAPDATA_H

#include <iostream>
#include <vector>

namespace maligner_maps {

  typedef std::vector<int> IntVec;

  ///////////////////////////////////////////
  // Store MetaData for a map.
  class MapData {
  public:

      MapData() = default;
      
      MapData(const std::string& map_name,
        size_t num_frags,
        int length,
        bool is_circular = false,
        bool is_bounded = true,
        bool is_random = false) :
        map_name_(map_name),
        num_frags_(num_frags),
        is_circular_(is_circular),
        is_bounded_(is_bounded),
        is_random_(is_random)
      {
        // num_constructs++;

      };

      MapData(const MapData& o) :
        num_frags_(o.num_frags_),
        length_(o.length_),
        map_name_(o.map_name_),
        is_circular_(o.is_circular_),
        is_bounded_(o.is_bounded_),
        is_random_(o.is_random_)
      {
        // num_copies++;
      };

      // static void print_debug_stats() {
      //   std::cerr << "MapData Stats:\n"
      //            << "\tnum_constructs: " << num_constructs << "\n"
      //            << "\tnum_copies: " << num_copies << "\n";
      // }

      std::string map_name_;
      size_t num_frags_; // The number of fragments in the original (non-circularized) map.
      int length_; // number of bp
      bool is_circular_; // Is the map circular?
      bool is_bounded_; // Are the leftmost/rightmost fragments bounded by sites? NOTE: This is not yet used in code, and duplicates an argument in AlignOpts
      bool is_random_; // Is this a random map?
      // static int num_constructs;
      // static int num_copies;
  };

}


#endif
