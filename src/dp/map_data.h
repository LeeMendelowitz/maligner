#ifndef MAPDATA_H
#define MAPDATA_H

#include <iostream>

namespace maligner_dp {

  // Store MetaData for a map.
  class MapData {
  public:

      MapData() = default;
      
      MapData(const std::string& map_name,
        size_t num_frags,
        bool is_bounded = false) :
        map_name_(map_name),
        num_frags_(num_frags),
        is_bounded_(is_bounded)
      {
        num_constructs++;
      };

      MapData(const MapData& o) :
        num_frags_(o.num_frags_),
        map_name_(o.map_name_),
        is_bounded_(o.is_bounded_)
      {
        num_copies++;
      };

      static void print_debug_stats() {
        std::cerr << "MapData Stats:\n"
                 << "\tnum_constructs: " << num_constructs << "\n"
                 << "\tnum_copies: " << num_copies << "\n";
      }

      std::string map_name_;
      size_t num_frags_; // Total fragments in the map. This may be larger than the number of
                         // fragments provided for alignment in the case that we are aligning a slice
                         // of the map.
      bool is_bounded_; // Are the leftmost/rightmost fragments bounded by sites? NOTE: This is not yet used in code, and duplicates an argument in AlignOpts

      static int num_constructs;
      static int num_copies;
  };

}


#endif