#ifndef MAP_WRAPPERS_H
#define MAP_WRAPPERS_H

#include "map.h"
#include "map_data.h"
#include "map_wrapper_base.h"
#include "partialsums.h"  

namespace maligner_dp {

  using namespace maligner_maps;

  struct QueryMapWrapper : public maligner_maps::MapWrapper {

    QueryMapWrapper(const Map& map, int num_missed_sites) :
        MapWrapper(map, false, false),
        ps_forward_(map_.frags_, num_missed_sites, PartialSums::forward_tag()),
        ps_reverse_(map_.frags_, num_missed_sites, PartialSums::reverse_tag())
      {

      }

      PartialSums ps_forward_;
      PartialSums ps_reverse_;

  };

  struct RefMapWrapper : public maligner_maps::MapWrapper {

    RefMapWrapper(const Map& map, bool is_circular, int num_missed_sites, double sd_rate, double min_sd) :
      MapWrapper(map, is_circular, false),
      ps_(map_.frags_, num_missed_sites),
      sd_inv_(ps_, sd_rate, min_sd ),
      frags_(map.frags_)
    {

      if(is_circular) {

        // Circularize the frags_ by doubling.
        // Copy fragments from the Map's fragment vector.
        frags_.insert(frags_.end(), map.frags_.begin(), map.frags_.end());

      }

    }

    PartialSums ps_;
    SDInv sd_inv_;
    FragVec frags_;

  };

}

#endif
