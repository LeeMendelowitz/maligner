#ifndef MAP_WRAPPERS_H
#define MAP_WRAPPERS_H

#include "map.h"
#include "map_data.h"
#include "partialsums.h"  

namespace maligner_dp {

  using namespace maligner_maps;

  ////////////////////////////////////////////////////////////
  // These wrappers around a maligner_maps::Map
  // compute the MapData (i.e. metadata)
  // and potentially PartialSums and SDInv2, 
  // depending on which wrapper is used.
  // This is a convenience to bundle all map data together

  struct MapWrapper {

    MapWrapper(const Map& m, int num_missed_sites) :
      m_(m),
      md_(m_.name_, m_.frags_.size())
    {

    }

    Map m_;
    MapData md_;

    size_t num_frags() const {
      return m_.frags_.size();
    }


  };

  struct QueryMapWrapper : public MapWrapper {

    QueryMapWrapper(const Map& m, int num_missed_sites) :
        MapWrapper(m, num_missed_sites),
        ps_forward_(m_.frags_, num_missed_sites, PartialSums::forward_tag()),
        ps_reverse_(m_.frags_, num_missed_sites, PartialSums::reverse_tag())
      {

      }

      PartialSums ps_forward_;
      PartialSums ps_reverse_;

  };

  struct RefMapWrapper : public MapWrapper {

    RefMapWrapper(const Map& m, int num_missed_sites, double sd_rate, double min_sd) :
      MapWrapper(m, num_missed_sites),
      ps_(m_.frags_, num_missed_sites),
      sd_inv_2_(ps_, sd_rate, min_sd )
    {

    }

    PartialSums ps_;
    SDInv2 sd_inv_2_;

  };

}

#endif