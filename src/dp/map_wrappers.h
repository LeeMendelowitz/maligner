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

    QueryMapWrapper(const QueryMapWrapper& o) = default;
    QueryMapWrapper& operator=(const QueryMapWrapper& o) = default;
    
    QueryMapWrapper(QueryMapWrapper&& o) = default;
    QueryMapWrapper& operator=(QueryMapWrapper&& o) = default;


      PartialSums ps_forward_;
      PartialSums ps_reverse_;

    const PartialSums& get_partial_sums_forward() const {
      return ps_forward_;
    }

    const PartialSums& get_partial_sums_reverse() const {
      return ps_reverse_;
    }

  };

  struct RefMapWrapper : public maligner_maps::MapWrapper {

    RefMapWrapper(const Map& map, bool is_circular, int num_missed_sites, double sd_rate, double min_sd) :
      MapWrapper(map, is_circular, false),
      ps_(map_.frags_, num_missed_sites),
      ps_reverse_(get_frags_reverse(), num_missed_sites),
      sd_inv_(ps_, sd_rate, min_sd),
      sd_inv_reverse_(ps_reverse_, sd_rate, min_sd)
    {

    }


    RefMapWrapper(const Map& map, const MapData& map_data,
      int num_missed_sites, double sd_rate, double min_sd) :
      MapWrapper(map, map_data),
      ps_(map_.frags_, num_missed_sites),
      ps_reverse_(get_frags_reverse(), num_missed_sites),
      sd_inv_(ps_, sd_rate, min_sd ),
      sd_inv_reverse_(ps_reverse_, sd_rate, min_sd)
    {

    }

    RefMapWrapper(const RefMapWrapper& o) = default;
    RefMapWrapper& operator=(const RefMapWrapper& o) = default;
    
    RefMapWrapper(RefMapWrapper&& o) = default;
    RefMapWrapper& operator=(RefMapWrapper&& o) = default;

    const PartialSums& get_partial_sums() const {
      return ps_;
    }

    const PartialSums& get_partial_sums_reverse() const {
      return ps_reverse_;
    }

    PartialSums ps_;
    PartialSums ps_reverse_;
    SDInv sd_inv_;
    SDInv sd_inv_reverse_;

  };


}

#endif
