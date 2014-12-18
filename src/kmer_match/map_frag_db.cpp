#include "map_frag_db.h"
#include "map_frag.h"
#include <algorithm>
#include <iostream>

using namespace std;


MapFragDB::MapFragDB(const MapVec& maps) :
frags_forward_(get_frags<MapFragForward>(maps)),
frags_reverse_(get_frags<MapFragReverse>(maps)) {
  sort_frags();
}

void MapFragDB::sort_frags() {
  std::sort(frags_forward_.begin(), frags_forward_.end());
  std::sort(frags_reverse_.begin(), frags_reverse_.end());
}


MapFragForwardVec  MapFragDB::query_forward(const IntPairVec& bounds) {
  return query_frags_start_largest<MapFragForward>(bounds, frags_forward_);
}

MapFragReverseVec  MapFragDB::query_reverse(const IntPairVec& bounds) {
  return query_frags_start_largest<MapFragReverse>(bounds, frags_reverse_);
}

void MapFragDB::dump_frags(std::ostream& os) {

  std::ostream_iterator<MapFragForward> out_it_forward(os,"\n");
  std::copy (frags_forward_.begin(), frags_forward_.end(), out_it_forward);
  std::ostream_iterator<MapFragReverse> out_it_rev(os,"\n");
  std::copy (frags_reverse_.begin(), frags_reverse_.end(), out_it_rev);

}
