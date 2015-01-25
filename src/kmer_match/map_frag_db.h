#ifndef MAP_FRAG_DB
#define MAP_FRAG_DB

#include <utility>
#include <vector>
#include <iostream>
#include <algorithm>

#include "map.h"
#include "map_frag.h"

namespace kmer_match {

  using namespace maligner_maps;
  
  typedef std::pair<int, int> IntPair;
  typedef std::vector<IntPair> IntPairVec;

  class MapFragDB {

  public:

    MapFragDB(const MapVec& maps);

    void sort_frags();

    template< class T>
    std::vector<T> query_frags(int lb, int ub, const std::vector<T>& frags);

    template< class T>
    std::vector<T> query_frags(const IntPair& ip, const std::vector<T>& frags);

    template< class T>
    std::vector<T> query_frags(const IntPairVec& bounds, const std::vector<T>& frags);

    template< class T>
    std::vector<T> query_frags_start_largest(const IntPairVec& bounds, const std::vector<T>& frags);


    // Query in the forwards direction
    MapFragForwardVec query_forward(const IntPairVec& bounds);

    // Query in the reverse direction
    MapFragReverseVec query_reverse(const IntPairVec& bounds);

    MapFragForwardVec frags_forward_;
    MapFragReverseVec frags_reverse_;

    void dump_frags(std::ostream& os);

  };

  template< class T>
  std::vector<T> MapFragDB::query_frags(int lb, int ub, const std::vector<T>& frags) {
    typename std::vector<T>::const_iterator lbi = lower_bound(frags.begin(), frags.end(), lb);
    typename std::vector<T>::const_iterator ubi = upper_bound(frags.begin(), frags.end(), ub);
    std::vector<T> ret(lbi, ubi);
    return ret;
  }

  template< class T>
  std::vector<T> MapFragDB::query_frags(const IntPair& ip, const std::vector<T>& frags) {
    return query_frags<T>(ip.first, ip.second, frags);
  }

  template< class T>
  std::vector<T> MapFragDB::query_frags(const IntPairVec& bounds, const std::vector<T>& frags) {

    if (bounds.empty()) {
      return std::vector<T>();
    }

    std::vector<T> active_set = query_frags<T>(bounds[0], frags);
    typename std::vector<T>::iterator s_iter(active_set.begin()), e_iter(active_set.end());

    for(int i = 1; i < int(bounds.size()); i++) {

      if (e_iter == s_iter) break;

      // For items that remain, mark as invalid if the size is not in range
      int lb = bounds[i].first;
      int ub = bounds[i].second;

      // Get next element for every item in active_set
      for(auto frag_iter = s_iter; frag_iter != e_iter; frag_iter++) {
        frag_iter->operator++();
        int size = frag_iter->get_size();
        frag_iter->valid_ = frag_iter->valid_ && (size >= lb && size <= ub);
      }

      // Filter out elements that are no longer valid
      //e_iter = std::remove_if(s_iter, e_iter, map_frag_is_invalid);
      e_iter = std::remove_if(s_iter, e_iter, FragIsInvalid<T>());
      
    }
    active_set.erase(e_iter, active_set.end());

    // "Rewind" the MapFrag's back to the beginning of the match
    for(auto frag_iter = s_iter; frag_iter != e_iter; frag_iter++) {
      frag_iter->operator+=(-(bounds.size() - 1));
    }
    return active_set;
  }

  template< class T>
  std::vector<T> MapFragDB::query_frags_start_largest(const IntPairVec& bounds, const std::vector<T>& frags) {

    if (bounds.empty()) {
      return std::vector<T>();
    }

    // Start the search with the largest fragment, which should have the smallest number of hits since
    // large fragments are less prevalent.
    size_t start_ind = 0;
    int max_bound = -1;
    for(auto i = bounds.begin(); i != bounds.end(); i++) {
      if (i->second > max_bound) {
        max_bound = i->second;
        start_ind = i - bounds.begin();
      }
    }

    std::vector<T> active_set = query_frags<T>(bounds[start_ind], frags);
    typename std::vector<T>::iterator s_iter(active_set.begin()), e_iter(active_set.end());

    // ******************************************
    // Search forward from the starting bound
    for(int i = start_ind + 1; i < int(bounds.size()); i++) {

      if (e_iter == s_iter) break;

      // For items that remain, mark as invalid if the size is not in range
      int lb = bounds[i].first;
      int ub = bounds[i].second;

      // Get next element for every item in active_set
      for(auto frag_iter = s_iter; frag_iter != e_iter; frag_iter++) {
        frag_iter->operator++();
        int size = frag_iter->get_size();
        frag_iter->valid_ = frag_iter->valid_ && (size >= lb && size <= ub);
      }

      // Filter out elements that are no longer valid
      //e_iter = std::remove_if(s_iter, e_iter, map_frag_is_invalid);
      e_iter = std::remove_if(s_iter, e_iter, FragIsInvalid<T>());
      
    }

    // Start by rewinding the MapFrags to the starting index of bounds.
    for(auto frag_iter = s_iter; frag_iter != e_iter; frag_iter++) {
      frag_iter->operator+=(-(bounds.size() - 1 - start_ind));
    }

    // Search backwards
    for(int i = start_ind - 1; i >= 0; i--) {

      if(e_iter == s_iter) break;

      int lb = bounds[i].first;
      int ub = bounds[i].second;

      // Get next element for every item in active_set
      for(auto frag_iter = s_iter; frag_iter != e_iter; frag_iter++) {
        frag_iter->operator--();
        int size = frag_iter->get_size();
        frag_iter->valid_ = frag_iter->valid_ && (size >= lb && size <= ub);
      }

      // Filter out elements that are no longer valid
      //e_iter = std::remove_if(s_iter, e_iter, map_frag_is_invalid);
      e_iter = std::remove_if(s_iter, e_iter, FragIsInvalid<T>());

    }

    // All of the MapFrag's that remain in the active_set should be at 0.

    active_set.erase(e_iter, active_set.end());
    return active_set;
    
  }
}

#endif
