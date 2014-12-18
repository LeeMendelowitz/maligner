#include <iostream>
#include <sstream>
#include <fstream>

#include "map_frag.h"
#include "map.h"

using namespace std;


bool operator<(int val, const MapFragForward& frag) {
  return val < frag.get_size();
}

bool operator<(int val, const MapFragReverse& frag) {
  return val < frag.get_size();
}

///////////////////////////////////////////////////////
// Print a MapFrag to ostream.
std::ostream& operator<<(std::ostream& os, const MapFragForward& m) {
  os << m.pMap_ << ": "
     << m.get_size() << "bp " 
     << m.index_ <<  " "
     << "Forward "
     << (m.valid_ ? "VALID" : "INVALID");
  return os;

}

std::ostream& operator<<(std::ostream& os, const MapFragReverse& m) {
  os << m.pMap_ << ": "
     << m.get_size() << "bp " 
     << m.index_ <<  " "
     << "Reverse "
     << (m.valid_ ? "VALID" : "INVALID");
  return os;

}

// Write the context of the map_frag
std::ostream& MapFragForward::write_context(std::ostream& os, size_t N) const {

  os << "(" << pMap_->name_ << "," << index_ << ","
     << "Forward" << ") ";
  
  std::ostream_iterator<int> out_it(os, " ");

  IntVec::const_iterator s = pMap_->frags_.begin() + index_;
  IntVec::const_iterator e = s + N;
  if (e > pMap_->frags_.end()) { 
    e = pMap_->frags_.end(); 

  }
  std::copy(s, e, out_it);

  return os;

}

std::ostream& MapFragReverse::write_context(std::ostream& os, size_t N) const {

  os << "(" << pMap_->name_ << "," << index_ << ","
     << "Reverse" << ") ";
  
  std::ostream_iterator<int> out_it(os, " ");
  IntVec::const_reverse_iterator s = pMap_->frags_.rend() - index_ - 1;
  IntVec::const_reverse_iterator e = s + N;
  if ( e > pMap_->frags_.rend() ) {
    e = pMap_->frags_.rend();
  }

  std::copy(s, e, out_it);


  return os;

}