#ifndef MAP_H
#define MAP_H

#include <string>
#include <vector>
#include <iostream>
#include <iterator>

namespace maligner_maps {
  
  typedef std::vector<int> IntVec;
  typedef int FragSize;
  typedef std::vector<FragSize> FragVec;

  /////////////////////////////////////////////////////////////////
  // The Map class is a reflection of a map in the .maps file
  // Additional meta-data for a map (such is if the fragments are bounded, or if the map is circular)
  // are provided by command line arguments. Such metadata is stored in MapData.
  class Map {

    public:

    ////////////////////////////////////////////////
    // METHODS

    // Default constructor  
    Map() {};
    Map(const std::string& name, int size, const IntVec& frags);
    Map(const std::string& line);

    Map(const Map& o) :
      name_(o.name_),
      size_(o.size_),
      frags_(o.frags_) {}

    // Move constructor
    Map(Map&& o) :
      name_(std::move(o.name_)),
      size_(std::move(o.size_)),
      frags_(std::move(o.frags_))
    {}
    
    // Copy assignment
    Map& operator=(const Map& o) {
      name_ = o.name_;
      size_ = o.size_;
      frags_ = o.frags_;
      return *this;
    }

    // Move assignment
    Map& operator=(Map&& o) {
      name_ = std::move(o.name_);
      size_ = o.size_;
      frags_ = std::move(o.frags_);
      return *this;
    }

    /////////////////////////////////////////////////
    // MEMBERS

    // name of the map
    std::string name_;

    // size in bp 
    unsigned int size_;

    // vector of fragments 
    FragVec frags_; 

    const FragVec& get_frags() const { return frags_; }

  };


  typedef std::vector<Map> MapVec;  
  typedef std::vector<Map*> MapPVec;

  std::ostream& operator<<(std::ostream& os, const Map& m);

  // MapVec read_maps(const std::string& file_name);

}



#endif