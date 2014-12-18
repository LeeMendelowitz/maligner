#ifndef MAP_H
#define MAP_H

#include <string>
#include <vector>
#include <iostream>
#include <iterator>

typedef std::vector<int> IntVec;
typedef int FragSize;
typedef std::vector<FragSize> FragVec;

class Map {

  public:

  ////////////////////////////////////////////////
  // METHODS

  // Default constructor  
  Map(const std::string& name, int size, const IntVec& frags);
  Map(const std::string& line);

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

std::ostream& operator<<(std::ostream& os, const Map& m);

MapVec read_maps(const std::string& file_name);



#endif