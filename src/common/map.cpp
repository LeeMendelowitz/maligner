#include <iostream>
#include <sstream>
#include <fstream>
#include <iterator>
#include <exception>

#include "map.h"
#include "map_reader.h"

using namespace std;


namespace maligner_maps {

  //////////////////////////////////////////////////////
  // Map default constructor
  Map::Map(const std::string& name, int size, const FragVec& frags) :
    name_(name),
    size_(size),
    frags_(frags) {

    }


  /////////////////////////////////////////////////////
  // Construct a map from a line.
  Map::Map(const std::string& line) {

    istringstream iss(line);
    
    iss >> name_;

    if (iss.fail()) throw ReadMapException(line);

    iss >> size_;

    if (iss.fail()) throw ReadMapException(line);

    unsigned int num_frags;
    iss >> num_frags;

    if (iss.fail()) throw ReadMapException(line);

    frags_.clear();
    frags_.reserve(num_frags);

    int frag_size;

    while(iss) {

      // The next operation may fail if we reach the end of the file
      iss >> frag_size;
      if(!iss.fail()) {
        frags_.push_back(frag_size);
      }
    }

    // Check if we failed for some reason other than eof.
    if (iss.fail() && !iss.eof()) {
      throw ReadMapException(line);
    }

    // Check that the number of fragments read matches num_frags:
    if (num_frags != frags_.size()) {
      throw ReadMapException(line);
    }

  }


  //////////////////////////////////////////////////////
  // Read maps from a file. 
  MapVec read_maps(const std::string& file_name) {
    MapReader reader(file_name);
    return reader.read_all_maps();
  }

  ///////////////////////////////////////////////////////
  // Print a map to ostream.
  std::ostream& operator<<(std::ostream& os, const Map& m) {

    os << m.name_ << " "
       << m.size_ << " "
       << m.frags_.size() << " "
       << "(";

    std::ostream_iterator<int> oo{os, " "};
    std::copy(m.frags_.begin(), m.frags_.end(), oo);
    os << ")";

    return os;

  }
}
