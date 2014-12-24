#include <iostream>
#include "map_reader.h"

using namespace std;
using namespace kmer_match;

bool MapReader::next(Map& map) {

  std::string line;

  while(std::getline(f_, line)) {

    try {
        map = std::move(Map(line));
        return true;
    } catch(ReadMapException& e) {
        cerr << e.what() << "\n";
    }

  }

  return false;

}

MapVec MapReader::read_all_maps() {
  MapVec maps;
  Map map;
  while(next(map)) {
    maps.push_back(std::move(map));
  }
  return maps;
}