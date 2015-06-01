#include <random>
#include <algorithm>

#include "map_sampler.h"

namespace maligner_maps {


  FragVec sample_with_replacement(const FragVec& frags, const size_t n) {
    
    FragVec ret(n);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, frags.size());

    for(size_t i = 0; i < n; i++) {
      ret[i] = frags[dis(gen)];
    }

    return ret;

  }


  FragVec permute(const FragVec& frags) {
    
    FragVec ret(frags);
    random_shuffle(ret.begin(), ret.end());
    return ret;

  }


}