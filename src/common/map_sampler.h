#ifndef MAP_SAMPLER_H
#define MAP_SAMPLER_H

#include <vector>
#include <random>
#include <algorithm>
#include <utility>

// Utilities for permuting maps and sampling fragments with replacement

// Sample n fragments with replacement
template< class T>
std::vector<T> sample_with_replacement(const std::vector<T>& frags, const size_t n) {

    std::vector<T> ret(n);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, frags.size());

    for(size_t i = 0; i < n; i++) {
      ret[i] = frags[dis(gen)];
    }

    return ret;
}

// Return a permuted copy of frags
template< class T>
std::vector<T> permute(const std::vector<T>& frags) {

    std::vector<T> ret(frags);
    std::random_device rd;
    std::mt19937 gen(rd());

    const size_t n = ret.size();
    for(size_t i = 0; i < n; i++) {
        std::uniform_int_distribution<> dis(i, n-1);
        std::swap(ret[i], ret[dis(gen)]);
    }
    return ret;
}


#endif