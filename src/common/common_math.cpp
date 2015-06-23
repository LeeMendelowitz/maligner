#include <vector>
#include <math.h>  
#include <algorithm>
#include <cmath>

#include "common_math.h"

// Compute the median absolute deviation
double mad(const std::vector<double>& in) {

  using std::size_t;
  if(in.empty()) return 0.0;

  std::vector<double> v(in);

  const size_t N = v.size();
  double m = median(v);

  for(size_t i = 0; i < N; i++) {
    v[i] = std::abs(v[i] - m);
  }

  return median(v);

}

// Compute the median absolute deviation, if the median m already known.
double mad(const std::vector<double>& in, double m) {

  using std::size_t;
  if(in.empty()) return 0.0;

  std::vector<double> v(in);

  const size_t N = v.size();

  for(size_t i = 0; i < N; i++) {
    v[i] = std::abs(v[i] - m);
  }

  return median(v);

}



// Compute the median
double median(const std::vector<double>& in) {
  
  using std::size_t;
  if(in.empty()) return 0.0;
  std::vector<double> v(in);
  std::sort(v.begin(), v.end());

  const size_t N = v.size();
  
  if( N % 2 == 0) {
    const size_t mid_plus {N/2};
    return 0.5*(v[mid_plus - 1] + v[mid_plus]);
  }

  return v[size_t(N/2)];

}

// Compute the median if the numbers are already sorted.
double median_sorted(const std::vector<double>& in) {

  using std::size_t;
  if(in.empty()) return 0.0;

  const size_t N = in.size();
  
  if( N % 2 == 0) {
    const size_t mid_plus {N/2};
    return 0.5*(in[mid_plus - 1] + in[mid_plus]);
  }

  return in[size_t(N/2)];
}
