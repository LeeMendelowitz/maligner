#include "timer.h"
#include <ostream>

namespace lmm_utils {
  
  std::ostream& operator<<(std::ostream& os, const Timer& t) {
    using namespace std::chrono;
    Timer::time_point end_time = t.has_ended_ ? t.end_time_ : steady_clock::now();
    auto elapsed = end_time - t.start_time_;
    os << "Timer " << t.name_ << ": "
       <<  duration_cast<milliseconds>(elapsed).count()
       << " ms elapsed.";
    return os;
  }

}