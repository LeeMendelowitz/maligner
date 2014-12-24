#ifndef LMM_TIMER_H
#define LMM_TIMER_H

#include <string>
#include <chrono>
#include <ostream>

namespace lmm_utils {

  class Timer {

  public:
    typedef std::chrono::steady_clock::time_point time_point;

    Timer() : has_ended_(false) {
      start();
    };

    Timer(const std::string& name) :
      name_(name),
      has_ended_(false)
    { 
        start();
    };

    void start() {
      start_time_ = std::chrono::steady_clock::now();
      end_time_ = start_time_;
      has_ended_ = false;
    }

    void end() {
      end_time_ = std::chrono::steady_clock::now();
      has_ended_ = true;
    }

    void reset() { // alias for start
      start();
    }

    void restart() { // alias
      start();
    }

  private:

    

    std::string name_;
    time_point start_time_;
    time_point end_time_;
    bool has_ended_; // Affects the printing of the Timer.

    friend std::ostream& operator<<(std::ostream&, const Timer& t);
 
  };


}

#endif
