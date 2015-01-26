#ifndef PARTIALSUMS_H
#define PARTIALSUMS_H

namespace maligner_dp {

  class PartialSums {
    // This precomputes chunk lengths to assist in speeding up dynamic programming.
    // It behaves as a 2D array with chunk lengths, indexed the following way:
    // ps(i, j): The length of the chunk of j fragments ending with fragment i.
    //
    // Example:
    //     -----------|-----------|-----------|-----|-----
    //      f[i-2]      f[i-1]     f[i]
    //
    // ps(i, 0) = f[i]
    // ps(i, 1) = f[i] + f[i-1]
    // ps(i, 2) = f[i] + f[i-1] + f[i-2]
    //
    // This indexing may be counterintuitive, but it works more seamlessly with the dynamic programming

    friend class SDInv;
    IntVec d_;
    const size_t m_;

  public:

    struct forward_tag {};
    struct reverse_tag {};

    PartialSums(const IntVec& d, int max_miss) : m_(max_miss + 1) {
      fill_forward(d, max_miss); 
    }

    PartialSums(const IntVec& d, int max_miss, forward_tag) : m_(max_miss + 1) {
      fill_forward(d, max_miss);
    }

    PartialSums(const IntVec& d, int max_miss, reverse_tag) : m_(max_miss + 1) {
      fill_reverse(d, max_miss);
    }

    int operator()(size_t i, size_t num_miss) const {
      return d_[i*m_ + num_miss];
    }

  private:

    void fill_forward(const IntVec& d, int max_miss) {
      const size_t n = d.size();
      d_ = IntVec(n*m_);

      for(size_t s = 0; s < n; s++) {

        size_t offset = s*m_;
        
        // Constraints: s - i >= 0 -> i <= s -> i < s + 1
        // i < m_
        const size_t L = min(m_, s + 1);
        int a = 0;
        for(size_t i = 0; i < L; i++) {
          a += d[s - i];
          d_[offset+i] = a;
        }

      }

    }

    void fill_reverse(const IntVec& d, int max_miss) {
      const size_t n = d.size();
      d_ = IntVec(n*m_);

      for(size_t s = 0; s < n; s++) {

        const size_t sr = n - s - 1; // Extract right to left from d
        size_t offset = s*m_; // but fill d_ left to right
      
        // Constraints: sr - i >= 0 -> i <= sr -> i < sr + 1
        // i < m_
        const size_t L = min(m_, sr + 1);
        int a = 0;
        for(size_t i = 0; i < L; i++) {
          a += d[sr - i];
          d_[offset+i] = a;
        }

      }

    }

  };

  class SDInv {

    // Precompute 1/(sd) for each reference chunk

    public:
      SDInv(const PartialSums& ps, double sd_rate, double min_sd) :
        d_(ps.d_.size()),
        m_(ps.m_) {

          const size_t n = d_.size();
          for(size_t i = 0; i < n; i++) {
            int ref_size = ps.d_[i];
            double sd = sd_rate * ref_size;
            if (sd < min_sd) { sd = min_sd; }
            double sd_1 = 1.0/sd;
            d_[i] = sd_1;
          }

      }     
      
      double operator()(size_t i, size_t num_miss) const {
        return d_[i*m_ + num_miss];
      }

    DoubleVec d_;
    const size_t m_;
  };
  
}

#endif
