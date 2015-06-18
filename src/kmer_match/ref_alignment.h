#ifndef REFALIGNMENT_H
#define REFALIGNMENT_H

#include <vector>
#include <exception>
#include <algorithm>

#include "map_chunk.h"

namespace kmer_match {
  
  // Define ReferenceAlignment Class 
  // This only represents the "reference" portion of the alignment.
  // The chunks are oriented in same orientation as the query.
  class ReferenceAlignment {

  public:

    // An alignment is given by chunks of the reference
    ReferenceAlignment() {};
    ReferenceAlignment(const ConstMapChunkPVec& chunks) : chunks_(chunks) { compute_summary(); }
    ReferenceAlignment(const ConstMapChunkPVec&& chunks) : chunks_(chunks) { compute_summary(); }

    // Constructor for an alignment comprising a single reference chunk.
    ReferenceAlignment(const MapChunk* chunk) : chunks_(1, chunk) { compute_summary(); }

    template <class ChunkIter>
    ReferenceAlignment(ChunkIter begin, ChunkIter end) : chunks_(begin, end) { compute_summary(); }

    // Get the reference map of the alignment.
    const Map * get_map() const {

      if (chunks_.size() > 0) {
        return chunks_[0]->get_map();
      }

      return nullptr;
    }

    // Alignment is to forward strand of reference
    bool is_forward() const {

      if (chunks_.empty()) {
        return true; // arbitrary response
      }

      const MapChunk * first = chunks_.front();
      const MapChunk * last = chunks_.back();
      return first->start_ < last->start_;

    }

    // Alignment is to reverse strand of reference
    bool is_reverse() const {
      return !is_forward();
    }

    // Reverse the alignment
    void reverse() {
      std::reverse(chunks_.begin(), chunks_.end());
      compute_summary();
    }

    // Count the number of reference misses, number of matches, and
    // the miss rate
    void compute_summary();

    template <class Iter>
    void prepend_chunks(Iter start, Iter end);

    void prepend_chunk(const MapChunk*);

    template <class Iter>
    void append_chunks(Iter start, Iter end);

    void append_chunk(const MapChunk*);

    void check_sane() const; // Check that the reference alignment is sane (i.e all chunks from the same map and in proper order?)

    ConstMapChunkPVec chunks_;
    int start_;
    int end_;
    int num_misses_;
    int num_matches_;
    double miss_rate_;

  };

  template <class Iter>
  inline void ReferenceAlignment::prepend_chunks(Iter start, Iter end) {

    int delta = end - start;

    if (delta > 0) {

      chunks_.resize(chunks_.size() + delta);

      // copy existing chunks right
      std::copy(chunks_.rbegin() + delta, chunks_.rend(), chunks_.rbegin());

      // copy new chunks to beginning
      std::copy(start, end, chunks_.begin());

      compute_summary();

    }
  }

  inline void ReferenceAlignment::prepend_chunk(const MapChunk* chunk) {

      chunks_.resize(chunks_.size() + 1);

      // copy existing chunks right
      std::copy(chunks_.rbegin() + 1, chunks_.rend(), chunks_.rbegin());

      chunks_[0] = chunk;

      compute_summary();

  }



  template <class Iter>
  inline void ReferenceAlignment::append_chunks(Iter start, Iter end) {
    
    int delta = end - start;

    if (delta > 0) {

      const int orig_size = chunks_.size();
      chunks_.resize(orig_size + delta);

      // copy new chunks to end
      std::copy(start, end, chunks_.begin() +  orig_size);

      compute_summary();

    }
  }

  inline void ReferenceAlignment::append_chunk(const MapChunk* chunk) {
    chunks_.push_back(chunk);
    compute_summary(); 
  }

  inline void ReferenceAlignment::compute_summary() {

    start_ = -1;
    end_ = -1;
    num_misses_ = -1;
    num_matches_ = -1;
    miss_rate_ = -1.0;

    if (chunks_.empty()) { return; }

    num_matches_ = int(chunks_.size() + 1);
    start_ = chunks_.front()->start_;
    end_ = chunks_.back()->end_;
    num_misses_ = int((end_ - start_) - chunks_.size());
    miss_rate_ = double(num_misses_)/double(num_misses_ + num_matches_);

  }

  class ReferenceAlignmentException: public std::exception
  {
    virtual const char* what() const throw()
    {
      return "Reference Alignment is not sane.";
    }
  };

  inline void ReferenceAlignment::check_sane() const {

    if (chunks_.empty()) return;

    bool aln_is_forward = is_forward();
    const Map * pMap = chunks_[0]->get_map();
    size_t last_start = chunks_[0]->start_;
    size_t last_end = chunks_[0]->end_;

    if(aln_is_forward) {

      for(auto i = chunks_.begin() + 1; i != chunks_.end(); i++) {
        const MapChunk* c = *i;
        if ( (c->get_map() != pMap) || (c->start_ != last_end) ) {
          throw ReferenceAlignmentException();
        }

      }

    } else {

      for(auto i = chunks_.begin() + 1; i != chunks_.end(); i++) {
        const MapChunk* c = *i;
        if ( (c->get_map() != pMap) || (c->end_ != last_start) ) {
          throw ReferenceAlignmentException();
        }

      }

    }

  }


  struct ReferenceAlignmentMissRateSort {
    bool operator()(const ReferenceAlignment& l, const ReferenceAlignment& r) {
       return l.miss_rate_ < r.miss_rate_;
    }
  };


  std::ostream& operator<<(std::ostream& os, const ReferenceAlignment& aln);

  typedef std::vector<ReferenceAlignment> RefAlignmentVec;

}

#endif