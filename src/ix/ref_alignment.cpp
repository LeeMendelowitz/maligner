#include "ref_alignment.h"


namespace kmer_match {
  
  std::ostream& operator<<(std::ostream& os, const ReferenceAlignment& aln) {

    const MapChunk * first = aln.chunks_.front();
    const MapChunk * last = aln.chunks_.back();
    bool is_forward = first->start_ < last->start_;
    
    os << first->get_map()->name_ << " ";

    if(is_forward) {
      os << "F " << first->start_ << " " << last->end_ << " ";
    } else {
      os << "R " << last->start_ << " " << first->end_ << " ";
    }

    for(ConstMapChunkPVec::const_iterator iter = aln.chunks_.begin(); iter != aln.chunks_.end(); iter++) {
      os << **iter << " ";
    }

    return os;
  }

}