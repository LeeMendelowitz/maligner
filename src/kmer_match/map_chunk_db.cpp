#include <utility>
#include <algorithm>
#include <iterator>
#include <cassert>

#include "map_chunk_db.h"
#include "ref_alignment.h"

using namespace std;



MapChunkDB::MapChunkDB(const MapVec& maps, size_t frags_per_chunk) :
frags_per_chunk_(frags_per_chunk) {
  
  // Let's just ensure we have enough space.
  size_t num_frags = 0;

  for(auto mi = maps.begin(); mi != maps.end(); mi++) {
    num_frags += mi->frags_.size();
  }

  // Reserve enough space to avoid reallocation.
  map_chunks_.reserve(num_frags*frags_per_chunk);

  // For each map, compute the chunks.
  for(auto mi = maps.begin(); mi != maps.end(); mi++) {

    const Map * pMap = &(*mi);

    ChunksAtIndex chunks_at_index_start_(mi->frags_.size());
    ChunksAtIndex chunks_at_index_end_(mi->frags_.size());
    
    const size_t num_frags = mi->frags_.size();
    for(size_t start = 0; start < num_frags; start++) {
      const size_t last = min(start + frags_per_chunk, num_frags - 1);
      for(size_t end = start + 1; end <= last; end++) {
        map_chunks_.emplace_back(pMap, start, end);
        chunks_at_index_start_[start].push_back(&map_chunks_.back());
        chunks_at_index_end_[end].push_back(&map_chunks_.back());
      } 
    }

    map_to_chunks_at_start_[pMap] = std::move(chunks_at_index_start_);
    map_to_chunks_at_end_[pMap] = std::move(chunks_at_index_end_);

  }

  sort_chunks();

}

class MapChunkPCmp {
public:
  bool operator()(const MapChunk* c1, const MapChunk* c2) {
    return c1->size_ < c2->size_;
  }
};

class IntMapChunkPCmp {
public:
  bool operator()(int s, const MapChunk* c2) const {
    return s < c2->size_;
  }
};

class MapChunkPIntCmp {
public:
  bool operator()(const MapChunk* c2, int s) const {
    return c2->size_ < s;
  }
};


void MapChunkDB::sort_chunks() {
  map_chunk_index_ = MapChunkConstPVec(map_chunks_.size());
  
  size_t j = 0;
  for(vector<MapChunk>::const_iterator i = map_chunks_.begin();
      i != map_chunks_.end();
      i++, j++) {
    map_chunk_index_[j] = &(*i);
  }

  sort(map_chunk_index_.begin(), map_chunk_index_.end(), MapChunkPCmp());

}

MapChunkVecConstIterPair MapChunkDB::query(int lb, int ub) const {
  MapChunkVecConstIter lbi = lower_bound(map_chunk_index_.begin(), map_chunk_index_.end(), lb, MapChunkPIntCmp());
  MapChunkVecConstIter ubi = upper_bound(map_chunk_index_.begin(), map_chunk_index_.end(), ub, IntMapChunkPCmp());
  return MapChunkVecConstIterPair(lbi, ubi);
}

MapChunkVecConstIterPair MapChunkDB::query(const IntPair& ip) const {
  return query(ip.first, ip.second);
}


int MapChunkDB::count_compatible_seeds(const IntPairVec& bounds) const {

  // See if there is a sequence of fragments in the database within the error bounds.
  // For now we are only matching in the forward direction.

  if (bounds.empty()) {
    return true; // or false?
  }

  // Find the largest fragment the bounds.
  int max_bound = -1;
  IntPairVec::const_iterator iter_max_bound;
  for(IntPairVec::const_iterator i = bounds.begin(); i != bounds.end(); i++) {
    if (i->second > max_bound) {
      max_bound = i->second;
      iter_max_bound = i;
    }
  }

  MapChunkVecConstIterPair matches = query(*iter_max_bound);

  // For each match to the largest fragment, test if it is compatible for the other
  // bounds.
  int num_matches = 0;
  typedef reverse_iterator<IntPairVec::const_iterator> RevIter;
  for(auto mi = matches.first; mi != matches.second; mi++) {
    
    // bool have_forward = is_compatible_forward(*mi, iter_max_bound + 1, bounds.end());
    RevIter rev_start(iter_max_bound - 1);
    RevIter rev_end(bounds.begin());
    // bool have_reverse = is_compatible_reverse(*mi, start, end);

    // Match in the forward direction
    int c1 = count_compatible_right_bfs(*mi, iter_max_bound + 1, bounds.end());
    int c2 = count_compatible_left_bfs(*mi, rev_start, rev_end);

    // Match in the reverse direction
    int c3 = count_compatible_right_bfs(*mi, rev_start, rev_end);
    int c4 = count_compatible_left_bfs(*mi, iter_max_bound + 1, bounds.end());

    num_matches += c1*c2;
    num_matches += c3*c4;

  }

  return num_matches;

}

///////////////////////////////////////////////////////////////////////////
// Helper function to merge the left alignments and the right alignments
void merge_alignments(const AlignmentVector& left_aln, const AlignmentVector& right_aln,
    const MapChunk* middle_chunk,
    AlignmentVector& merged) {

  using Size_T = AlignmentVector::size_type;

  if (left_aln.empty() && right_aln.empty()) return;

  size_t nl = max(left_aln.size(), Size_T(1));
  size_t nr = max(right_aln.size(), Size_T(1));
  merged.reserve(merged.size() + nl*nr);

  if(left_aln.empty()) {

    // Only merge the middle chunk with each right alignment.
    for(const auto& aln: right_aln) {

      const Size_T N(1 + aln.size());
      Alignment new_aln(N);
      new_aln[0] = middle_chunk;
      copy(aln.begin(), aln.end(), new_aln.begin() + 1);
      merged.push_back(move(new_aln));

    }
    return;

  }

  if(right_aln.empty()) {

    // Only merge the left alignment with the middle chunk.
    // Orient forwards.
    for(const auto& aln: left_aln) {

      Alignment new_aln(aln.rbegin(), aln.rend());
      new_aln.push_back(middle_chunk);
      merged.push_back(move(new_aln));

    }
    return;
  }

  // Both left_aln and right_aln are not empty, so we must merge them.
  for(const auto& l : left_aln) {

    Alignment new_aln_left(l.rbegin(), l.rend()); // copy and reverse
    new_aln_left.push_back(middle_chunk);
    const Size_T new_aln_size(new_aln_left.size());
    for(auto& r: right_aln) {
      Alignment new_aln(new_aln_left);
      new_aln.resize(new_aln_size + r.size());
      copy(r.begin(), r.end(), new_aln.begin() + new_aln_size);
      merged.push_back(move(new_aln));
    }
  }
}

///////////////////////////////////////////////////////////////////////////
// Helper function to merge the left alignments and the right alignments
// Merge only the best left alignment with the best right alignment.
ReferenceAlignment merge_alignments_select_best(AlignmentVector&& left_aln, AlignmentVector&& right_aln,
    const MapChunk* middle_chunk) {

  using Size_T = AlignmentVector::size_type;

  if (left_aln.empty() && right_aln.empty()) return ReferenceAlignment();

  size_t nl = max(left_aln.size(), Size_T(1));
  size_t nr = max(right_aln.size(), Size_T(1));

  ReferenceAlignmentMissRateSort ref_aln_sorter;

  if(left_aln.empty()) {

    // Select the best right alignment, merge the middle chunk, and return.

    RefAlignmentVec right_ref_alns;
    right_ref_alns.reserve(right_aln.size());
    for(AlignmentVector::iterator iter = right_aln.begin(); iter != right_aln.end(); iter++){
      right_ref_alns.emplace_back(std::move(*iter));
    }

    // Sort the right alignments by miss rate
    std::sort(right_ref_alns.begin(), right_ref_alns.end(), ref_aln_sorter);

    // Construct a new reference alignment by prepending the middle chunk to the chunks of the best right alignment.
    ReferenceAlignment& best = right_ref_alns[0];
    best.prepend_chunk(middle_chunk);

    return best;
  }

  if(right_aln.empty()) {

    // Select the best left alignment, merge the middle chunk, and return.

    RefAlignmentVec left_ref_alns;
    left_ref_alns.reserve(left_aln.size());
    for(AlignmentVector::iterator iter = left_aln.begin(); iter != left_aln.end(); iter++){
      left_ref_alns.emplace_back(std::move(*iter));
    }

    // Sort the left alignments by miss rate
    std::sort(left_ref_alns.begin(), left_ref_alns.end(), ref_aln_sorter);

    // Construct a new reference alignment by prepending the middle chunk to the chunks of the best right alignment.
    ReferenceAlignment& best = left_ref_alns[0];
    best.reverse(); // Reverse to orient forward
    best.append_chunk(middle_chunk);
    
    return best;

  }

  // Select the best right alignment
  RefAlignmentVec right_ref_alns;
  right_ref_alns.reserve(right_aln.size());
  for(AlignmentVector::iterator iter = right_aln.begin(); iter != right_aln.end(); iter++){
    right_ref_alns.emplace_back(std::move(*iter));
  }
  std::sort(right_ref_alns.begin(), right_ref_alns.end(), ref_aln_sorter);
  ReferenceAlignment& best_right = right_ref_alns[0];

  // Select the best left alignment
  RefAlignmentVec left_ref_alns;
  left_ref_alns.reserve(left_aln.size());
  for(AlignmentVector::iterator iter = left_aln.begin(); iter != left_aln.end(); iter++){
    left_ref_alns.emplace_back(std::move(*iter));
  }
  std::sort(left_ref_alns.begin(), left_ref_alns.end(), ref_aln_sorter);
  ReferenceAlignment& best_left = left_ref_alns[0];
  best_left.reverse(); // reverse the left alignment to orient forward with respect to the reference.

  // Merge the middle chunk and right alignment with the left alignment
  ReferenceAlignment& best_aln = best_left;
  best_aln.append_chunk(middle_chunk);
  best_aln.append_chunks(best_right.chunks_.begin(), best_right.chunks_.end());

  return best_aln;

}


AlignmentVector MapChunkDB::get_compatible_alignments(const IntPairVec& bounds, size_t max_unmatched) const {

  // See if there is a sequence of fragments in the database within the error bounds.
  // For now we are only matching in the forward direction.

  if (bounds.empty()) {
    return AlignmentVector();
  }

  // Find the largest fragment of the bounds.
  // We seed on the largest fragment.
  int max_bound = -1;
  IntPairVec::const_iterator iter_max_bound;
  for(IntPairVec::const_iterator i = bounds.begin(); i != bounds.end(); i++) {
    if (i->second > max_bound) {
      max_bound = i->second;
      iter_max_bound = i;
    }
  }


  bool must_search_left = (iter_max_bound != bounds.begin());
  bool must_search_right = (iter_max_bound != bounds.end() - 1);

  // For reverse searches, the query is oriented in reverse direction.
  // Left refers to left in the reference, right refers to right in the reference.
  bool must_search_left_reverse = must_search_right;
  bool must_search_right_reverse = must_search_left;

  MapChunkVecConstIterPair matches = query(*iter_max_bound);

  if(bounds.size() == 1) {
    // In this case, the hits are simply the range of MapChunk *'s returned by the query function
    AlignmentVector ret;
    ret.reserve(matches.second - matches.first);
    for(auto mi = matches.first; mi != matches.second; ++mi)  {
      const MapChunk* middle_chunk = *mi;
      ret.emplace_back(1, middle_chunk);
    }
    return ret;
  }

  // For each match to the largest fragment, test if it is compatible for the other
  // bounds.
  int num_matches = 0;
  typedef reverse_iterator<IntPairVec::const_iterator> RevIter;

  AlignmentVector all_alignments;
  
  bool have_alignments_forward;
  bool have_alignments_left;

  for(auto mi = matches.first; mi != matches.second; ++mi) {
    
    const MapChunk* middle_chunk = *mi;

    if (middle_chunk->num_unmatched() > max_unmatched) continue;
    size_t extension_max_unmatched = max_unmatched - middle_chunk->num_unmatched();

    RevIter rev_start(iter_max_bound); // Starts 1 to left of iter_max_bound
    RevIter rev_end(bounds.begin()); // Ends 1 to left of begin (exclusive).


    ///////////////////////////////////////////////////////////////////////////////////
    // Match in the forward direction
    {
      AlignmentVector aln_right_forward, aln_left_forward; 
      bool have_alignments_forward = true;

      if(must_search_right) {
        aln_right_forward = get_compatible_right_dfs(middle_chunk, iter_max_bound + 1, bounds.end(), extension_max_unmatched);
        have_alignments_forward = !aln_right_forward.empty();
      }

      if(have_alignments_forward && must_search_left) {
        aln_left_forward = get_compatible_left_dfs(middle_chunk, rev_start, rev_end, extension_max_unmatched);
        have_alignments_forward = !aln_left_forward.empty();
      }

      if(have_alignments_forward) {
        merge_alignments(aln_left_forward, aln_right_forward, middle_chunk, all_alignments);
      }
    } // end match in the forward direction

    ///////////////////////////////////////////////////////////////////////////////////
    // Match in the reverse direction
    {

      bool have_alignments_reverse = true;
      AlignmentVector aln_right_reverse, aln_left_reverse; 

      if(must_search_right_reverse) {
        aln_right_reverse = get_compatible_right_dfs(*mi, rev_start, rev_end, extension_max_unmatched);
        have_alignments_reverse = !aln_right_reverse.empty();
      }

      if(have_alignments_reverse && must_search_left_reverse) {
        aln_left_reverse = get_compatible_left_dfs(*mi, iter_max_bound + 1, bounds.end(), extension_max_unmatched);
        have_alignments_reverse = !aln_left_reverse.empty();
      }

      if(have_alignments_reverse) {
          merge_alignments(aln_right_reverse, aln_left_reverse, middle_chunk, all_alignments);
      }

    } // end match in the reverse direction

  }

  return all_alignments;

}

/////////////////////////////////////////////////////////////////////////////////////
// Get the best alignment in the forward and reverse direction for each seed hit.
// Unlike get_compatbile_alignments, which returns all alignments in forward/reverse direction for a seed,
// this will return at most two alignments per seed (one for forward direction, the other for reverse).
RefAlignmentVec MapChunkDB::get_compatible_alignments_best(const IntPairVec& bounds, size_t max_unmatched) const {

  // See if there is a sequence of fragments in the database within the error bounds.

  if (bounds.empty()) {
    return RefAlignmentVec();
  }

  // Find the largest fragment of the bounds.
  // We seed on the largest fragment.
  int max_bound = -1;
  IntPairVec::const_iterator iter_max_bound;
  for(IntPairVec::const_iterator i = bounds.begin(); i != bounds.end(); i++) {
    if (i->second > max_bound) {
      max_bound = i->second;
      iter_max_bound = i;
    }
  }


  bool must_search_left = (iter_max_bound != bounds.begin());
  bool must_search_right = (iter_max_bound != bounds.end() - 1);

  // For reverse searches, the query is oriented in reverse direction.
  // Left refers to left in the reference, right refers to right in the reference.
  bool must_search_left_reverse = must_search_right;
  bool must_search_right_reverse = must_search_left;

  MapChunkVecConstIterPair matches = query(*iter_max_bound);

  if(bounds.size() == 1) {

    // In this case, the hits are simply the range of MapChunk *'s returned by the query function
    RefAlignmentVec ret;
    ret.reserve(matches.second - matches.first);
    for(auto mi = matches.first; mi != matches.second; ++mi)  {
      const MapChunk* middle_chunk = *mi;
      ret.emplace_back(middle_chunk);
    }
    return ret;
  }

  // For each match to the largest fragment, test if it is compatible for the other
  // bounds.
  int num_matches = 0;
  typedef reverse_iterator<IntPairVec::const_iterator> RevIter;

  RefAlignmentVec all_alignments;
  
  bool have_alignments_forward;
  bool have_alignments_left;

  for(auto mi = matches.first; mi != matches.second; ++mi) {
    
    const MapChunk* middle_chunk = *mi;

    if (middle_chunk->num_unmatched() > max_unmatched) continue;
    size_t extension_max_unmatched = max_unmatched - middle_chunk->num_unmatched();

    RevIter rev_start(iter_max_bound); // Starts 1 to left of iter_max_bound
    RevIter rev_end(bounds.begin()); // Ends 1 to left of begin (exclusive).


    ///////////////////////////////////////////////////////////////////////////////////
    // Match in the forward direction
    {
      AlignmentVector aln_right_forward, aln_left_forward; 
      bool have_alignments_forward = true;

      if(must_search_right) {
        aln_right_forward = get_compatible_right_dfs(middle_chunk, iter_max_bound + 1, bounds.end(), extension_max_unmatched);
        have_alignments_forward = !aln_right_forward.empty();
      }

      if(have_alignments_forward && must_search_left) {
        aln_left_forward = get_compatible_left_dfs(middle_chunk, rev_start, rev_end, extension_max_unmatched);
        have_alignments_forward = !aln_left_forward.empty();
      }

      if(have_alignments_forward) {
        ReferenceAlignment best_forward = merge_alignments_select_best(std::move(aln_left_forward), std::move(aln_right_forward), middle_chunk);
        all_alignments.push_back(std::move(best_forward));
      }
    } // end match in the forward direction

    ///////////////////////////////////////////////////////////////////////////////////
    // Match in the reverse direction
    {

      bool have_alignments_reverse = true;
      AlignmentVector aln_right_reverse, aln_left_reverse; 

      if(must_search_right_reverse) {
        aln_right_reverse = get_compatible_right_dfs(middle_chunk, rev_start, rev_end, extension_max_unmatched);
        have_alignments_reverse = !aln_right_reverse.empty();
      }

      if(have_alignments_reverse && must_search_left_reverse) {
        aln_left_reverse = get_compatible_left_dfs(middle_chunk, iter_max_bound + 1, bounds.end(), extension_max_unmatched);
        have_alignments_reverse = !aln_left_reverse.empty();
      }

      if(have_alignments_reverse) {
          ReferenceAlignment best_reverse = merge_alignments_select_best(std::move(aln_right_reverse), std::move(aln_left_reverse), middle_chunk);
          all_alignments.push_back(std::move(best_reverse));
      }

    } // end match in the reverse direction

  }

  return all_alignments;

}




////////////////////////////////////////////////////////////////////
// NOTE: If we know the orientation of the alignment, we can avoid
// some of the conditional branching in this function by writing
// separate print functions for forward and reverse alignments.
std::ostream& operator<<(std::ostream& os, const Alignment& aln) {
  
  // if(aln.empty()) {
  //   return os;
  // }

  const MapChunk * first = aln.front();
  const MapChunk * last = aln.back();
  bool is_forward = first->start_ < last->start_;
  
  os << first->pMap_->name_ << " ";

  if(is_forward) {
    os << "F " << first->start_ << " " << last->end_ << " ";
  } else {
    os << "R " << last->start_ << " " << first->end_ << " ";
  }

  for(Alignment::const_iterator iter = aln.begin(); iter != aln.end(); iter++) {
    os << **iter << " ";
  }

  return os;
}