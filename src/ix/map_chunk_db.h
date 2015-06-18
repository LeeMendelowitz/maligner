#ifndef MAP_CHUNK_DB
#define MAP_CHUNK_DB

#include <unordered_map>
#include <utility>
#include <cassert>

#include "map.h"
#include "map_chunk.h"
#include "ref_alignment.h"

namespace kmer_match {

  using namespace maligner_maps;

  typedef std::pair<int, int> IntPair;
  typedef std::vector<IntPair> IntPairVec;

  typedef std::vector<MapChunk*> MapChunkPVec;
  typedef std::vector<const MapChunk*> MapChunkConstPVec;
  typedef std::vector<MapChunkConstPVec> ChunksAtIndex;
  typedef std::vector<const MapChunk*> Alignment;
  typedef std::vector<Alignment> AlignmentVector;
  typedef std::vector<const MapChunk*>::const_iterator MapChunkVecConstIter;
  typedef std::pair<MapChunkVecConstIter, MapChunkVecConstIter> MapChunkVecConstIterPair;

  class MapChunkDB {
    
  public:

    MapChunkDB(const MapWrapperPVec& maps, size_t frags_per_chunk);

    void sort_chunks();

    MapChunkVecConstIterPair query(int lb, int ub) const;
    MapChunkVecConstIterPair query(const IntPair&) const;

    int count_compatible_seeds(const IntPairVec& bounds) const;
    
    AlignmentVector get_compatible_alignments(const IntPairVec& bounds,
      size_t max_unmatched = std::numeric_limits<size_t>::max()) const;

    RefAlignmentVec get_compatible_alignments_best(const IntPairVec& bounds,
      size_t max_unmatched = std::numeric_limits<size_t>::max()) const;

    template <typename BoundIter >
    int count_compatible_right_bfs(const MapChunk * p_start_chunk, BoundIter bound_start, BoundIter bound_end) const;

    template <typename BoundIter >
    int count_compatible_left_bfs(const MapChunk * p_start_chunk, BoundIter bound_start, BoundIter bound_end) const;
    
    template <typename BoundIter >
    AlignmentVector get_compatible_right_bfs(const MapChunk * p_start_chunk, BoundIter bound_start, BoundIter bound_end) const;

    template <typename BoundIter >
    AlignmentVector get_compatible_left_bfs(const MapChunk * p_start_chunk,
      BoundIter bound_start, BoundIter bound_end) const;

    template <typename BoundIter >
    AlignmentVector get_compatible_right_dfs(const MapChunk * p_start_chunk,
      BoundIter bound_start, BoundIter bound_end,
      size_t max_unmatched = std::numeric_limits<size_t>::max()) const;
    
    template <typename BoundIter >
    AlignmentVector get_compatible_left_dfs(const MapChunk * p_start_chunk,
      BoundIter bound_start, BoundIter bound_end, 
      size_t max_unmatched = std::numeric_limits<size_t>::max()) const;
      

    // This is a cache friendly layout of the MapChunks.
    // MapChunks that belong to the same map with the same start index
    // are adjacent.
    std::vector<MapChunk> map_chunks_;

    // This is a query friendly layout of MapChunks, sorted by size.
    std::vector<const MapChunk *> map_chunk_index_;

    std::unordered_map<const Map *, ChunksAtIndex> map_to_chunks_at_start_;
    std::unordered_map<const Map *, ChunksAtIndex> map_to_chunks_at_end_;
    size_t frags_per_chunk_;

  };

  ////////////////////////////////////
  // Search node for DFS
  template< typename BoundIter>
  class SearchNode {
  public:

      SearchNode(const MapChunk * chunk,
        MapChunkVecConstIter chunk_next,
        MapChunkVecConstIter chunk_next_end,
        BoundIter b,
        int num_unmatched) :
          chunk_(chunk),
          chunk_next_(chunk_next),
          chunk_next_end_(chunk_next_end),
          bound_iter_next_(b),
          num_unmatched_(num_unmatched)
        { };
      
      const MapChunk * chunk_; // The current chunk we are at
      MapChunkVecConstIter chunk_next_; // the next chunk to consider
      const MapChunkVecConstIter chunk_next_end_; // the end of the next chunk vector (so we know when we are done trying to extend)
      BoundIter bound_iter_next_;  // The next bound to consider
      size_t num_unmatched_; // The number of misses incurred so far, not counting current chunk.
  };


  template <typename BoundIter >
  int MapChunkDB::count_compatible_right_bfs(const MapChunk * p_start_chunk, BoundIter bound_start, BoundIter bound_end) const {
    
    using std::vector;

    // Get all chunks which follow this one and test which are compatible with the rest of the bounds
    if(bound_start >= bound_end) {
      return 1;
    }

    const size_t num_bounds = bound_end - bound_start;
    
    MapChunkConstPVec active_set;
    active_set.reserve(num_bounds);
    active_set.push_back(p_start_chunk);

    const Map* p_map = p_start_chunk->get_map();
    const ChunksAtIndex& chunks_at_start = map_to_chunks_at_start_.at(p_map);

    for(BoundIter bi = bound_start; bi != bound_end; bi++) {
      
      MapChunkConstPVec new_active_set;

      auto lb = bi->first;
      auto ub = bi->second;

      for(MapChunkVecConstIter ci = active_set.begin(); ci != active_set.end(); ci++) {

        const MapChunk* p_chunk = *ci;

        // Get chunks that follow this one from the database.
        const Map* p_map = p_chunk->get_map();
        if(p_chunk->end_ == p_map->frags_.size()) {

          // We've reached the end of the map!
          continue;

        }

        const MapChunkConstPVec& next_chunks = chunks_at_start[p_chunk->end_];
        for(size_t i = 0; i < next_chunks.size(); i++) {
          const MapChunk * p_next = next_chunks[i];
          if(p_next->size_ >= lb && p_next->size_ <= ub) {
            new_active_set.push_back(p_next);
          }
        }

      }

      if (new_active_set.empty()) { return 0; }
      active_set = std::move(new_active_set);

    }

    return active_set.size();

  }

  template <typename BoundIter >
  int MapChunkDB::count_compatible_left_bfs(const MapChunk * p_start_chunk, BoundIter bound_start, BoundIter bound_end) const {
    
    using std::vector;

    // Get all chunks which follow this one and test which are compatible with the rest of the bounds
    if(bound_start >= bound_end) {
      return 1;
    }

    const size_t num_bounds = bound_end - bound_start;
    
    MapChunkConstPVec active_set;
    active_set.reserve(num_bounds);
    active_set.push_back(p_start_chunk);

    const Map* p_map = p_start_chunk->get_map();
    const ChunksAtIndex& chunks_at_end = map_to_chunks_at_end_.at(p_map);

    for(BoundIter bi = bound_start; bi != bound_end; bi++) {
        
      vector<const MapChunk *> new_active_set;

      auto lb = bi->first;
      auto ub = bi->second;

      for(MapChunkVecConstIter ci = active_set.begin(); ci != active_set.end(); ci++) {

        // Get chunks that follow this one from the database.
        const MapChunk* p_chunk = *ci;

        if(p_chunk->start_ == 0) {
          // We've reached the start of the map!
          continue;
        }

        const MapChunkConstPVec& next_chunks = chunks_at_end[p_chunk->start_]; // Get chunks that end where the current chunk starts.
        for(size_t i = 0; i < next_chunks.size(); i++) {
          const MapChunk * p_next = next_chunks[i];
          if(p_next->size_ >= lb && p_next->size_ <= ub) {
            new_active_set.push_back(p_next);
          }
        }

      }

      if (new_active_set.empty()) { return 0; }
      active_set = std::move(new_active_set);

    }

    return active_set.size();

  }

  //////////////////////////////////////////////////////////////////////////////////////
  // This uses BFS to get the compatible alignments to the right.
  // This is not memory efficient, since it maintains all of the active alignments.
  template <typename BoundIter >
  AlignmentVector MapChunkDB::get_compatible_right_bfs(const MapChunk * p_start_chunk, BoundIter bound_start, BoundIter bound_end) const {
    
    using std::vector;

    // Get all chunks which follow this one and test which are compatible with the rest of the bounds
    if(bound_start >= bound_end) {
      return AlignmentVector();
    }

    const size_t num_bounds = bound_end - bound_start;
    
    AlignmentVector active_set;
    active_set.reserve(num_bounds);

    // Write the current alignment to the active set
    Alignment cur_alignment;
    cur_alignment.push_back(p_start_chunk);
    active_set.push_back(cur_alignment);

    const Map* p_map = p_start_chunk->get_map();
    const ChunksAtIndex& chunks_at_start = map_to_chunks_at_start_.at(p_map);

    for(BoundIter bi = bound_start; bi != bound_end; bi++) {
        
      AlignmentVector new_active_set;

      auto lb = bi->first;
      auto ub = bi->second;

      for(auto ci = active_set.begin(); ci != active_set.end(); ci++) {

        // Get chunks that follow this one from the database.
        const Alignment& aln = *ci;
        size_t aln_end = aln.back()->end_;

        if(aln_end == p_map->frags_.size()) {
          // We've reached the end of the map!
          continue;
        }

        const MapChunkConstPVec&  next_chunks = chunks_at_start[aln_end]; // Get chunks that start where current alignment ends
        for(size_t i = 0; i < next_chunks.size(); i++) {

          const MapChunk * p_next = next_chunks[i];

          if(p_next->size_ >= lb && p_next->size_ <= ub) {
            Alignment new_aln(aln); // unnecessary copy if this doesnt get anywhere
            new_aln.push_back(p_next);
            new_active_set.push_back(std::move(new_aln));
          }

        }
      }

      if (new_active_set.empty()) { return AlignmentVector(); }

      active_set = std::move(new_active_set);

    }

    return active_set;

  }


  //////////////////////////////////////////////////////////////////////////////////////
  // This uses BFS to get the compatible alignments to the left.
  // This is not memory efficient, since it maintains all of the active alignments.
  template <typename BoundIter >
  AlignmentVector MapChunkDB::get_compatible_left_bfs(const MapChunk * p_start_chunk, BoundIter bound_start, BoundIter bound_end) const {
    
    using std::vector;

    // Get all chunks which follow this one and test which are compatible with the rest of the bounds
    if(bound_start >= bound_end) {
      return AlignmentVector();
    }

    if(p_start_chunk->start_ == 0) {
      //Cannot extend left!
      return AlignmentVector();
    }

    const size_t num_bounds = bound_end - bound_start;
    
    AlignmentVector active_set;
    active_set.reserve(num_bounds);


    /////////////////////////////////////////////////////////////
    // Find the first chunk which extends left.
    const Map* p_map = p_start_chunk->get_map();
    const ChunksAtIndex& chunks_at_end = map_to_chunks_at_end_.at(p_map);
    BoundIter bi = bound_start;
    auto lb = bi->first;
    auto ub = bi->second;
    const MapChunkConstPVec& next_chunks = chunks_at_end[p_start_chunk->start_]; // Get chunks that end where the current chunk starts.
    for(size_t i = 0; i < next_chunks.size(); i++) {

      const MapChunk * p_next = next_chunks[i];

      if(p_next->size_ >= lb && p_next->size_ <= ub) {
        Alignment new_aln = { p_next };
        new_aln.push_back(p_next);
        active_set.push_back(std::move(new_aln));
      }
    }


    for(BoundIter bi = bound_start; bi != bound_end; bi++) {
        
      AlignmentVector new_active_set;

      auto lb = bi->first;
      auto ub = bi->second;

      for(auto ci = active_set.begin(); ci != active_set.end(); ci++) {

        // Get chunks that follow this one from the database.
        const Alignment& aln = *ci;
        size_t aln_start = aln.back()->start_;
        if (aln_start == 0) {
          // We've reached the start of the map!
          continue;
        }

        const MapChunkConstPVec& next_chunks = chunks_at_end[aln_start]; // Get chunks that end where the current chunk starts.
        for(size_t i = 0; i < next_chunks.size(); i++) {

          const MapChunk * p_next = next_chunks[i];

          if(p_next->size_ >= lb && p_next->size_ <= ub) {
            Alignment new_aln(aln);
            new_aln.push_back(p_next);
            new_active_set.push_back(std::move(new_aln));
          }

        }
      }

      if (new_active_set.empty()) { return AlignmentVector(); }

      active_set = std::move(new_active_set);

    }

    return active_set;

  }

  template <typename BoundIter >
  AlignmentVector MapChunkDB::get_compatible_right_dfs(
    const MapChunk * p_start_chunk,
    BoundIter bound_start,
    BoundIter bound_end,
    size_t max_unmatched) const {
    // Search from the start_chunk in the forward direction (with respect to the map containing the chunk)
    // The start chunk is the "middle" chunk of the alignment. We look to it's right.
    // Return a vector of alignments (i.e. a vector of vector<MapChunk *> ) where each vector<MapChunk*> is compatible with the bounds
    // from bound_start to bound_end.
    // The Vector<MapChunk *> is oriented the same way as the bounds given by bound_start and bound_end.
    
    using std::vector;

    // std::cerr << "get_compatible_right: " << p_start_chunk.get_map()->name_ << " " << *p_start_chunk << std::endl;

    // Get all chunks which follow this one and test which are compatible with the rest of the bounds
    if(bound_start >= bound_end) {
      return AlignmentVector();
    }

    const Map* p_map = p_start_chunk->get_map();
    const ChunksAtIndex& chunks_at_start = map_to_chunks_at_start_.at(p_map);

    if(p_start_chunk->end_ >= p_map->frags_.size()) {
      // Can't extend off the right end of the map.
      return AlignmentVector();
    }

    AlignmentVector alignments;

    typedef SearchNode<BoundIter> Node;
    std::vector<Node> stack;
    stack.reserve(bound_end - bound_start);

    stack.emplace_back(p_start_chunk,
         chunks_at_start[p_start_chunk->end_].begin(),
         chunks_at_start[p_start_chunk->end_].end(),
         bound_start, 0); 

    Alignment cur_aln; // Don't include the start chunk in the alignment.

    while(!stack.empty()) {

      process_stack:

        Node& cur = stack.back();

        // If there is no next bound, we've reached the end of the search bounds and our currently alignment is compatible.
        // Save this alignment and backtrack.
        if(cur.bound_iter_next_ == bound_end) {
          alignments.push_back(cur_aln);
          stack.pop_back();
          if(!cur_aln.empty()) { cur_aln.pop_back(); }
          // std::cerr << "stack pop!" << std::endl;
          continue;
        }

        // std::cerr << "Stack size: " << stack.size()
        //       << " cur node: " << cur.chunk_.get_map()->name_ << " " << *cur.chunk_
        //       << " lb: " << cur.bound_iter_next_->first
        //       << " ub: " << cur.bound_iter_next_->second
        //       << std::endl;

        const int lb = cur.bound_iter_next_->first;
        const int ub = cur.bound_iter_next_->second;

        // Get the next successor node, and if it is compatible with the bound, put it on the stack.
        bool pushed = false;
        while(cur.chunk_next_ != cur.chunk_next_end_) {
          const MapChunk* cur_chunk = cur.chunk_;
          const MapChunk* next_chunk = *cur.chunk_next_++;
          if (next_chunk->end_ >= p_map->frags_.size()) continue; // Don't extend off of the map.
          if (cur.num_unmatched_ >= max_unmatched) continue;
          if (next_chunk->size_ >= lb && next_chunk->size_ <= ub) {
            
            stack.emplace_back(next_chunk,
                chunks_at_start[next_chunk->end_].begin(),
                chunks_at_start[next_chunk->end_].end(),
                cur.bound_iter_next_ + 1,
                cur.num_unmatched_ + cur_chunk->num_unmatched());

            cur_aln.push_back(next_chunk);
            // std::cerr << "stack push!" << std::endl;
            pushed = true;
            goto process_stack; // Wow, I think this is a perfectly good use of goto!
            break;
          }
          
        }

      // if (pushed) continue;
      assert(!pushed);
      // If we did not find a successor, pop this node.
      if (!pushed) {
        stack.pop_back();
        if(!cur_aln.empty()) { cur_aln.pop_back(); } 
        // std::cerr << "stack pop!" << std::endl;
      } 

    }

    return alignments;

  }

  template <typename BoundIter >
  AlignmentVector MapChunkDB::get_compatible_left_dfs(const MapChunk * p_start_chunk,
    BoundIter bound_start,
    BoundIter bound_end,
    size_t max_unmatched) const {

    // Search from the start_chunk in the reverse direction (with respect to the map containing the chunk)
    // Return a vector of alignments (i.e. a vector of vector<MapChunk *> ) where each vector<MapChunk*> is compatible with the bounds
    // from bound_start to bound_end.
    // The Vector<MapChunk *> is oriented the same way as the bounds given by bound_start and bound_end.

    using std::vector;

    // std::cerr << "get_compatible_left2: " << p_start_chunk.get_map()->name_ << " " << *p_start_chunk << std::endl;

    // Get all chunks which follow this one and test which are compatible with the rest of the bounds
    if(bound_start >= bound_end) {
      return AlignmentVector();
    }

    const Map* p_map = p_start_chunk->get_map();
    const ChunksAtIndex& chunks_at_end = map_to_chunks_at_end_.at(p_map);

    if(p_start_chunk->start_ == 0) {
      // Can't extend off the left end of the map.
      return AlignmentVector();
    }

    AlignmentVector alignments;

    typedef SearchNode<BoundIter> Node;
    std::vector<Node> stack;
    stack.reserve(bound_end - bound_start);
    stack.emplace_back(p_start_chunk,
        chunks_at_end[p_start_chunk->start_].begin(),
        chunks_at_end[p_start_chunk->start_].end(),
        bound_start, 0); 

    Alignment cur_aln;

    while(!stack.empty()) {
      process_stack:
        Node& cur = stack.back();

        // If there is no next bound, save this alignment and backtrack.
        if(cur.bound_iter_next_ == bound_end) {
          alignments.push_back(cur_aln);
          stack.pop_back();
          if(!cur_aln.empty()) { cur_aln.pop_back(); }
          // std::cerr << "stack pop!" << std::endl;
          continue;
        }

        const int lb = cur.bound_iter_next_->first;
        const int ub = cur.bound_iter_next_->second;

        // Get the next successor node, and if it is compatible with the bound, put it on the stack.
        bool pushed = false;
        while(cur.chunk_next_ != cur.chunk_next_end_) {
          const MapChunk* cur_chunk = cur.chunk_;
          const MapChunk* next_chunk = *cur.chunk_next_++;

          if (next_chunk->start_ == 0 ) continue; // Don't extend off of the reference map.
          if (cur.num_unmatched_ >= max_unmatched) continue;
          if (next_chunk->size_ >= lb && next_chunk->size_ <= ub) {
            
            stack.emplace_back(next_chunk,
                chunks_at_end[next_chunk->start_].begin(),
                chunks_at_end[next_chunk->start_].end(),
                cur.bound_iter_next_ + 1,
                cur.num_unmatched_ + cur_chunk->num_unmatched());
            cur_aln.push_back(next_chunk);
            pushed = true;
            // std::cerr << "stack push!" << std::endl;
            goto process_stack;
            break;
          }
        }

      // if (pushed) continue;
      assert(!pushed);

      // If we did not find a successor, pop this node.
      if (!pushed) {
        stack.pop_back();
        if(!cur_aln.empty()) { cur_aln.pop_back(); }
        // std::cerr << "stack pop!" << std::endl;
      } 

    }

    return alignments;

  }

  std::ostream& operator<<(std::ostream& os, const Alignment& );

}

#endif

