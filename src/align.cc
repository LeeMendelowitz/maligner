#include <vector>
#include <cassert>
#include <iostream>
#include <utility>
#include <algorithm>

using std::cerr;

#include "utils.h"
#include "types.h"
#include "align.h"
#include "globals.h"
#include "ScoreMatrix.h"
#include "ScoreCell.h"

using Constants::INF;

#define DEBUG 0
#define GET_BEST_DEBUG 0
#define BUILD_TRAIL_DEBUG 0
#define FILL_DEBUG 0
#define BREAKS_DEBUG 0
#define RESCALE_DEBUG 0 

inline double sizing_penalty(int query_size, int ref_size, const AlignOpts& align_opts) {

  /* TODO: This can be baked into the dynamic programming routine */
  /* For each reference fragment, compute the standard deviation ahead of time */
  
  double delta = query_size - ref_size;
  double sd = align_opts.sd_rate * ref_size;

  if (sd < align_opts.min_sd) {
    sd = align_opts.min_sd;
  }

  double sd_1 = 1.0 / sd;
  double penalty = delta * sd_1;
  penalty = penalty * penalty;
  return penalty;

}

// rescale the query chunks using the query_scaling_factor, and
// recompute the sizing error for those chunks.
void Alignment::rescale_matched_chunks(AlignOpts& align_opts) {
  rescaled_score = score;
  rescaled_score.sizing_score = 0.0;
  const size_t l = rescaled_matched_chunks.size();
  for (size_t i = 0; i < l; i++) {
    MatchedChunk& mc = rescaled_matched_chunks[i];
    int old_query_size = mc.query_chunk.size;
    int new_query_size = query_scaling_factor*old_query_size;
    double old_sizing_score = mc.score.sizing_score;
    double new_sizing_score = sizing_penalty(new_query_size, mc.ref_chunk.size, align_opts);
    mc.query_chunk.size = new_query_size;
    if (!mc.query_chunk.is_boundary && !mc.ref_chunk.is_boundary) {
      mc.score.sizing_score = new_sizing_score;
      rescaled_score.sizing_score += new_sizing_score;
      #if RESCALE_DEBUG > 0
      std::cerr << "old_q: " << old_query_size << " new_q: " << new_query_size
              << " old_sizing_score: " << old_sizing_score << " new: " << new_sizing_score << "\n";
      #endif
    }
  }
}


/*
  Populate a score matrix using dynamic programming for ungapped alignment.

  The score matrix should have the same number of columns as the reference.

  The ScoreMatrix should already have the same nubmer of columns as the reference,
  and should have enough rows to accomodate the query.
*/
void fill_score_matrix(const AlignTask& align_task) {

  // Unpack the alignment task
  const IntVec& query = *align_task.query;
  const IntVec& ref = *align_task.ref;
  ScoreMatrix& mat = *align_task.mat;
  AlignOpts& align_opts = *align_task.align_opts;

  mat.reset();

  const int m = query.size() + 1;
  const int n = ref.size() + 1;

  // Note: Number of rows may be different from m if matrix is padded with extra rows.
  const int num_rows = mat.getNumRows();

  assert((int) mat.getNumCols() == n);
  assert((int) mat.getNumRows() >= m);

  #if FILL_DEBUG > 0
  cerr << "m: " << m
       << " n: " << n
       << " num_rows: " << num_rows
       << " num_cols: " << mat.getNumCols()
       << "\n";
  #endif

  // Initialize the first row
  for (int j = 0; j < n; j++) {
    ScoreCell* pCell = mat.getCell(0,j);
    pCell->score_ = 0.0;
    pCell->backPointer_ = nullptr;
  }

  // Initialize the first column
  for (int i = 1; i < m; i++ ) {
    ScoreCell* pCell = mat.getCell(i,0);
    pCell->score_ = -INF;
    pCell->backPointer_ = nullptr;
  }

  // Initialize the body of the matrix.
  for (int j = 1; j < n; j++) {
    
    // Matrix is column major ordered.
    int offset = j*num_rows;

    for (int i = 1; i < m; i++) {
      ScoreCell* pCell = mat.getCell(offset + i);
      pCell->score_ = -INF;
      pCell->backPointer_ = nullptr;
    }

  }

  #if BREAKS_DEBUG > 0
  int num_breaks = 0;
  #endif


  for (int j = 1; j < n; j++) {
    
    int l0 = (j > align_opts.ref_max_misses + 1) ? j - align_opts.ref_max_misses - 1 : 0;
    
    const int offset = num_rows*j;
    
    for (int i = 1; i < m; i++) {

      ScoreCell* pCell = mat.getCell(offset + i);

      // Try all allowable extensions

      ScoreCell* backPointer = nullptr;
      double best_score = -INF;
      int k0 = (i > align_opts.query_max_misses) ? i - align_opts.query_max_misses - 1 : 0;
      int ref_size = 0;
      for(int l = j-1; l >= l0; l--) {

        const bool is_ref_boundary = l == 0 || j == n - 1;
        
        ref_size += ref[l];

        int ref_miss = j - l - 1; // sites in reference unaligned to query
        double ref_miss_score = ref_miss * align_opts.ref_miss_penalty;

        int query_size = 0;
        const int offset_back = num_rows*l;

        for(int k = i-1; k >= k0; k--) {

          const bool is_query_boundary =  k == 0 || k == m - 1;

          #if FILL_DEBUG > 0
          cerr << "i: " << i
               << " j: " << j
               << " k: " << k
               << " l: " << l
               << "\n";
          #endif


          query_size += query[k];
          ScoreCell* pTarget = mat.getCell(offset_back + k);
          if (pTarget->score_ == -INF) continue;

          int query_miss = i - k - 1; // sites in query unaligned to reference
          double query_miss_score = query_miss * align_opts.query_miss_penalty;

          // Add sizing penalty only if this is not a boundary fragment.
          double size_penalty = 0.0;
          if (!is_ref_boundary && !is_query_boundary) {
            size_penalty = sizing_penalty(query_size, ref_size, align_opts);
          }

          if (query_size > ref_size && size_penalty > align_opts.max_chunk_sizing_error) {
            // The query chunk is already too big for the reference chunk.
              #if BREAKS_DEBUG > 0
                num_breaks++;
              #endif
            break;
          }

          // If the sizing penalty is too large, continue and do not populate matrix.
          if (size_penalty > align_opts.max_chunk_sizing_error) {
            continue;
          }

          double chunk_score = -size_penalty - query_miss_score - ref_miss_score;

          if (chunk_score + pTarget->score_ > best_score) {
            backPointer = pTarget;
            best_score = chunk_score + pTarget->score_;
          }

        } // for int k
      } // for int l

      // Assign the backpointer and score to pCell
      if (backPointer) {
        pCell->backPointer_ = backPointer;
        pCell->score_ = best_score;
      }

    } // for int i
  } // for int j
  
  #if BREAKS_DEBUG > 0
    std::cout << "num breaks: " << num_breaks << "\n";
  #endif

} // fill_score_matrix

void fill_score_matrix_using_partials(const AlignTask& align_task) {
  /*
  Fill score matrix using partial sums
  */

  // Unpack the alignment task
  const IntVec& query = *align_task.query;
  const IntVec& ref = *align_task.ref;
  const PartialSums& query_partial_sums = *align_task.query_partial_sums;
  const PartialSums& ref_partial_sums = *align_task.ref_partial_sums;

  ScoreMatrix& mat = *align_task.mat;
  mat.reset();
  AlignOpts& align_opts = *align_task.align_opts;

  // Compute the miss penalties
  IntVec query_miss_penalties(align_opts.query_max_misses+1, align_opts.query_miss_penalty);
  for (int i = 0; i < align_opts.query_max_misses+1; i++) {
    query_miss_penalties[i] *= (double) i;
  }

  IntVec ref_miss_penalties(align_opts.ref_max_misses+1, align_opts.ref_miss_penalty);
  for (int i = 0; i < align_opts.ref_max_misses+1; i++) {
    ref_miss_penalties[i] *= (double) i;
  }

  const int m = query.size() + 1;
  const int n = ref.size() + 1;

  // Note: Number of rows may be different from m if matrix is padded with extra rows.
  const int num_rows = mat.getNumRows();

  assert((int) mat.getNumCols() == n);
  assert((int) mat.getNumRows() >= m);

  #if FILL_DEBUG > 0
  cerr << "m: " << m
       << " n: " << n
       << " num_rows: " << num_rows
       << " num_cols: " << mat.getNumCols()
       << "\n";
  #endif

  // Initialize the first row
  for (int j = 0; j < n; j++) {
    ScoreCell* pCell = mat.getCell(0,j);
    pCell->score_ = 0.0;
    pCell->backPointer_ = nullptr;
  }

  // Initialize the first column
  for (int i = 1; i < m; i++ ) {
    ScoreCell* pCell = mat.getCell(i,0);
    pCell->score_ = -INF;
    pCell->backPointer_ = nullptr;
  }

  // Initialize the body of the matrix.
  for (int j = 1; j < n; j++) {
    
    // Matrix is column major ordered.
    int offset = j*num_rows;

    for (int i = 1; i < m; i++) {
      ScoreCell* pCell = mat.getCell(offset + i);
      pCell->score_ = -INF;
      pCell->backPointer_ = nullptr;
    }

  }

  #if BREAKS_DEBUG > 0
  int num_breaks = 0;
  #endif

  
  for (int j = 1; j < n; j++) {

    const IntVec& ref_ps = ref_partial_sums[j-1]; // reference partial sum

    int l0 = (j > align_opts.ref_max_misses + 1) ? j - align_opts.ref_max_misses - 1 : 0;
    const int offset = num_rows*j;
    
    for (int i = 1; i < m; i++) {

      const IntVec& query_ps = query_partial_sums[i-1]; // query partial sum

      ScoreCell* pCell = mat.getCell(offset + i);

      // Try all allowable extensions

      ScoreCell* backPointer = nullptr;
      double best_score = -INF;
      int k0 = (i > align_opts.query_max_misses) ? i - align_opts.query_max_misses - 1 : 0;

      for(int l = j-1; l >= l0; l--) {

        const bool is_ref_boundary = l == 0 || j == n - 1;

        int ref_miss = j - l - 1; // sites in reference unaligned to query
        double ref_miss_score = ref_miss_penalties[ref_miss];
        int ref_size = ref_ps[ref_miss];

        const int offset_back = num_rows*l;

        for(int k = i-1; k >= k0; k--) {

          const bool is_query_boundary =  (k == 0 || k == m - 1) && !align_opts.query_is_bounded;

          #if FILL_DEBUG > 0
          cerr << "i: " << i
               << " j: " << j
               << " k: " << k
               << " l: " << l
               << "\n";
          #endif

          ScoreCell* pTarget = mat.getCell(offset_back + k);
          if (pTarget->score_ == -INF) continue;

          int query_miss = i - k - 1; // sites in query unaligned to reference
          double query_miss_score = query_miss_penalties[query_miss];
          //double query_miss_score = query_miss * align_opts.query_miss_penalty;
          int query_size = query_ps[query_miss];

          // Add sizing penalty only if this is not a boundary fragment.
          double size_penalty = 0.0;
          if (!is_ref_boundary && !is_query_boundary) {

            size_penalty = sizing_penalty(query_size, ref_size, align_opts);
      
          }

          // Break if the query chunk is already too big for the reference
          if (query_size > ref_size && size_penalty > align_opts.max_chunk_sizing_error) {
              #if BREAKS_DEBUG > 0
                num_breaks++;
              #endif
            break;
          }

          // If the sizing penalty is too large, continue and do not populate matrix.
          if (size_penalty > align_opts.max_chunk_sizing_error) {
            continue;
          }

          double chunk_score = -size_penalty - query_miss_score - ref_miss_score;

          if (chunk_score + pTarget->score_ > best_score) {
            backPointer = pTarget;
            best_score = chunk_score + pTarget->score_;
          }

        } // for int k
      } // for int l

      // Assign the backpointer and score to pCell
      if (backPointer) {
        pCell->backPointer_ = backPointer;
        pCell->score_ = best_score;
      }

    } // for int i
  } // for int j
  
  #if BREAKS_DEBUG > 0
    std::cout << "num breaks: " << num_breaks << "\n";
  #endif

} // fill_score_matrix


bool get_best_alignment(const AlignTask& task, ScoreCellPVec& trail) {

  // Go to the last row of the ScoreMatrix and identify the best score.
  IntVec& query = *task.query;
  IntVec& ref = *task.ref;
  ScoreMatrix& mat = *task.mat;

  const int m = query.size() + 1;
  const int n = ref.size() + 1;
  const int num_rows = mat.getNumRows();
  const int num_cols = mat.getNumCols();
  const int last_row = m - 1;

  double best_score = -INF;
  ScoreCell * p_best_cell = nullptr;
  int index = last_row;

  // Get the cell with the best score in the last row.
  for (int i = 0; i < n; i++, index += num_rows) {
    
    ScoreCell * pCell = mat.getCell(index);

    #if GET_BEST_DEBUG > 0
    cerr << "index: " << index << ", ";
    cerr << "cell: " << *pCell << "\n";
    #endif

    if (pCell->score_ > best_score) {
      p_best_cell = pCell;
      best_score = pCell->score_;
    }
  }

  #if GET_BEST_DEBUG > 0
    if (p_best_cell) {
        cerr << "\np_best_cell: " << *p_best_cell <<  "\n";
    } else {
        cerr << "\np_best_cell: " << p_best_cell <<  "\n";
    }
  #endif
 
  if (p_best_cell && p_best_cell->score_ > -INF) {
    // Get the traceback
    trail.clear();
    trail.reserve(m);
    build_trail(p_best_cell, trail);
    return true;
  }

  return false;

}

/////////////////////////////////////////////////////////////////////
// Get the best n alignments in the task.
// The score matrix should already be filled out.
//
// Append trail seeds to the trail_seeds in the alignment task.
// Return True if alignments were found.
int get_best_alignments(const AlignTask& task) {

  // Go to the last row of the ScoreMatrix and identify the best score.
  AlignOpts& align_opts = *task.align_opts;
  IntVec& query = *task.query;
  IntVec& ref = *task.ref;
  ScoreMatrix& mat = *task.mat;
  AlignmentPVec& alignments = *task.alignments;


  const int m = query.size() + 1;
  const int n = ref.size() + 1;
  const int num_rows = mat.getNumRows();
  const int num_cols = mat.getNumCols();
  const int last_row = m - 1;

  int num_alignments_found = 0;
  alignments.reserve(align_opts.alignments_per_reference);
  for (int r = 0; r < align_opts.alignments_per_reference; r++) {

    double best_score = -INF;
    int best_cell_col = -1;
    ScoreCell * p_best_cell = nullptr;

    // Get the cell with the best score in the last row that is
    // still in play.
    int index = last_row;
    for (int i = 0; i < n; i++, index += num_rows) {
      
      if (!mat.cell_in_play(i)) continue;

      ScoreCell * pCell = mat.getCell(index);

      #if GET_BEST_DEBUG > 0
      cerr << "index: " << index << ", ";
      cerr << "cell: " << pCell << "\n";
      #endif

      if (pCell->backPointer_ && pCell->score_ > best_score ) {
        p_best_cell = pCell;
        best_score = pCell->score_;
        best_cell_col = i;
      } 
    }

    #if GET_BEST_DEBUG > 0
      if (p_best_cell) {
          cerr << "\np_best_cell: " << *p_best_cell <<  "\n";
      } else {
          cerr << "\np_best_cell: " << p_best_cell <<  "\n";
      }
    #endif
   
    bool have_alignment = p_best_cell && p_best_cell->score_ > -INF && p_best_cell->backPointer_;
    if (!have_alignment) break;

    // Build the alignment
    alignments.push_back(alignment_from_cell(task, p_best_cell));

    num_alignments_found++;

    // Mark neighboring cells as out of play.
    int lb = best_cell_col - align_opts.min_alignment_spacing + 1;
    int ub = best_cell_col + align_opts.min_alignment_spacing;
    if (lb < 0) lb = 0;
    if (ub > n) ub = n;
    for (int col = lb; col < ub; col++) {
      mat.mark_cell_in_play(col, false);
    }
  }

  return num_alignments_found;
}


void build_trail(ScoreCell* pCell, ScoreCellPVec& trail) {

  #if BUILD_TRAIL_DEBUG > 0
  std::cerr << "Building trail from " << pCell << std::endl;
  #endif

  ScoreCell* pCur = pCell;
  while (pCur != nullptr) {
    #if BUILD_TRAIL_DEBUG > 0
    std::cerr << "pushing back: " << pCur << std::endl;
    #endif
    trail.push_back(pCur);
    pCur = pCur->backPointer_;
  }
}


// trail: starts from end of alignment
void build_chunk_trail(const AlignTask& task, ScoreCellPVec& trail, ChunkVec& query_chunks, ChunkVec& ref_chunks) {

  const int ts = trail.size();
  if (ts == 0) { return; }

  const IntVec& query = *task.query;
  const IntVec& ref = *task.ref;
  const AlignOpts& align_opts = *task.align_opts;

  query_chunks.clear();
  query_chunks.reserve(ts-1);
  ref_chunks.clear();
  ref_chunks.reserve(ts-1);
  //assert(query_chunks.size() == ts - 1);
  //assert(ref_chunks.size() == ts - 1);

  ScoreCell* pLast = trail[0];
  int ml = pLast->q_;

  #if DEBUG > 0
  std::cerr << "build_chunk_trail for: " << trail << std::endl;
  #endif

  int nl = pLast->r_;

  auto qi(query_chunks.begin());
  auto ri(ref_chunks.begin());
  for(int i = 1; i < ts; i++, qi++, ri++) {
    
    const ScoreCell* pCell = trail[i];
    int m = pCell->q_; // index of query site, one based
    int n = pCell->r_; // index of ref size, one based

    #if DEBUG > 0
    cerr << "cell: " << pCell << " " << *pCell << "\n"
         << " i: " << i 
         << " m: " << m << " ml: " << ml
         << " n: " << n << " nl: " << nl
         << "\n";
    #endif

    assert(m < ml);
    assert(n < nl);

    bool is_query_boundary = (m == 0 || ml == query.size()) && !align_opts.query_is_bounded;
    bool is_ref_boundary = (n == 0) || (nl == ref.size());

    // Build chunks
    int q_size = sum(query, m, ml);
    int r_size = sum(ref, n, nl);
    /* USE C++11 instead
    Chunk q_chunk(m, ml, q_size, is_query_boundary); 
    Chunk r_chunk(n, nl, r_size, is_ref_boundary);
    query_chunks.push_back(q_chunk);
    ref_chunks.push_back(r_chunk);
    */
    //query_chunks.emplace(qi, m, ml, q_size, is_query_boundary);
    //ref_chunks.emplace(ri, n, nl, r_size, is_ref_boundary);
    query_chunks.emplace_back(m, ml, q_size, is_query_boundary);
    ref_chunks.emplace_back(n, nl, r_size, is_ref_boundary);

    ml = m;
    nl = n;
  }
}

void print_chunk_trail(const ChunkVec& query_chunks, const ChunkVec& ref_chunks) {

  assert(query_chunks.size() == ref_chunks.size());
  int cs = query_chunks.size();
  cerr << "\n\ntrail: ";
  for(int i = 0; i < cs; i++) {

    const Chunk& qc = query_chunks[i];
    const Chunk& rc = ref_chunks[i];

    #if DEBUG > 0
    cerr << "q: " << qc << " r: " << rc << "\n";
    #endif

  }

}

// Make and return an alignment from the trail through the
// score matrix.
Alignment * alignment_from_trail(const AlignTask& task, ScoreCellPVec& trail) {

    const AlignOpts& align_opts = *task.align_opts;
    const IntVec& query = *task.query;

    const size_t chunk_vec_size = trail.size() - 1;

    ChunkVec query_chunks, ref_chunks;
    MatchedChunkVec matched_chunks;

    build_chunk_trail(task, trail, query_chunks, ref_chunks);

    // The trail and the chunks are oriented from the end of the alignment.
    // Build MatchChunk's in the forward direction.
    assert(query_chunks.size() == ref_chunks.size());

    const size_t n = query_chunks.size();
    const int num_query_frags = query.size();

    Score total_score(0.0, 0.0, 0.0);
    matched_chunks.reserve(n);
    for (int i = n-1; i >= 0; i--) {
      
        Chunk& qc = query_chunks[i];
        Chunk& rc = ref_chunks[i];
        int ref_misses = rc.end - rc.start - 1; // sites in reference that are unaligned to query
        int query_misses = qc.end - qc.start - 1; // sites in query that are unaligned to reference

        double query_miss_score = align_opts.query_miss_penalty * query_misses;
        double ref_miss_score = align_opts.ref_miss_penalty * ref_misses;
        double sizing_score = 0.0;

        const bool query_is_boundary = (qc.start == 0 || qc.end == num_query_frags) && !align_opts.query_is_bounded;

        if (!query_is_boundary) {
          sizing_score = sizing_penalty(qc.size, rc.size, align_opts);
        }

        Score score(query_miss_score, ref_miss_score, sizing_score);

        total_score.query_miss_score += query_miss_score;
        total_score.ref_miss_score += ref_miss_score;
        total_score.sizing_score += sizing_score;

        // MatchedChunk mc(query_chunks[i], ref_chunks[i], score);
        // Try fancy C++11:
        matched_chunks.emplace_back(query_chunks[i], ref_chunks[i], score);
        
    }

    return new Alignment(std::move(matched_chunks), total_score);
}

Alignment * alignment_from_cell(const AlignTask& task, ScoreCell* p_cell) {
  
  const size_t m = task.query->size() + 1;

  ScoreCellPVec trail;
  trail.reserve(m);
  build_trail(p_cell, trail);

  Alignment * aln = alignment_from_trail(task, trail);
  return aln;
}


// Make partial sums of the preceeding fragment sizes, up to (missed_sites + 1) fragments.
PartialSums make_partial_sums(const IntVec& frags, const int missed_sites) {
/*

  Consider fragments with indices i and fragment sizes f
  ...|---------|----------|--------|...
       i-2        i-1        i
       f_(i-2)    f_(i-1)   f_(i)

  The partial sums fragment i for the case missed_sites = 2 will be:
   [ f_i, f_i + f_(i-1), f_i + f_(i-1) + f(i-2)] 

 */
  const IntVec zero_sums(missed_sites + 1, 0);
  const int num_frags = frags.size();

  PartialSums partial_sums(num_frags, zero_sums);

  for (int i = 0; i < num_frags; i++) {
    IntVec& ps = partial_sums[i];
    const int lower_index = i - missed_sites; // inclusive
    int ind = 0;
    int cur_sum = 0;
    for (int j = i, ind = 0; j >= lower_index; j--, ind++) {
      if (j < 0) break;
      cur_sum += frags[j];
      ps[ind] = cur_sum;
    }
  }

  return partial_sums;
}

// Make partial sums of the preceeding fragment sizes, up to (missed_sites + 1) fragments.
PartialSums* make_partial_sums_new(const IntVec& frags, const int missed_sites) {
/*

  Consider fragments with indices i and fragment sizes f
  ...|---------|----------|--------|...
       i-2        i-1        i
       f_(i-2)    f_(i-1)   f_(i)

  The partial sums fragment i for the case missed_sites = 2 will be:
   [ f_i, f_i + f_(i-1), f_i + f_(i-1) + f(i-2)] 

 */
  const IntVec zero_sums(missed_sites + 1, 0);
  const int num_frags = frags.size();

  PartialSums* p_partial_sums = new PartialSums(num_frags, zero_sums);

  for (int i = 0; i < num_frags; i++) {
    IntVec& ps = (*p_partial_sums)[i];
    const int lower_index = i - missed_sites; // inclusive
    int ind = 0;
    int cur_sum = 0;
    for (int j = i, ind = 0; j >= lower_index; j--, ind++) {
      if (j < 0) break;
      cur_sum += frags[j];
      ps[ind] = cur_sum;
    }
  }

  return p_partial_sums;
}


//////////////////////////////////////////////////////////
// Fill score matrix, find best alignment, and return it.
Alignment * make_best_alignment(const AlignTask& task) {

  // populate the score matrix
  fill_score_matrix(task);

  // get the best alignment.
  ScoreCellPVec trail;
  bool have_alignment = get_best_alignment(task, trail);
  if (!have_alignment) {
    return nullptr;
  }

  return alignment_from_trail(task, trail);
}

// Fill score matrix, find best alignment, and return it.
Alignment * make_best_alignment_using_partials(const AlignTask& task) {

  // populate the score matrix
  fill_score_matrix_using_partials(task);

  // get the best alignment.
  ScoreCellPVec trail;
  bool have_alignment = get_best_alignment(task, trail);
  if (!have_alignment) {
    return nullptr;
  }

  return alignment_from_trail(task, trail);
}

//////////////////////////////////////////////////////////
// Fill score matrix, find best n alignments, and return them in a vector..
int make_best_alignments_using_partials(const AlignTask& task) {

  // Unpack the alignment task
  const IntVec& query = *task.query;
  const IntVec& ref = *task.ref;
  const AlignOpts& align_opts = *task.align_opts;

  // populate the score matrix
  fill_score_matrix_using_partials(task);

  // get the best alignments and store the results in the task.
  int num_alignments = get_best_alignments(task);
  return num_alignments;
}


std::ostream& operator<<(std::ostream& os, const Chunk& chunk) {

  os << "([" << chunk.start << ", " << chunk.end << "], " 
     << chunk.size << ")";

  return os;

}

std::ostream& operator<<(std::ostream& os, const MatchedChunk& chunk) {

  os << "q: " << chunk.query_chunk
     << " r: " << chunk.ref_chunk
     << " score: " << chunk.score << "\n";

  return os;
}

std::ostream& operator<<(std::ostream& os, const Score& score) {

  os << "(" << score.query_miss_score
     << ", " << score.ref_miss_score
     << ", " << score.sizing_score
     << ")";

  return os;
}

std::ostream& operator<<(std::ostream& os, const Alignment& aln) {

  os << "Alignment: " << aln.matched_chunks << "\n";
  return os;

}

std::ostream& print_align_task(std::ostream& os, const AlignTask& task) {

  os << "align_task:\n"
     << "\tquery: " << task.query << "\n"
     << "\tref: " << task.ref << "\n"
     << std::endl;

  return os;
}
