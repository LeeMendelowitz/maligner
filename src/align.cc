#include <vector>
#include <cassert>
#include <iostream>
using std::cerr;

#include "utils.h"
#include "types.h"
#include "align.h"
#include "globals.h"
#include "ScoreMatrix.h"
#include "ScoreCell.h"

using Constants::INF;

#define DEBUG 1

double sizing_penalty(int query_size, int ref_size, const AlignOpts& align_opts) {
  /* TODO: This can be baked into the dynamic programming routine */
  double delta = query_size - ref_size;
  double sd = 0.1*ref_size;
  double penalty = delta*delta/(sd*sd);
  return penalty;
}

/*
  Populate a score matrix using dynamic programming for ungapped alignment.

  The score matrix should have the same number of columns as the reference.

  The ScoreMatrix should already have the same nubmer of columns as the reference,
  and should have enough rows to accomodate the query.
*/
void fill_score_matrix(AlignTask& align_task) {

  // Unpack the alignment task
  const IntVec& query = align_task.query;
  const IntVec& ref = align_task.ref;
  ScoreMatrix& mat = align_task.mat;
  AlignOpts& align_opts = align_task.align_opts;


  const int m = query.size() + 1;
  const int n = ref.size() + 1;

  // Note: Number of rows may be different from m if matrix is padded with extra rows.
  const int num_rows = mat.getNumRows();

  assert((int) mat.getNumCols() == n);
  assert((int) mat.getNumRows() >= m);

  #if DEBUG > 0
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


  for (int j = 1; j < n; j++) {
    
    bool is_query_boundary = (j==1) || (j==n-1);
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

        bool is_ref_boundary = (l == 0) || (j == n);
        int ref_miss = j - l - 1;
        ref_size += ref[l];
        double ref_miss_score = ref_miss * align_opts.ref_miss_penalty;

        int query_size = 0;
        const int offset_back = num_rows*l;

        for(int k = i-1; k >= k0; k--) {

          #if DEBUG > 0
          cerr << "i: " << i
               << " j: " << j
               << " k: " << k
               << " l: " << l
               << "\n";
          #endif


          query_size += query[k];
          ScoreCell* pTarget = mat.getCell(offset_back + k);
          if (pTarget->score_ == -INF) continue;

          int query_miss = i - k - 1;
          double query_miss_score = query_miss * align_opts.query_miss_penalty;

          // Score the extension
          double chunk_score = -sizing_penalty(query_size, ref_size, align_opts) - query_miss_score - ref_miss_score;

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
  
} // fill_score_matrix



bool get_best_alignment(AlignTask& task, ScoreCellPVec& trail) {

  // Go to the last row of the ScoreMatrix and identify the best score.
  const int m = task.query.size() + 1;
  const int n = task.ref.size() + 1;
  const int num_rows = task.mat.getNumRows();
  const int num_cols = task.mat.getNumCols();
  const int last_row = m - 1;

  double best_score = -INF;
  ScoreCell * p_best_cell = nullptr;
  int index = last_row;

  // Get the cell with the best score in the last row.
  for (int i = 0; i < n; i++, index += num_rows) {
    
    ScoreCell * pCell = task.mat.getCell(index);

    #if DEBUG > 0
    cerr << "index: " << index << ", ";
    cerr << "cell: " << *pCell << "\n";
    #endif

    if (pCell->score_ > best_score) {
      p_best_cell = pCell;
      best_score = pCell->score_;
    }
  }

  #if DEBUG > 0
    if (p_best_cell) {
        cerr << "\np_best_cell: " << *p_best_cell <<  "\n";
    } else {
        cerr << "\np_best_cell: " << p_best_cell <<  "\n";
    }
  #endif
 
  if (p_best_cell && p_best_cell->score_ > -INF) {
    // Get the traceback
    trail.reserve(m);
    build_trail(p_best_cell, trail);
    return true;
  }

  return false;

}

void build_trail(ScoreCell* pCell, ScoreCellPVec& trail) {
  ScoreCell* pCur = pCell;
  while (pCur != nullptr) {
    trail.push_back(pCur);
    pCur = pCur->backPointer_;
  }
}


// trail: starts from end of alignment
void build_chunk_trail(AlignTask& task, ScoreCellPVec& trail, ChunkVec& query_chunks, ChunkVec& ref_chunks) {

  const size_t ts = trail.size();
  if (ts == 0) { return; }

  query_chunks.clear();
  ref_chunks.clear();
  query_chunks.reserve(trail.size());
  ref_chunks.reserve(trail.size());

  const IntVec& query = task.query;
  const IntVec& ref = task.ref;

  ScoreCell* pLast = trail[0];
  int ml = pLast->q_;
  int nl = pLast->r_;

  for(size_t i = 1; i < ts; i++) {
    
    ScoreCell* pCell = trail[i];
    int m = pCell->q_; // index of query site, one based
    int n = pCell->r_; // index of ref size, one based

    #if DEBUG > 0
    cerr << "cell: " << *pCell << "\n"
         << " i: " << i 
         << " m: " << m << " ml: " << ml
         << " n: " << n << " nl: " << nl
         << "\n";
    #endif

    assert(m < ml);
    assert(n < nl);

    // Build chunks
    int q_size = sum(query, m, ml);
    int r_size = sum(ref, n, nl);
    Chunk q_chunk(m, ml, q_size); 
    Chunk r_chunk(n, nl, r_size);
    query_chunks.push_back(q_chunk);
    ref_chunks.push_back(r_chunk);
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



std::ostream& operator<<(std::ostream& os, const Chunk& chunk) {

  os << "([" << chunk.start << ", " << chunk.end << "], " 
     << chunk.size << ")";

  return os;

}