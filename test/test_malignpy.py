import sys, os
from os.path import split, abspath

cwd = os.getcwd()
root_dir = split(abspath(cwd))[0]
lib_dir = os.path.join(root_dir, 'lib')
sys.path.append(lib_dir)

import malignpy as mp

def run_task(align_task):
  mp.fill_score_matrix(align_task)
  trail = mp.ScoreCellVec()
  mp.get_best_alignment(align_task, trail)
  qchunks = mp.ChunkVec()
  rchunks = mp.ChunkVec()
  mp.build_chunk_trail(align_task, trail, qchunks, rchunks)
  #mp.print_chunk_trail(qchunks, rchunks);
  return (qchunks, rchunks)

def run1():
  q = [1, 4, 5, 6, 7, 8, 10, 12, 4, 7]
  q = [v*1000 for v in q]
  r = [7000, 9000, 40000] + q + [36000, 1000]

  # Make a score matrix
  q2 = mp.IntVec()
  q2.assign(q)
  r2 = mp.IntVec()
  r2.assign(r)


  sm = mp.ScoreMatrix(len(q2) + 1, len(r2) + 1)
  align_opts = mp.AlignOpts(3, 5, 2, 2)
  align_task = mp.AlignTask(q2, r2, sm, align_opts)
  res = run_task(align_task)
  return res
