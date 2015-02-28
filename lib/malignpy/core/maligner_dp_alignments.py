# Utilities for parsing a maligner_dp alignments file
from itertools import izip
import numpy as np
from numpy import array

HEADER = "query_map ref_map is_forward  num_matched_chunks  query_misses  ref_misses  query_miss_rate ref_miss_rate total_score total_rescaled_score  sizing_score  sizing_score_rescaled query_scaling_factor  chunk_string"

TYPES = [str, str, str, int, int, int, float, float, float, float, float, float, float, str]
FIELDS = HEADER.strip().split()
OUTPUT_FIELDS = FIELDS + ["ref_start", "ref_end", "ql", "rl", "max_chi2", "median_chi2"]


class Chunk(object):

  def __init__(self, t):
    """Initialize from a tuple"""
    self.qs = int(t[0])
    self.qe = int(t[1])
    self.ql = int(t[2])
    self.rs = int(t[3])
    self.re = int(t[4])
    self.rl = int(t[5])
    self.is_boundary = False

  def __str__(self):
    ks = ['qs', 'qe', 'ql','rs','re','rl']
    return str(tuple(getattr(self, k) for k in ks))
  def __repr__(self):
    return str(self)

class Alignment(object):
  """A MalignerDP Alignment"""

  def __init__(self, line, has_boundary_chunks = True):
    """Initialize Alignment object from a line in the maligner_dp output"""

    fields = line.strip().split()

    for k,t,f in izip(FIELDS, TYPES, fields):
      setattr(self, k, t(f))

    # Parse the chunk_string
    chunks = (c for c in self.chunk_string.split(';') if c)
    self.chunks = [Chunk(c.split(',')) for c in chunks]

    self.ref_start = -1
    self.ref_end = -1

    if self.chunks:

      self.ref_start = self.chunks[0].rs
      self.ref_end = self.chunks[1].re

      if self.is_forward == "F":
        self.query_start = self.chunks[0].qs
        self.query_end = self.chunks[-1].qe
      else:
        self.query_start = self.chunks[-1].qs
        self.query_end = self.chunks[0].qe


  def compute_stats(self, has_boundary_chunks = True, min_sd = 500, sd_rate = 0.05):
    """Compute additional alignment characteristics"""

    if self.chunks and has_boundary_chunks:
      self.chunks[0].is_boundary = True
      self.chunks[-1].is_boundary = True

    # Compute other characteristics of the alignment
    self.qls = array([c.ql for c in self.chunks])
    self.rls = array([c.rl for c in self.chunks])

    if has_boundary_chunks:
      self.ql = sum(self.qls[1:-1])
      self.rl = sum(self.rls[1:-1])
    else:
      self.ql = sum(self.qls)
      self.rl = sum(self.rls)

    # Absolute error
    self.abs_error = np.abs(self.qls - self.rls)

    # Compute the chi2 score
    sd = sd_rate * self.rls
    sd = np.where(sd < min_sd, min_sd, sd)
    self.chi2 = (self.abs_error/sd)**2
    if has_boundary_chunks:
      self.chi2 = self.chi2[1:-1]

    self.median_chi2 = np.median(self.chi2)
    self.max_chi2 = np.max(self.chi2)


  def get_data(self):
    return {k: getattr(self, k) for k in FIELDS}

  def get_output_data(self):
    return {k: getattr(self, k, None) for k in OUTPUT_FIELDS}






