# Utilities for parsing a maligner_dp alignments file
# And adding additional summary statistics about each alignment
from itertools import izip
import numpy as np
from numpy import array

HEADER = "query_map ref_map is_forward  num_matched_chunks  query_misses  ref_misses  query_miss_rate ref_miss_rate total_score total_rescaled_score m_score sizing_score  sizing_score_rescaled query_scaling_factor chunk_string"
TYPES = [str, str, str, int, int, int, float, float, float, float, float, float, float, str]
FIELDS = HEADER.strip().split()
assert(len(FIELDS) == len(TYPES))
OUTPUT_FIELDS = FIELDS + ["ref_start", "ref_end", "ql", "rl", "query_misses_interior", "ref_misses_interior",
  "query_miss_rate_interior", "ref_miss_rate_interior", 
  "query_sites_total", "query_sites_interior",
  "ref_sites_total", "ref_sites_interior",
  "max_chi2", "median_chi2", "query_miss_penalty", "ref_miss_penalty"]

FIELDS_TYPES = zip(FIELDS, TYPES)


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

    self.chi2 = None

  def __str__(self):
    ks = ['qs', 'qe', 'ql','rs','re','rl']
    return str(tuple(getattr(self, k) for k in ks))

  def __repr__(self):
    return str(self)

  @property
  def query_misses(self):
    return self.qe - self.qs - 1

  @property
  def ref_misses(self):
    return self.re - self.rs - 1

  def flip_query_coords(self, n):
    """Flip the coordinates with respect to the query with n fragments"""
    qs = self.qs
    self.qs = n - self.qe
    self.qe = n - qs

  def compute_chi2(self, min_sd, sd_rate):
    if self.is_boundary:
      self.chi2 = 0
      return
    delta = self.ql - self.rl
    sd = max(sd_rate*self.rl, min_sd)
    self.chi2 = (delta/sd)**2

class Alignment(object):
  """A MalignerDP Alignment"""

  def __init__(self, line, has_boundary_chunks = True):
    """Initialize Alignment object from a line in the maligner_dp output"""

    self.has_boundary_chunks = has_boundary_chunks

    fields = line.strip().split()

    for k,t,f in izip(FIELDS, TYPES, fields):
      setattr(self, k, t(f))


    # Parse the chunk_string
    chunks = (c for c in self.chunk_string.split(';') if c)
    self.chunks = [Chunk(c.split(',')) for c in chunks]
    assert(self.num_matched_chunks) == len(self.chunks)


    if self.has_boundary_chunks:

      self.chunks[0].is_boundary = True
      self.chunks[-1].is_boundary = True
      self.interior_chunks = self.chunks[1:-1]
      self.num_interior_chunks = len(self.interior_chunks)
      self.query_sites_total = self.num_matched_chunks - 1 + sum(c.query_misses for c in self.chunks)
      self.ref_sites_total = self.num_matched_chunks - 1 + sum(c.ref_misses for c in self.chunks)
      self.query_sites_interior = self.num_interior_chunks + 1 + sum(c.query_misses for c in self.chunks if not c.is_boundary)
      self.ref_sites_interior = self.num_interior_chunks + 1 + sum(c.ref_misses for c in self.chunks if not c.is_boundary)

    else:
      self.interior_chunks = self.chunks
      self.num_interior_chunks = len(self.interior_chunks)
      self.query_sites_total = self.num_matched_chunks + 1 + sum(c.query_misses for c in self.chunks)
      self.ref_sites_total = self.num_matched_chunks + 1 + sum(c.ref_misses for c in self.chunks)
      self.query_sites_interior = self.num_matched_chunks + 1 + sum(c.query_misses for c in self.chunks if not c.is_boundary)
      self.ref_sites_interior = self.num_matched_chunks + 1 + sum(c.ref_misses for c in self.chunks if not c.is_boundary)


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


  def compute_stats(self, has_boundary_chunks = True, query_miss_penalty = 18.0, ref_miss_penalty = 3.0,
    min_sd = 500, sd_rate = 0.05):
    """Compute additional alignment characteristics"""

    # Compute other characteristics of the alignment
    self.qls = array([c.ql for c in self.chunks])
    self.rls = array([c.rl for c in self.chunks])

    # Compute interior query length and reference length,
    # ignoring boundary chunks.
    if has_boundary_chunks:
      self.ql = np.sum(self.qls[1:-1])
      self.rl = np.sum(self.rls[1:-1])
    else:
      self.ql = np.sum(self.qls)
      self.rl = np.sum(self.rls)

    # Absolute error
    self.abs_error = np.abs(self.qls - self.rls)

    # Interior missed sites
    self.query_misses_interior = sum(c.query_misses for c in self.interior_chunks)
    self.ref_misses_interior = sum(c.ref_misses for c in self.interior_chunks)
    self.query_miss_rate_interior = float(self.query_misses_interior)/self.query_sites_interior
    self.ref_miss_rate_interior = float(self.ref_misses_interior)/self.ref_sites_interior

    # Compute the chi2 score
    sd = sd_rate * self.rls
    sd = np.where(sd < min_sd, min_sd, sd)
    self.chi2 = (self.abs_error/sd)**2
    if has_boundary_chunks:
      self.chi2 = self.chi2[1:-1]

    self.median_chi2 = np.median(self.chi2)
    self.max_chi2 = np.max(self.chi2)

    self.query_miss_penalty = query_miss_penalty * sum(c.query_misses for c in self.chunks)
    self.ref_miss_penalty = ref_miss_penalty * sum(c.ref_misses for c in self.chunks)

    # Compute chunk chi2
    for c in self.chunks:
      c.compute_chi2(min_sd, sd_rate)


  def get_data(self):
    return {k: getattr(self, k) for k in FIELDS}

  def get_output_data(self):
    return {k: getattr(self, k, None) for k in OUTPUT_FIELDS}






