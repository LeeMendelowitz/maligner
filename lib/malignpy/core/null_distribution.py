"""
Build a distribution of alignment scores under the null model 
that alignments are random. This is done by taking two patterns
of fragments of length n, and computing the score of the alignment between them.

The fragments are selected in such a way that there are missed sites in the query
and potentially the reference a simple Bernoulli model.
"""
import numpy as np
from itertools import izip
import logging
logger = logging.getLogger()
logger.setLevel(logging.INFO)

class NullModelControls(object):
  """Controls for running a null model simulation
  """
  query_miss_penalty = 30
  ref_miss_penalty = 3
  query_miss_probability = 0.0 # Probability of having an unmatched site in the query in an alignment
  ref_miss_probability = 0.25 # Probability of having an unmatched site in reference in an alignment
  min_sd = 500 # minimum standard deviation in basepairs
  N = 1000 # number of random maps to use in computation of null distribution
  min_chunks = 1 # minimum number of fragments to simulate in the null distribution.
  max_chunks = 20  # maximum nubmer of chunks to simulate in the null distribution.
  sd_rate = 0.1 # The standard deviation model for sizing error: sd = sd_rate * ref_chunk_size

  def __init__(self, **kwargs):
    for k, v in kwargs.iteritems():
      setattr(self, k, v)

class NullModelDatabase(object):
  pass

class NullModelScorer(object):

  def __init__(self, controls):
    self.controls = controls

  def compute_sds(self, chunks):
    sds = self.controls.sd_rate * chunks
    sds = np.maximum(self.controls.min_sd, sds)
    return sds

  def sizing_score(self, query_chunks, ref_chunks):
    """
    Score the alignment of query_chunks to ref_chunks. Treat ref_chunks as the reference,
    use it to compute the standard deviations.
    """
    assert(len(query_chunks) == len(ref_chunks))
    deltas = ref_chunks - query_chunks
    sds = self.compute_sds(ref_chunks)
    chunk_scores = deltas/sds
    return np.sum(chunk_scores * chunk_scores)

  def score(self, query_chunks, ref_chunks, query_misses = 0, ref_misses = 0):
    s = self.sizing_score(query_chunks, ref_chunks) + self.controls.query_miss_penalty * query_misses + \
            self.controls.ref_miss_penalty * ref_misses
    return s

class NullModelMap(object):

  def __init__(self, chunks, num_misses_per_chunk):
    """
      chunks: list of chunk lengths (bp)
      num_misses_per_chunk: number of unmatched sites in each chunk
    """
    self.chunks = chunks
    self.num_misses_per_chunk = num_misses_per_chunk

class NullModelSimulator(object):

  def __init__(self, frags, controls):
    self.frags = np.array(frags)
    self.controls = controls

  def _simulate_one_map(self, miss_probability):
    """Return a NullModelMap populated with chunks.
       Sample fragments from self.frags using parameters in self.controls
    """

    # Generate a random pattern from frags with self.controls.max_chunks chunks.
    n = self.controls.max_chunks
    n_to_simulate = 10*(n/(1 - miss_probability) + 1)

    # Generate a random pattern until we have enough chunks
    while True:
      bit_pattern = np.random.rand(n_to_simulate) > miss_probability
      if ( np.sum(bit_pattern) >= n + 1 ):
        break
      n_to_simulate = n_to_simulate * 2

    # Generate random fragments from the fragment distribution
    frags = self.frags
    frags_ind = np.random.randint(0, len(frags), len(bit_pattern))
    map_frags = np.array([frags[ind] for ind in frags_ind])
    matched_sites_ind = np.where(bit_pattern)[0]
    if matched_sites_ind[0] != 0:
      matched_sites_ind = np.concatenate(([0], matched_sites_ind))

    start_end = zip(matched_sites_ind[:-1], matched_sites_ind[1:])
    chunks = np.array([ np.sum(map_frags[s:e]) for s, e in start_end ])
    num_misses = np.array([e - s - 1 for s, e in start_end])

    return NullModelMap(chunks, num_misses)

  def simulate_query_map(self):
    miss_prob = self.controls.query_miss_probability
    return self._simulate_one_map(miss_prob)

  def simulate_ref_map(self):
    miss_prob = self.controls.ref_miss_probability
    return self._simulate_one_map(miss_prob)

  def simulate_query_maps(self, n):
    self.full_query_maps = [ self.simulate_query_map() for i in range(n) ]

  def simulate_ref_maps(self, n):
    self.full_ref_maps = [ self.simulate_ref_map() for i in range(n)]

  def run(self):
    
    logging.info("Generating query maps")
    self.simulate_query_maps(self.controls.N)
    
    logging.info("Generating reference maps")
    self.simulate_ref_maps(self.controls.N)

    # Score alignments of random subsets of query to random subsets of reference
    # TODO!




