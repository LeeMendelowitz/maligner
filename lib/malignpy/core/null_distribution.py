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

###################################################################
class NullModelControls(object):
  """Controls for running a null model simulation
  """
  query_miss_penalty = 30 # Penalty for having a site in query that is unmatched in an alignment.
  ref_miss_penalty = 3 # Penalty for having a site in reference that is unmatched in an alignment.
  query_miss_probability = 0.25 # Probability of missing a nick site when sampling an RMap read from reference.
  min_sd = 500 # minimum standard deviation in basepairs
  N = 1000 # number of random maps to use in computation of null distribution
  min_chunks = 1 # minimum number of fragments to simulate in the null distribution.
  max_chunks = 20  # maximum nubmer of chunks to simulate in the null distribution.
  sd_rate = 0.1 # The standard deviation model for sizing error: sd = sd_rate * ref_chunk_size

  def __init__(self, **kwargs):
    for k, v in kwargs.iteritems():
      setattr(self, k, v)

###################################################################
class NullModelDatabase(object):
  pass


###################################################################
class NullModelScorer(object):

  def __init__(self, controls):
    self.controls = controls

  def compute_sds(self, chunk_lengths):
    sds = self.controls.sd_rate * chunk_lengths
    sds = np.maximum(self.controls.min_sd, sds)
    return sds

  def sizing_score(self, query_chunk_lengths, ref_chunk_lengths):
    """
    Score the alignment of query_chunks to ref_chunks. Treat ref_chunks as the reference,
    use it to compute the standard deviations.
    """
    assert(len(query_chunk_lengths) == len(ref_chunk_lengths))
    deltas = ref_chunk_lengths - query_chunk_lengths
    sds = self.compute_sds(ref_chunk_lengths)
    chunk_scores = deltas/sds
    return np.sum(chunk_scores * chunk_scores)

  def score(self, query_chunk_lengths, ref_chunk_lengths, query_misses = 0, ref_misses = 0):
    s = self.sizing_score(query_chunk_lengths, ref_chunk_lengths) + \
        self.controls.query_miss_penalty * query_misses + \
        self.controls.ref_miss_penalty * ref_misses
    return s

  def score_map(self, query_map, ref_map):
    query_misses = np.sum(query_map.misses_per_chunk)
    ref_misses = np.sum(ref_map.misses_per_chunk)
    return self.score(query_map.chunk_lengths, ref_map.chunk_lengths,
        query_misses, ref_misses)

  def score_at_position(self, query_map, ref_map, position):
    lr = len(ref_map)
    lq = len(query_map)
    assert(position + lq <= lr)
    adjusted_ref = ref_map[position:position+lq]
    return self.score_map(query_map, adjusted_ref)


###################################################################
class NullModelMap(object):
  """
  A very simple representation of a random Rmap sampled from
  a distribution of restriction fragment lengths.

  The NullModelMap samples missed sites, but not sizing error.
  Each chunk is a length, in bp, with a certain number of interior missed sites.
  """

  def __init__(self, chunk_lengths, misses_per_chunk):
    """
    """
    self.chunk_lengths = chunk_lengths
    self.misses_per_chunk = misses_per_chunk
    assert(len(self.chunk_lengths) == len(self.misses_per_chunk))

  def __getitem__(self, item):
    """This make Map sliceable!"""

    if isinstance(item, slice):
      chunk_lengths = self.chunk_lengths[item]
      misses_per_chunk = self.misses_per_chunk[item]
      return NullModelMap(chunk_lengths, misses_per_chunk)

    raise TypeError("Expected slice!")

  def __len__(self):
    return len(self.chunk_lengths)

  def __str__(self):
    s = "<NullModelMap with {num_chunks} chunks at 0x{mem:x}>"
    return s.format(num_chunks = len(self),
                    mem = id(self))


###################################################################
class NullModelSimulator(object):

  def __init__(self, frags, controls):
    self.frags = np.array(frags)
    self.controls = controls

  def _simulate_one_map(self, miss_probability, num_chunks):
    """Return a NullModelMap populated with chunks.
       Sample fragments from self.frags using parameters in self.controls.
       Sampling includes missed sites (but no sizing error)
    """

    # Generate a random pattern from frags with num_chunks
    hit_prob = 1 - miss_probability
    g = np.random.geometric(hit_prob, num_chunks)
    randint = np.random.randint
    sum = np.sum
    frags = self.frags
    N = len(frags)

    inds = (randint(0, N, num_frags) for num_frags in g)
    chunk_lengths = np.fromiter((sum(frags[i]) for i in inds), np.int)
    misses_per_chunk = g - 1

    return NullModelMap(chunk_lengths, misses_per_chunk)

  def simulate_query_map(self, num_chunks = None):
    miss_prob = self.controls.query_miss_probability
    if not num_chunks:
      num_chunks = self.controls.max_chunks
    return self._simulate_one_map(miss_prob, num_chunks)

  def simulate_ref_map(self, num_chunks = None):
    miss_prob = 0.0
    if not num_chunks:
      num_chunks = self.controls.max_chunks
    return self._simulate_one_map(miss_prob, num_chunks)

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


def all_scores_forward(scorer, query_map, ref_map):
  lq = len(query_map)
  lr = len(ref_map)
  last_pos = lr - lq
  scores = [scorer.score_at_position(query_map, ref_map, p) for p in xrange(last_pos)]
  return scores



