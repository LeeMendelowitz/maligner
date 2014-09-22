"""
Wrap malignpy C++ library to make it more python friendly.
"""

from .malignpy import *
from utils import is_iterable
from collections import OrderedDict

import numpy as np
log10 = np.log10
log = np.log

LOG_1_OVER_SQRT_2PI = -0.5*log(2 * np.pi)

# Make the IntVec more friendly
def _IntVec_create(iterable=None):
  """
  Initialize integer vector. If iterable is 
  provided, copy the values into the vector.
  """
  res = make_int_vec()
  if iterable:
    vals = (int(i + 0.5) for i in iterable)
    res.extend(vals)
  return res

IntVec.create = staticmethod(_IntVec_create)
############################################
def get_data(obj):
  d = getattr(obj, 'get_data', None)
  if d is not None:
    return d()
  return obj

def get_data_from_obj(obj, k):
  obj = getattr(obj, k, None)
  if obj is None:
    return None
  return get_data(obj)

#################################################################
# Define data functions which return a dictionary representation
# for various objects.
def _chunk_get_data(self):
    keys = ['start', 
    'end',
    'size',
    'is_boundary']
    return { k: getattr(self, k) for k in keys }
Chunk.get_data = _chunk_get_data


def _score_get_data(self):
  keys = ['query_miss_score',
          'ref_miss_score',
          'sizing_score']
  return { k: getattr(self, k) for k in keys }
Score.get_data = _score_get_data


def _matched_chunk_data(self):
  data = {}
  data['query_chunk'] = self.query_chunk.get_data()
  data['ref_chunk'] = self.ref_chunk.get_data()
  data['score'] = self.score.get_data()
  return data
MatchedChunk.get_data = _matched_chunk_data


#############################################################################################################
def _alignment_set_data(self, query_id, ref_id, query_num_frags, query_num_sites, ref_is_forward=True, **kwargs):
  """
  Set alignment metadata.
  """
  self.query_id = query_id
  self.ref_id = ref_id
  self.ref_is_forward = ref_is_forward
  self.query_num_frags = query_num_frags
  self.query_num_sites = query_num_sites
  for k,v in kwargs:
    setattr(self,k,v)

def _alignment_get_data(self):
  """Return alignment data as a python dictionary.
  """

  data = OrderedDict()

  # Add these attributes
  keys = [ 'query_id',
          'ref_id', 'ref_is_forward', 
          'score', 'rescaled_score',
          'query_num_sites', 'query_num_frags',
          'num_matched_sites', 'query_misses', 'ref_misses',
          'query_miss_rate', 'ref_miss_rate', 'total_miss_rate',
          'query_interior_size', 'ref_interior_size', 'interior_size_ratio',
          'query_scaling_factor',
          'query_id', 'ref_id', 'ref_is_forward',
          'hit_log_likelihood',
          'miss_log_likelihood',
          'total_log_likelihood',
          'chunk_log_likelihood',
          'chunk_log_likelihoods']

  data.update((k, get_data_from_obj(self, k)) for k in keys)
  data['matched_chunks'] = [m.get_data() for m in self.matched_chunks]
  data['rescaled_matched_chunks'] = [m.get_data() for m in self.rescaled_matched_chunks]
  data['total_score'] = self.score.total
  data['total_score_rescaled'] = self.rescaled_score.total
  data['max_chunk_sizing_score'] = max(c['score']['sizing_score'] for c in data['matched_chunks'])


  return data


def _compute_log_likelihood(self, nick_rate, min_sd, sd_rate, sd_scale = 1000.0):
  """Compute the log likelihood score of an alignment.

  Assumes a Bernoulli process for nicking, and a sizing error model
  where sigma(r_i) = sd_rate * r_i. Sizing error is normally distributed
  with mean 0.

  Parameters:
     - nick_rate
     - min_sd
     - sd_rate
     - sd_scale: Scaling on the sd. in the likelihood evaluation. 
           if sd_scale = 1000.0, then sd is in kb units.
  """

  num_misses = self.ref_misses
  num_hits = self.num_matched_sites
  
  query_chunk_sizes = np.array([c.query_chunk.size for c in self.rescaled_matched_chunks if not c.is_boundary])
  ref_chunk_sizes = np.array([c.ref_chunk.size for c in self.rescaled_matched_chunks if not c.is_boundary])

  size_deltas = query_chunk_sizes - ref_chunk_sizes
  sds = sd_rate * ref_chunk_sizes
  sds[sds < min_sd] = min_sd

  chunk_log_likelihoods = LOG_1_OVER_SQRT_2PI - log(1.0/sd_scale * sds) + ( - size_deltas * size_deltas / (2.0 * sds * sds) )
  
  hit_log_likelihood = num_hits * log(nick_rate)
  miss_log_likelihood = num_misses * log(1.0 - nick_rate)
  chunk_log_likelihood = np.sum(chunk_log_likelihoods)
  total_log_likelihood = hit_log_likelihood + miss_log_likelihood + chunk_log_likelihood

  self.hit_log_likelihood = hit_log_likelihood
  self.miss_log_likelihood = miss_log_likelihood
  self.chunk_log_likelihood = chunk_log_likelihood
  self.chunk_log_likelihoods = [cll for cll in chunk_log_likelihoods]
  self.total_log_likelihood = total_log_likelihood


def _AlignOpts__str__(self):
  attrs = ['query_miss_penalty',
           'ref_miss_penalty',
           'query_max_misses',
           'ref_max_misses',
           'sd_rate',
           'max_chunk_sizing_error',
           'min_sd',
           'alignments_per_reference',
           'min_alignment_spacing',
           'neighbor_delta',
           'query_is_bounded',
           'sd_likelihood_scale'
          ]

  outs = ''
  for a in attrs:
    outs += '%s: %s\n'%(a, getattr(self, a, 'N/A'))
  return outs


# Attach functions to the Alignment class
Alignment.set_data = _alignment_set_data
Alignment.get_data = _alignment_get_data
Alignment.compute_log_likelihood = _compute_log_likelihood


# Add docstrings
AlignOpts.__doc__ = "Alignment Options"
AlignOpts.__init__.__func__.__doc__ = \
"""Initialize a set of alignment options.
  
   Args:

      - query_miss_penalty: Cost of missing have an extra unmatched
         site in query
      - ref_miss_penalty: Cost of having an extra unmatched site
        in reference.
      - query_max_misses: Maximum consecutive unmatched sites in query. 
      - ref_max_misses: Maximum consecutive unmatched sites in reference.
      - sd_rate: Fraction of reference chunk size considered to be 1 s.d. in sizing error model.
      - min_sd: Minimum standard deviation imposed by sizing error model, in bp.
      - max_chunk_sizing_error: Maximum sizing error score allowed
       for a reference chunk.
      - alignments_per_reference: Number of alignments to store per reference.
      - min_alignment_spacing: Minimum number of fragments between rightmost reference
        fragment to accept a secondary alignment. Alignments are selecting in descending
        order of score, so the best alignment is selected first, then the second best alignment
        provided it is min_alignment_spacing fragments away in reference.
      - neighbor_delta: Select alignments within +/- neighbor_delta reference fragments of the best scoring
                      alignments.
      - query_is_bounded: True if query leftmost/rightmost fragment is bounded on both
        ends by a nick site.

"""

AlignOpts.__str__ = _AlignOpts__str__




