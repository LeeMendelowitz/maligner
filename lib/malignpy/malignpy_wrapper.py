"""
Wrap malignpy C++ library to make it more python friendly.
"""

from .malignpy import *
from utils import is_iterable
from collections import OrderedDict

# Make the IntVec more friendly
def _IntVec_create(iterable=None):
  """
  Initialize integer vector. If iterable is 
  provided, copy the values into the vector.
  """
  res = make_int_vec()
  if iterable:
    res.extend(iterable)
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
  data = OrderedDict()

  # Add these attributes
  keys = [ 'query_id', 'ref_id', 'ref_is_forward', 
          'score', 'rescaled_score',   'query_num_sites', 'query_num_frags', 'num_matched_sites', 'query_misses', 'ref_misses',
          'query_miss_rate', 'ref_miss_rate', 'total_miss_rate',
          'query_interior_size', 'ref_interior_size', 'interior_size_ratio', 'query_scaling_factor',
          'query_id', 'ref_id', 'ref_is_forward']
  data.update((k, get_data_from_obj(self, k)) for k in keys)
  data['matched_chunks'] = [m.get_data() for m in self.matched_chunks]
  data['rescaled_matched_chunks'] = [m.get_data() for m in self.rescaled_matched_chunks]
  data['total_score'] = self.score.total
  data['max_chunk_sizing_score'] = max(c['score']['sizing_score'] for c in data['matched_chunks'])


  return data

def _AlignOpts__str__(self):
  attrs = ['query_miss_penalty',
           'ref_miss_penalty',
           'query_max_misses',
           'ref_max_misses',
           'max_chunk_sizing_error',
           'min_sd',
           'alignments_per_reference',
           'min_alignment_spacing'
          ]

  outs = ''
  for a in attrs:
    outs += '%s: %s\n'%(a, getattr(self, a, 'N/A'))
  return outs



Alignment.set_data = _alignment_set_data
Alignment.get_data = _alignment_get_data

AlignOpts.__str__ = _AlignOpts__str__




