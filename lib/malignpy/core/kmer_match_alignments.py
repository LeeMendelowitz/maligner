# Objects for parsing output of kmer_match alignments
import numpy as np
from numpy import array, sum, zeros, sqrt, mean
from collections import defaultdict
from operator import attrgetter

def parse_frag(s):
  s = s.strip('()')
  fields = s.split(',')
  return tuple(int(f) for f in fields)

class Alignment(object):

  def __init__(self, line):
    f = line.strip().split()
    self.query = f[0]
    self.reference = f[1]
    self.orientation = f[2]
    self.start = int(f[3])
    self.end = int(f[4])

    frags_raw = f[5:]
    self.frags = [parse_frag(f) for f in frags_raw]
    self._frags_raw = frags_raw

    self.misses = array([f[1] - f[0] - 1 for f in self.frags])
    self.frag_lengths = array([f[2] for f in self.frags])

    self.total_misses = sum(self.misses)
    self.num_frags = len(self.frags)
    self.num_chunks = self.num_frags
    self.num_ref_frags = self.end - self.start

    self.miss_rate = float(self.total_misses) / (self.total_misses + self.num_frags + 1)
    
    # Placeholders
    self.log_H0 = None
    self.log_HA = None
    self.log_likelihood_ratio = None

  def set_query_data(self, query_map):
    """The alignment only stores information about the reference and assumes the entire query
    is aligned. Here we attach a query_map to the alignment so we can compute frag sizing error
    """
    self.query_map = query_map
    self.query_frags = array(query_map.frags)
    inner_frags = self.query_frags[1:-1]
    # self.frag_lengths: reference lengths of the alignment chunks
    self.frag_sizing_error = inner_frags - self.frag_lengths
    self.rel_frag_sizing_error = self.frag_sizing_error / self.frag_lengths.astype(float)
    self.mean_rel_frag_sizing_error = np.mean(self.rel_frag_sizing_error)
    self.mean_abs_rel_frag_sizing_error = np.mean(np.abs(self.rel_frag_sizing_error))
    self.sd_rel_frag_sizing_error = np.std(self.rel_frag_sizing_error)

  match_string_header = "\t".join([
    "query",
    "reference",
    "orientation",
    "start",
    "end",
    "num_ref_frags",
    "num_chunks",
    "total_misses",
    "miss_rate",
    "log_HA",
    "log_H0",
    "log_likelihood_ratio",
    "E_H0",
    "p_H0",
    "p_H0_alt",
    "mean_rel_frag_sizing_error",
    "mean_abs_rel_frag_sizing_error",
    "sd_rel_frag_sizing_error",
    "frag_string"
  ])

  def match_string(self):
    """Return a string representation of the alignment, including total_misses,
    miss_rate, log_H0, log_HA, and log_likelihood_ratio"""

    def float_formatter(s):
      try:
        return '%.6f'%s
      except TypeError:
        return "NA"

    def float_scientific(s):
      try:
        return '%.6E'%s
      except TypeError:
        return "NA"

    def g(k):
      return getattr(self, k, "NA")

    ret = "\t".join([self.query,
                     self.reference,
                     self.orientation,
                     str(self.start),
                     str(self.end),
                     str(self.num_ref_frags),
                     str(self.num_chunks),
                     str(self.total_misses),
                     float_formatter(self.miss_rate),
                     float_formatter(g("log_HA")),
                     float_formatter(g("log_H0")),
                     float_formatter(g("log_likelihood_ratio")),
                     float_scientific(g("E_H0")),
                     float_scientific(g("p_H0")),
                     float_scientific(g("p_H0_alt")),
                     float_formatter(g("mean_rel_frag_sizing_error")),
                     float_formatter(g("mean_abs_rel_frag_sizing_error")),
                     float_formatter(g("sd_rel_frag_sizing_error")),
                     ' '.join(self._frags_raw)])
    return ret


class ReferenceCover(object):
  """Mark reference fragments that are covered"""

  def __init__(self):
    self.ref_to_covered = defaultdict(lambda: zeros(0))

  def reset(self):
    for k, a in self.ref_to_covered.iteritems():
      self.ref_to_covered[k] = zeros(a.shape[0])

  def _get(self, ref, start_ind, end_ind):
    """Get the array for the ref. resize if necessary."""
    a = self.ref_to_covered[ref]
    if end_ind > a.shape[0]:
      a.resize(end_ind, refcheck=False)
      self._save(ref, a)
    return a

  def _save(self, ref, a):
    """Save the array for the ref"""
    self.ref_to_covered[ref] = a

  def mark(self, ref, start_ind, end_ind):
    a = self._get(ref, start_ind, end_ind)
    a[start_ind:end_ind] = 1
    self._save(ref, a)

  def is_marked(self, ref, start_ind, end_ind):
    a = self._get(ref, start_ind, end_ind)
    return np.any(a[start_ind:end_ind])

  def mark_and_check_overlap(self, ref, start_ind, end_ind):
    """Check if the interval is already covered, then mark the interval.
    return True if the interval was already covered"""
    a = self._get(ref, start_ind, end_ind)
    was_covered = np.any(a[start_ind:end_ind])
    a[start_ind:end_ind] = 1
    self._save(ref, a)
    return was_covered

def iter_query_groups(f):
  """
  Return lists of alignments by query.
  This assumes that alignments for the same query appear consecutively in the alignments file.
  """
  
  lines = (l.strip() for l in f)
  lines = (l for l in lines if l)
  alns = (Alignment(l) for l in lines)

  cur_alns = []
  cur_query = None
  for a in alns:
    
    if a.query == cur_query:
      cur_alns.append(a)
      continue

    if cur_alns:
      yield cur_alns

    cur_alns = [a]
    cur_query = a.query

  if cur_alns:
    yield cur_alns


def iter_non_overlapping_alns(f):
  """Return non-overlapping alignments by query.
     If there are overlapping conflicts, return the 
     alignment with the lower miss_rate. Break ties arbitrarily.

     TODO: Break ties by likelihood?
  """

  key_func = attrgetter('miss_rate')
  ref_cover = ReferenceCover()

  for alns in iter_query_groups(f):
    
    if len(alns) == 1:
      yield alns[0]
      continue

    
    alns = sorted(alns, key = key_func)

    ref_cover.reset()
    for a in alns:
      was_covered = ref_cover.mark_and_check_overlap(a.reference, a.start, a.end)
      if not was_covered:
        yield a








