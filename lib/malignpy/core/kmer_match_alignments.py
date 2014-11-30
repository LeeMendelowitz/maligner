# Objects for parsing output of kmer_match alignments
import numpy as np
from numpy import array, sum

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
    self.start = f[3]
    self.end = f[4]

    frags_raw = f[5:]
    self.frags = [parse_frag(f) for f in frags_raw]
    self._frags_raw = frags_raw

    self.misses = array([f[1] - f[0] - 1 for f in self.frags])
    self.frag_lengths = array([f[2] for f in self.frags])

    self.total_misses = sum(self.misses)
    self.num_frags = len(self.frags)

    self.miss_rate = float(self.total_misses / (self.total_misses + self.num_frags))

