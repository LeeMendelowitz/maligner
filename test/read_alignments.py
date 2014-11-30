# Some code to test the null distribution on alignments.

if __name__ == "__main__":

  import sys, os
  from os.path import split, abspath

  file_dir = split(abspath(__file__))[0]
  root_dir = split(file_dir)[0]
  lib_dir = os.path.join(root_dir, 'lib')
  sys.path.append(lib_dir)

import malignpy.core.null_distribution
from malignpy.core.null_distribution import *
from malignpy.maps.utils import read_maps
import numpy as np
import pandas
from pandas import DataFrame
import sys


# Read human reference maps
print 'reading fragments'
map_d = read_maps("data/human.NtBspQI.smoothed.map")
all_frags = []
for m in map_d.itervalues():
  all_frags.extend(m.frags)
all_frags = np.array(all_frags)

controls = NullModelControls()
print 'building simulator'
simulator = NullModelSimulator(all_frags, controls)
scorer = NullModelScorer(controls)
print 'simulating null distribution'
simulator.simulate_null_distribution()
print 'making interpolators'
simulator.make_interpolators()


print 'reading 1000 alignments'
from malignpy.core.kmer_match_alignments import Alignment
aln_file = open('/media/sf_Lee/workspace/maligner/kmer_match/parrot.chunk.all.new_output2.aln')
num_lines = 1000
lines = (aln_file.next() for i in xrange(num_lines))
alns = [Alignment(l) for l in lines]

for a in alns:
  a.log_likelihood_ratio = simulator.log_likelihood_ratio(a.frag_lengths, a.misses)

repeats_file = open('/media/sf_Lee/workspace/maligner/kmer_match/parrot.chunk.all.new_output2.repeats.aln')
lines = (repeats_file.next() for i in xrange(num_lines))
repeats_alns = [Alignment(l) for l in lines]
for i,r in enumerate(repeats_alns):
  r.log_likelihood_ratio = simulator.log_likelihood_ratio(r.frag_lengths, r.misses)

for i,r in enumerate(repeats_alns):
  print i,r,r.log_likelihood_ratio

from collections import defaultdict
query_to_alns = defaultdict(list)
for r in repeats_alns:
  query_to_alns[r.query].append(r)