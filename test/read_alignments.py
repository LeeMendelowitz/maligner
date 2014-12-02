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
from malignpy.core.kmer_match_alignments import \
  (Alignment, iter_non_overlapping_alns, iter_query_groups)


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
ALN_FILE = 'kmer_match/MM52.all.rel_error_10p0.aln.head'
aln_file = open(ALN_FILE)
# num_lines = 10000
# lines = (aln_file.next() for i in xrange(num_lines))
# alns = [Alignment(l) for l in lines]
alns = iter_non_overlapping_alns(aln_file)

fout = open(ALN_FILE + '.filtered', 'w')
fout.write(Alignment.match_string_header + "\n")
for i,a in enumerate(alns):
  simulator.assign_likelihoods(a)
  # a.log_HA = simulator.log_pattern_probability(a.misses)
  # a.log_H0 = np.log(simulator.null_probability_generated(a.frag_lengths, a.misses))
  # a.log_likelihood_ratio = simulator.log_likelihood_ratio(a.frag_lengths, a.misses)
  fout.write(a.match_string() + '\n')
fout.close()



