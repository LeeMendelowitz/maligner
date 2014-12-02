#!/usr/bin/env python
"""Filter alignments by selected the best non-overlapping alignments for each query 
and computing alignment likelihoods.
"""

import argparse
from itertools import islice
import os, sys
import numpy as np
import pandas
from pandas import DataFrame

import malignpy
from malignpy.common import wrap_file_function
import malignpy.core.null_distribution
from malignpy.core.null_distribution import *
from malignpy.maps.utils import read_maps
from malignpy.core.kmer_match_alignments import \
  (Alignment, iter_non_overlapping_alns, iter_query_groups)


parser = argparse.ArgumentParser(description="""Filter alignments by selecting the best non-overlapping alignments and compute 
  alignment probabilities. Write output alignments file as tab delimited file.""")
parser.add_argument('alns_file', metavar='ALIGNMENTS_FILE', type=str,
                   help='Alignments file.')
parser.add_argument('ref_maps', metavar='REF_MAP_FILE', type=str,
                   help='Reference map file.')
parser.add_argument('-o', '--output', metavar='OUTPUT_FILE', help = "Output file name. (Default: STDOUT)")


def serr(msg):
  sys.stderr.write(msg + '\n')

@wrap_file_function('r', 'w')
def filter_alignments_file(aln_file, output_file, simulator):

  alns = iter_non_overlapping_alns(aln_file)

  fout = output_file
  fout.write(Alignment.match_string_header + "\n")
  for i,a in enumerate(alns):
    simulator.assign_likelihoods(a)
    fout.write(a.match_string() + '\n')
  fout.close()


if __name__ == '__main__':

  args = parser.parse_args()

  if args.output is None:
    args.output = sys.stdout

  serr('reading fragments')
  map_d = read_maps(args.ref_maps)
  all_frags = []
  for m in map_d.itervalues():
    all_frags.extend(m.frags)
  all_frags = np.array(all_frags)

  controls = NullModelControls()
  serr('building simulator')
  simulator = NullModelSimulator(all_frags, controls)
  scorer = NullModelScorer(controls)
  serr('simulating null distribution')
  simulator.simulate_null_distribution()
  serr('making interpolators')
  simulator.make_interpolators()
  serr('filtering alignments')
  filter_alignments_file(args.alns_file, args.output, simulator)