#!/usr/bin/env python
"""
Produce summary statistics for a maps file.
"""

import argparse, sys, os
import numpy as np
from numpy import empty, array, concatenate, mean, percentile, min, max
from malignpy.maps.SOMAMap import SOMAMap
from malignpy.maps.utils import gen_maps
from malignpy.common import wrap_file_function

@wrap_file_function('r')
def run(input_maps_file):
  map_gen = gen_maps(input_maps_file)
  frags_per_map = []
  frags = []
  lengths = []

  sys.stderr.write("reading maps...\n")
  for map in map_gen:
    frags_per_map.append(map.numFrags)
    frags.extend(map.frags)
    lengths.append(map.length)

  sys.stderr.write("computing summary...\n\n")

  frags_per_map = array(frags_per_map)
  frags = array(frags)
  lengths = array(lengths)
  probs = [25, 50, 75]

  def print_stats(vals, pfx):
    i = len(probs)
    p = percentile(vals, probs)
    for i in range(i):
      row_name = '{}_Q{:02d}'.format(pfx, probs[i])
      print '{}\t{:,.2f}'.format(row_name, p[i])

    print '{}_mean\t{:,.2f}'.format(pfx, mean(vals))
    print '{}_min\t{:,.2f}'.format(pfx, min(vals))
    print '{}_max\t{:,.2f}'.format(pfx, max(vals))



  print_stats(lengths, 'map_length_bp')
  print '\n'
  print_stats(frags_per_map, 'frags_per_map')
  print '\n'
  print_stats(frags, 'frag_length_bp')
  print '\n'



if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="""Produce summary statistics for a maps file.""")
  parser.add_argument('maps_file', metavar='MAPS_FILE', type=str,
                   help='Input maps file in the maligner maps format.')
  args = parser.parse_args()
  run(args.maps_file)
