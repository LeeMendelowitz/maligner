# Some code to test the null distribution
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
print 'making interpolaters'
simulator.make_interpolaters()

# Generate null distribution for maps
# from malignpy.maps.utils import gen_maps
# maps_file = 'data/MM52_Human_NtBspQI.1000.nmaps'
# maps_gen = gen_maps(maps_file)
# fout = open('MM52_Human_NtBspQI.1000.nmaps.unmatched_probs', 'w')
# for m in maps_gen:
#   res = simulator.simulate_miss_rates(m.frags)
#   fout.write('%s\t%s\n'%(m.mapId, '\t'.join('%6.3f'%f for f in res)))
# fout.close()
