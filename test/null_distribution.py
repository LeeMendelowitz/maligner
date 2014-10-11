# Some code to test the null distributionself.

import malignpy.core.null_distribution
reload(malignpy.core.null_distribution)
from malignpy.core.null_distribution import *
from malignpy.maps.utils import read_maps
import numpy as np



# Read human reference maps
map_d = read_maps("test/human.NtBspQI.smoothed.map")
all_frags = []
for m in map_d.itervalues():
  all_frags.extend(m.frags)
all_frags = np.array(all_frags)

controls = NullModelControls()
simulator = NullModelSimulator(all_frags, controls)

map1 = simulator.simulate_query_map()

scorer = NullModelScorer(controls)

#####################################################
ref_map = simulator.simulate_ref_map(num_chunks = 1000)
query_map = simulator.simulate_query_map(num_chunks = 20)