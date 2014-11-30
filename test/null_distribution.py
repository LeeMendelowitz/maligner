# Some code to test the null distributionself.

import malignpy.core.null_distribution
reload(malignpy.core.null_distribution)
from malignpy.core.null_distribution import *
from malignpy.maps.utils import read_maps
import numpy as np
import sys



# Read human reference maps
map_d = read_maps("data/human.NtBspQI.smoothed.map")
all_frags = []
for m in map_d.itervalues():
  all_frags.extend(m.frags)
all_frags = np.array(all_frags)

controls = NullModelControls()
simulator = NullModelSimulator(all_frags, controls)

scorer = NullModelScorer(controls)

#####################################################
num_frags = len(all_frags)

sys.stdout.write("simulating map 1....")
ref_map1 = simulator.simulate_query_map(num_chunks = num_frags)
sys.stdout.write("done!\n")

sys.stdout.write("simulating map 2....")
ref_map2 = simulator.simulate_query_map(num_chunks = num_frags)
sys.stdout.write("done!\n")

# Simulate 1000 query maps, align them to the permuated reference
# to determine the null distribution scores.
num_queries = 1000
query_maps = [simulator.simulate_query_map(num_chunks = 20) for i in range(num_queries)]

N = num_queries
N = 1000

fout = open('scores2.csv', 'w')
for i in range(N):
  s = datetime.now()
  ret = scorer.chunk_score_matrix(query_maps[i], ref_map1)
  line = ",".join(str(s) for s in ret.best_scores)
  fout.write(line + "\n")
  e = datetime.now()
  elapsed = (e-s).total_seconds()
  print i, elapsed

fout.close()
