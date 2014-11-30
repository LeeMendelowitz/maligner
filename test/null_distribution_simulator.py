# Some code to test the null distribution. We read in the human
# reference map and simulate random chunks with missed sites.
# We output the probability distribution of unmatched sites over random
# chunks of a given size.
#
######################################################################

import malignpy.core.null_distribution
from malignpy.core.null_distribution import *
from malignpy.maps.utils import read_maps
import numpy as np
import pandas
from pandas import DataFrame
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

N = 10000000
sim_results = simulator.generate_random_frags(N)

o = np.argsort(sim_results.sizes)
sim_results.sizes = sim_results.sizes[o]
sim_results.misses = sim_results.misses[o]
max_misses = np.max(sim_results.misses)

# Determine the distribution of misses sites among randomly 
# generated fragments by the size of the fragment.
start = 2000
end = 100000 + 0.1
by = 2000
rel_error = 0.05
min_error = 1000 # bp
test_points = np.arange(start, end, by)
num_test_points = len(test_points)
lb = np.maximum(test_points - np.maximum(rel_error*test_points, min_error), 0)
ub = test_points + np.maximum(rel_error*test_points, min_error)
lbi = np.searchsorted(sim_results.sizes, lb, side = 'left')
ubi = np.searchsorted(sim_results.sizes, ub, side = 'right')
misses_at_point = [sim_results.misses[lbi[i]:ubi[i]] for i in xrange(num_test_points)]
misses_at_point_bincount = np.array([np.bincount(m, minlength = max_misses + 1) for m in misses_at_point])
nr,nc = misses_at_point_bincount.shape

num_in_bin = np.array([ubi[i] - lbi[i] for i in xrange(num_test_points)])

# normalize the bincount to make fractions
denom = num_in_bin.reshape(len(num_in_bin), 1).repeat(nc, axis=1).astype(np.float64)
misses_at_point_frac = misses_at_point_bincount / denom

# Create a pandas data frame summarizing the results of the simulation
df_data = {'midpt' : test_points,
           'n' : num_in_bin}
df = pandas.DataFrame(data = df_data)
ncol = misses_at_point_bincount.shape[1]
for i in range(ncol):
    k = 'miss%i'%i
    df[k] = misses_at_point_bincount[:,i]
    k = 'frac%i'%i
    df[k] = misses_at_point_frac[:,i]

df.to_csv('sim.results.csv', foat_format = "%.6f", index = False)

# Write a pandas array to summarize the distribution.

# We can use the scipy rv_discrete distribution to simulate a distribution of misses site rate for 
# generating a random map compatible with the query:
# http://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.rv_discrete.html#scipy.stats.rv_discrete


