"""
Build a distribution of alignment scores under the null model 
that alignments are random. This is done by taking two patterns
of fragments of length n, and computing the score of the alignment between them.

The fragments are selected in such a way that there are missed sites in the query
and potentially the reference a simple Bernoulli model.
"""

import numpy as np
from numpy import bincount, cumsum, log
from scipy.interpolate import interp1d
from scipy.stats import rv_discrete
from itertools import izip
import sys
import logging
import pandas
from datetime import datetime

from malignpy.maps.SOMAMap import SOMAMap
logger = logging.getLogger()
logger.setLevel(logging.INFO)

def stdout(msg):
  sys.stdout.write(msg)
  sys.stdout.flush()

class ReturnObj(object):
  pass

###################################################################
class NullModelControls(object):
  """Controls for running a null model simulation
  """
  query_miss_penalty = 30 # Penalty for having a site in query that is unmatched in an alignment.
  ref_miss_penalty = 3 # Penalty for having a site in reference that is unmatched in an alignment.
  query_miss_probability = 0.25 # Probability of missing a nick site when sampling an RMap read from reference.
  min_sd = 500 # minimum standard deviation in basepairs

  N = 1000 # number of random maps to use in computation of null distribution
  min_chunks = 1 # minimum number of fragments to simulate in the null distribution.
  max_chunks = 20  # maximum number of chunks to simulate in the null distribution.
  sd_rate = 0.1 # The standard deviation model for sizing error: sd = sd_rate * ref_chunk_size

  num_frags = 10000000 # Number of random fragments to generate for null distribution
  max_misses = 15 # The number of maximum consecutive misses to generate for the null distribution.
  rel_error = 0.05 # The relative error allowed for fragment compatability
  min_error = 1000 # The minimum error considered for fragment compatibility.

  sim_results_csv = 'sim.results.csv'

  def __init__(self, **kwargs):
    for k, v in kwargs.iteritems():
      setattr(self, k, v)

###################################################################
class NullModelDatabase(object):
  pass


###################################################################
class NullModelScorer(object):

  def __init__(self, controls):
    self.controls = controls

  def compute_sds(self, chunk_lengths):
    sds = self.controls.sd_rate * chunk_lengths
    sds = np.maximum(self.controls.min_sd, sds)
    return sds

  def sizing_score(self, query_chunk_lengths, ref_chunk_lengths):
    """
    Score the alignment of query_chunks to ref_chunks. Treat ref_chunks as the reference,
    use it to compute the standard deviations.
    """
    assert(len(query_chunk_lengths) == len(ref_chunk_lengths))
    deltas = ref_chunk_lengths - query_chunk_lengths
    sds = self.compute_sds(ref_chunk_lengths)
    chunk_scores = deltas/sds
    return np.sum(chunk_scores * chunk_scores)

  def score(self, query_chunk_lengths, ref_chunk_lengths, query_misses = 0, ref_misses = 0):
    s = self.sizing_score(query_chunk_lengths, ref_chunk_lengths) + \
        self.controls.query_miss_penalty * query_misses + \
        self.controls.ref_miss_penalty * ref_misses
    return s

  def score_map(self, query_map, ref_map):
    query_misses = np.sum(query_map.misses_per_chunk)
    ref_misses = np.sum(ref_map.misses_per_chunk)
    return self.score(query_map.chunk_lengths, ref_map.chunk_lengths,
        query_misses, ref_misses)

  def score_at_position(self, query_map, ref_map, position):
    lr = len(ref_map)
    lq = len(query_map)
    assert(position + lq <= lr)
    adjusted_ref = ref_map[position:position+lq]
    return self.score_map(query_map, adjusted_ref)


  def chunk_score_matrix(self, query_map, ref_map):
    """Return an np.array with the scores of matching query chunk q_i with
    reference chunk r_j.
    """
    ref_chunk_lengths = ref_map.chunk_lengths
    ref_misses_per_chunk = ref_map.misses_per_chunk
    query_chunk_lengths = query_map.chunk_lengths
    ref_sds = self.compute_sds(ref_chunk_lengths)
    chunk_score_mat = np.zeros((len(query_chunk_lengths), len(ref_chunk_lengths)))
    nrow = len(query_chunk_lengths)
    ncol = len(ref_chunk_lengths)
    sd_mat = np.tile(ref_sds, (nrow, 1))
    r_mat = np.tile(ref_chunk_lengths, (nrow, 1))
    q_mat = np.tile(query_chunk_lengths.reshape(nrow, 1), (1, ncol))
    delta = r_mat - q_mat

    chunk_score_matrix = delta*delta / (sd_mat*sd_mat)

    # Compute the alignment score matrix,
    # where entry (i,j) has the chunk score of the ith chunk of the alignment starting
    # at reference location j
    m,n = chunk_score_matrix.shape
    nr, nc = m, n - m + 1
    r = np.arange(nr).repeat(nc)
    c = np.tile(np.arange(nc), nr) + np.arange(nr).repeat(nc)
    alignment_score_matrix = chunk_score_matrix[r,c].reshape(nr, nc)

    ret = ReturnObj()
    ret.chunk_score_matrix = chunk_score_matrix
    ret.alignment_score_matrix = alignment_score_matrix
    ret.alignment_score_matrix_cumsum = np.cumsum(alignment_score_matrix, axis=0)
    ret.best_scores = np.min(ret.alignment_score_matrix_cumsum, axis = 1)


    return ret






###################################################################
class NullModelMap(object):
  """
  A very simple representation of a random Rmap sampled from
  a distribution of restriction fragment lengths.

  The NullModelMap samples missed sites, but not sizing error.
  Each chunk is a length, in bp, with a certain number of interior missed sites.
  """

  def __init__(self, chunk_lengths, misses_per_chunk):
    """
    """
    self.chunk_lengths = chunk_lengths
    self.misses_per_chunk = misses_per_chunk
    assert(len(self.chunk_lengths) == len(self.misses_per_chunk))

  def __getitem__(self, item):
    """This make Map sliceable!"""

    if isinstance(item, slice):
      chunk_lengths = self.chunk_lengths[item]
      misses_per_chunk = self.misses_per_chunk[item]
      return NullModelMap(chunk_lengths, misses_per_chunk)

    raise TypeError("Expected slice!")

  def __len__(self):
    return len(self.chunk_lengths)

  def __str__(self):
    s = "<NullModelMap with {num_chunks} chunks at 0x{mem:x}>"
    return s.format(num_chunks = len(self),
                    mem = id(self))


###################################################################
class NullModelSimulator(object):

  def __init__(self, frags, controls):
    """Create a simulator for a reference genome with restriction fragments
    given by frags.

    The controls object is an instance of NullModelControls which controls the
    simulation. 
    """
    self.frags = np.array(frags)
    self.controls = controls
    self.N  = len(self.frags)
    self.num_reference_frags = len(self.frags)

    self.query_miss_probability = self.controls.query_miss_probability
    self.query_hit_probability = 1.0 - self.query_miss_probability
    self.log_query_miss_probability = np.log(self.query_miss_probability)
    self.log_query_hit_probability = np.log(self.query_hit_probability)

  def generate_random_reference(self):
    """Generate a random reference in-silico map where we sample from replacement
    from the catalogue of reference fragments.
    """
    inds = np.random.randint(0, self.num_reference_frags, self.num_reference_frags)
    frags = self.frags[inds]
    return frags

  def generate_random_references(self, N):
    frags = (self.generate_random_reference() for i in xrange(N))
    maps = [SOMAMap(frags = f, mapId = '') for f in frags]
    for i,m in enumerate(maps):
      m.mapId = 'random_map_%i'%i
    return maps

  def _simulate_one_map(self, miss_probability, num_chunks):
    """Return a NullModelMap populated with chunks.
       Sample fragments from self.frags using parameters in self.controls.
       Sampling includes missed sites (but no sizing error)
    """

    # Generate a random pattern from frags with num_chunks
    hit_prob = 1 - miss_probability
    g = np.random.geometric(hit_prob, num_chunks)
    randint = np.random.randint
    sum = np.sum
    frags = self.frags
    N = len(frags)

    inds = (randint(0, N, num_frags) for num_frags in g)
    chunk_lengths = np.fromiter((sum(frags[i]) for i in inds), np.int)
    misses_per_chunk = g - 1

    return NullModelMap(chunk_lengths, misses_per_chunk)

  def simulate_query_map(self, num_chunks = None):
    miss_prob = self.controls.query_miss_probability
    if not num_chunks:
      num_chunks = self.controls.max_chunks
    return self._simulate_one_map(miss_prob, num_chunks)

  def simulate_ref_map(self, num_chunks = None):
    miss_prob = 0.0
    if not num_chunks:
      num_chunks = self.controls.max_chunks
    return self._simulate_one_map(miss_prob, num_chunks)

  def simulate_query_maps(self, n):
    self.full_query_maps = [ self.simulate_query_map() for i in range(n) ]

  def simulate_ref_maps(self, n):
    self.full_ref_maps = [ self.simulate_ref_map() for i in range(n)]

  def run(self):
    
    logging.info("Generating query maps")
    self.simulate_query_maps(self.controls.N)
    
    logging.info("Generating reference maps")
    self.simulate_ref_maps(self.controls.N)


  def generate_random_frags(self, N):
    """
    Generate random restriction fragments through simulation, using the 
    nick probability in self.controls.
    """
    miss_prob = self.controls.query_miss_probability
    hit_prob = 1.0 - miss_prob
    frags = self.frags
    _generate_random_frags = self._generate_random_frags

    # To avoid too much sampling at once, only sample MAX_COLS
    # at once.
    MAX_COLS = 1000000
    num_left = N

    n = min(num_left, MAX_COLS)
    ret = _generate_random_frags(frags, hit_prob, n)
    num_left = num_left - n

    while num_left > 0:
      n = min(num_left, MAX_COLS)
      res = _generate_random_frags(frags, hit_prob,  n)
      ret.misses = np.concatenate((ret.misses, res.misses))
      ret.sizes = np.concatenate((ret.sizes, res.sizes))
      num_left = num_left - n

    return ret

  @staticmethod
  def _generate_random_frags(frags, hit_prob, N):
    """
    Generate random fragments by generating patterns by sampling with replacement from frags,
    with prob. of nick hit_prob.

    Return N N-map restriction fragments generated by this process. 
    """
    num_misses = np.random.geometric(hit_prob, N) - 1
    nrow = np.max(num_misses) + 1
    ncol = N
    inds = np.random.randint(0, len(frags), nrow*ncol) # Sampling way more than we need.
    frag_mat = frags[inds].reshape(nrow, ncol)
    chunk_lengths = np.cumsum(frag_mat, axis = 0)
    ret = ReturnObj()
    ret.sizes = chunk_lengths[num_misses, np.arange(ncol)]
    ret.misses = num_misses
    return ret

  @staticmethod
  def generate_all_random_frags(frags, N, max_frags):
    """
    Generate random fragments by generating patterns by sampling max_misses + 1 consecutive
    fragments with replacement from frags.
    Return all of the cumulative sums of N such samples.

    Returns an ReturnObject where sizes is an (max_misses + 1) x (N) numpy array
    where the (i,j) cells is the jth sample with i interior misses.
    """
    nrow = max_frags
    ncol = N
    ncells = nrow*ncol
    inds = np.random.randint(0, len(frags), ncells)
    frag_mat = frags[inds].reshape(nrow, ncol)
    chunk_lengths = np.cumsum(frag_mat, axis = 0)
    del inds
    del frag_mat

    # Sort the rows in ascending order.
    chunk_lengths.sort(axis = 1)

    ret = ReturnObj()
    ret.sizes = chunk_lengths
    return ret

  def simulate_null_distribution(self):
    """
    Simulate random N-map restriction fragments from the reference fragment distribution.

    Save a csv with the distribution of interior unmatched sites by chunk size.
    """

    sim_results = self.generate_random_frags(self.controls.num_frags)

    o = np.argsort(sim_results.sizes)
    sim_results.sizes = sim_results.sizes[o]
    sim_results.misses = sim_results.misses[o]
    max_misses = np.max(sim_results.misses)

    start = 0
    end = 500000 + 0.1
    by = 1000
    rel_error = self.controls.rel_error
    min_error = self.controls.min_error
    test_points = np.arange(start, end, by)
    num_test_points = len(test_points)

    # Perform binary search to determine the indices of the lower bound and upper bound
    # compatible with each fragment size.
    lb = np.maximum(test_points - np.maximum(rel_error*test_points, min_error), 0)
    ub = test_points + np.maximum(rel_error*test_points, min_error)
    lbi = np.searchsorted(sim_results.sizes, lb, side = 'left')
    ubi = np.searchsorted(sim_results.sizes, ub, side = 'right')

    misses_at_point = [sim_results.misses[lbi[i]:ubi[i]] for i in xrange(num_test_points)]
    misses_at_point_bincount = np.array([np.bincount(m, minlength = max_misses + 1) for m in misses_at_point])
    nr,nc = misses_at_point_bincount.shape

    num_in_bin = np.array([ubi[i] - lbi[i] for i in xrange(num_test_points)])

    # Make the probability interpolator for sampling compatible chunk of a given size
    frac_in_bin = num_in_bin.astype(float) / self.controls.num_frags

    self.frac_in_bin = frac_in_bin # save for debugging
    self.test_points = test_points # save for debugging
    self.frac_compatible_interpolator = interp1d(test_points, frac_in_bin, bounds_error = False)

    # normalize the bincount to make fractions (i.e. empirical probabilities)
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

    df.to_csv(self.controls.sim_results_csv, foat_format = "%.6f", index = False)
    self.null_distribution = df

    return df

  def simulate_null_distribution_ver2_helper(self, lb, ub, N, max_n_per_iter = 100000):
    """
    Helper function to generate counts of how many chunks fall in the sizes
    given by lowerbounds lb and upperbounds ub.

    It will draw N samples in increments of max_n_per_iter.
    """
    nrow = self.controls.max_misses + 1
    ncol = len(lb)

    num_in_bin = np.zeros((nrow, ncol))

    num_samples = 0
    while num_samples < N:
      n_to_do = min(max_n_per_iter, N - num_samples)
      sim_results = self.generate_all_random_frags(self.frags,
          n_to_do,
          self.controls.max_misses + 1)
      simulated_chunk_sizes = sim_results.sizes
      lbi = np.array([np.searchsorted(row, lb, side = 'left') for row in simulated_chunk_sizes])
      ubi = np.array([np.searchsorted(row, ub, side = 'right') for row in simulated_chunk_sizes])
      nib = ubi - lbi
      num_in_bin = num_in_bin + nib
      num_samples += n_to_do

    return num_in_bin


  def simulate_null_distribution_ver2(self):
    """
    Simulate random N-map restriction fragments from the reference fragment distribution.
    Random N-map restriction fragments are built by concatenating 0 to self.controls.max_misses + 1 fragments
    together.
    Save a csv with the distribution of interior unmatched sites by chunk size.
    """

    # sim_results = self.generate_all_random_frags(self.frags,
    #     self.controls.num_frags,
    #     self.controls.max_misses)

    # simulated_chunk_sizes = sim_results.sizes
    # nrow, ncol = sim_results.sizes.shape
    # ncell = nrow*ncol

    # del sim_results

    start = 0
    end = 500000 + 0.1
    by = 1000
    rel_error = self.controls.rel_error
    min_error = self.controls.min_error
    test_points = np.arange(start, end, by)
    num_test_points = len(test_points)

    # Perform binary search to determine the indices of the lower bound and upper bound
    # compatible with each fragment size.
    lb = np.maximum(test_points - np.maximum(rel_error*test_points, min_error), 0)
    ub = test_points + np.maximum(rel_error*test_points, min_error)

    # lbi = np.array([np.searchsorted(row, lb, side = 'left') for row in simulated_chunk_sizes])
    # ubi = np.array([np.searchsorted(row, ub, side = 'right') for row in simulated_chunk_sizes])
    # num_in_bin = ubi - lbi
    num_in_bin = self.simulate_null_distribution_ver2_helper(lb, ub, self.controls.num_frags)

    # Make the probability interpolator for sampling compatible chunk of a given size
    frac_in_bin = num_in_bin.astype(float) / self.controls.num_frags

    self.num_in_bin = num_in_bin #save for debugging
    self.frac_in_bin = frac_in_bin # save for debugging
    self.test_points = test_points # save for debugging
    self.prob_compatible_by_frag_count_interpolator = interp1d(test_points, frac_in_bin, axis = 1, bounds_error = False)

    # Create a pandas data frame summarizing the results of the simulation
    nrow, ncol = frac_in_bin.shape
    ncell = nrow*ncol
    df_data = {'midpt' : np.tile(test_points, nrow),
               'num_frag' : np.repeat(np.arange(1, nrow +1), ncol),
               'prob' : frac_in_bin.reshape(ncell,) }

    df = pandas.DataFrame(data = df_data)
    df.to_csv(self.controls.sim_results_csv, foat_format = "%.6f", index = False)
    self.null_distribution = df

    return df



  def make_interpolators(self):
    """Make interpolators for computing number of unmatched sites in 
    compatible random alignments
    """
    nd = self.null_distribution
    midpt = nd['midpt']
    frac_cols = [col for col in nd.columns if 'frac' in col]
    num_unmatched = [int(col[4:]) for col in frac_cols]
    interpolators = {}
    for i in range(len(frac_cols)):
      unmatched = num_unmatched[i]
      col_name = frac_cols[i]
      interpolators[unmatched] = interp1d(midpt, nd[col_name], assume_sorted = True, bounds_error = False)
    self.interpolators = interpolators
    return self.interpolators

  def interpolate_compatible_probs(self, frags):
    return self.frac_compatible_interpolator(frags)
  
  def interpolate_unmatched_probs(self, frags):
    """Use self.interpolators to interpolate the 
    distribution of unmatched sites for an array of fragment sizes frags.
    """
    frags = np.array(frags)
    df = pandas.DataFrame()
    interpolators = self.interpolators
    num_interpolators = len(interpolators)
    data = {i:interpolators[i](frags) for i in xrange(num_interpolators)}
    df = pandas.DataFrame(data = data)
    # for i in range(num_interpolators):
    #   df[i] = interpolators[i](frags)
    df.index = frags
    return df

  def simulate_miss_rates(self, frags, N = 5000):
    """Simulate the distribution of unmatched sites for random compatible alignments
    to the sequence of fragment sizes given by N. Use N random samples.
    Return the cumulative distribution of unmatched sites.
    """
    
    unmatched_null_distribution = self.interpolate_unmatched_probs(frags)
    num_frags, max_misses = unmatched_null_distribution.shape
    rvs = [rv_discrete(values = (np.arange(max_misses), p)) for k, p in unmatched_null_distribution.iterrows()]

    # Make random samples for number of misses for each fragment
    mat = np.zeros((num_frags, N))
    for row, rv in enumerate(rvs):
      mat[row,:] = rv.rvs(size=N)

    # Take column sums for the number of misses
    miss_counts = mat.sum(axis=0).astype(int)
    miss_count_binned = cumsum(bincount(miss_counts)).astype(float)/N
    return miss_count_binned

  def null_probability(self, frag_sizes, interior_unmatched):
    """Compute the null probability of generating a random sequence
    of fragments of sizes frag_sizes with the given interior unmatched sites.
    """
    prob_selected = self.interpolate_compatible_probs(frag_sizes)
    prob_unmatched = self.interpolate_unmatched_probs(frag_sizes).as_matrix()
    unmatched_probs = prob_unmatched[range(len(frag_sizes)), interior_unmatched.astype(int)]
    return (prob_selected, unmatched_probs)

  def null_probability_ver2(self, frag_sizes, interior_unmatched):
    """Compute the null probability of generating a random sequence of chunks compatible
    with the given chunk sizes and the interior number of unmatched sites in a each chunk.
    """
    prob_selected = self.prob_compatible_by_frag_count_interpolator(frag_sizes)
    return np.prod(prob_selected[interior_unmatched, np.arange(len(frag_sizes))])

  def log_null_probability_ver2(self, frag_sizes, interior_unmatched):
    """Compute the null probability of generating a random sequence of chunks compatible
    with the given chunk sizes and the interior number of unmatched sites in a each chunk.
    """
    prob_selected = self.prob_compatible_by_frag_count_interpolator(frag_sizes)
    return np.sum(np.log(prob_selected[interior_unmatched, np.arange(len(frag_sizes))]))

  def null_probability_generated(self, frag_sizes, interior_unmatched):
    """Compute the log probability of generating a compatible alignment with these miss characteristics
    to a random genome"""
    prob_chunk_selected, unmatched_probs = self.null_probability(frag_sizes, interior_unmatched)
    prob_generated = np.prod(prob_chunk_selected*unmatched_probs)
    return prob_generated

  def log_pattern_probability(self, interior_unmatched):
    """Compute the probability of producing the restriction pattern.
    It's not clear if we should include the bounding matched sites or not.
    For a probability odds ratio, this will end up being a constant factor shift.
    
    For now, include the boundaries.
    """
    num_unmatched = np.sum(interior_unmatched)
    num_matched = 1 + len(interior_unmatched)
    return num_matched * self.log_query_hit_probability + num_unmatched * self.log_query_miss_probability

  def log_likelihood_ratio(self, frag_sizes, interior_unmatched):
    """Return the log likelihood_ratio that the alignment has the observed restriction pattern
    to the probability of generating the sequence of chunks randomly.

    Note: The probability of randomly generating the sequence of chunks as they appear in the alignment will
    be small, but this is not the same probability as asking "what is the probability that a random genome has this 
    chunk sequence in one or more locations". The latter probability is larger, because allows for the alignment chunk
    pattern to appear anywhere in either the forward or reverse orientation.
    """
    L_HA = self.log_pattern_probability(interior_unmatched)
    L_H0 = log(self.null_probability_generated(frag_sizes, interior_unmatched))
    return L_HA - L_H0

  def assign_likelihoods(self, aln):
    """
    Assign likelihoods to an alignment object aln.
    The alignment object should have misses and frag_lengths attributes which describe
    the number of misses and total length of each chunk.
    """
    aln.log_HA = self.log_pattern_probability(aln.misses)
    aln.p_H0 = self.null_probability_ver2(aln.frag_lengths, aln.misses)
    aln.log_H0 = self.log_null_probability_ver2(aln.frag_lengths, aln.misses)
    aln.p_H0_alt = np.exp(aln.log_H0)

    aln.log_likelihood_ratio = aln.log_HA - aln.log_H0
    aln.E_H0 = aln.p_H0*2*self.num_reference_frags # Expected number of times an alignment compatible to this one could appear in a random genome.


def all_scores_forward(scorer, query_map, ref_map):
  lq = len(query_map)
  lr = len(ref_map)
  last_pos = lr - lq
  sf = scorer.score_at_position
  scores = [sf(query_map, ref_map, p) for p in xrange(last_pos)]
  return scores







