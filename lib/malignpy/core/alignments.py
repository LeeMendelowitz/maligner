"""
Utilities for working with alignment objects
"""
from collections import defaultdict

# Sift through the alignments and select the best non-overlapping alignments.
def sift_alignments(alns):
  """
  Sift through the alignments and select non-overlapping aligments.

  Give preference to alignments in the order in which they appear.
  """

  ref_to_covered_intervals = defaultdict(list)

  def make_key(aln):
    return (aln.ref_id, aln.ref_is_forward)

  def make_value(aln):
    firstRefChunk = aln.matched_chunks[0].ref_chunk
    lastRefChunk = aln.matched_chunks[-1].ref_chunk
    return (firstRefChunk.start, lastRefChunk.end)

  def overlaps(i1, i2):
    s1, e1 = i1
    s2, e2 = i2
    return not ((e1 <= s2) or (e2 <= s1))

  def overlaps_any(i, ilist):
    return any(overlaps(i,iother) for iother in ilist)

  sifted = []
  for aln in alns:
    k = make_key(aln)
    interval = make_value(aln)
    if overlaps_any(interval, ref_to_covered_intervals[k]):
      continue
    ref_to_covered_intervals[k].append(interval)
    sifted.append(aln)

  return sifted
