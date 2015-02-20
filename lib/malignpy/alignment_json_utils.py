"""
Utilities for working with alignment json files.

Functions for indexing and parsing alignments files.
"""

import json
from json import loads
from collections import OrderedDict

from common import wrap_file_function

@wrap_file_function('r', 'w')
def index_alignment_json(fin, fout):
  """
  Index an alignment json.
  """

  index = {}
  last_query_id = None
  loads = json.loads
  count = 0
  while True:

    pos = fin.tell()
    line = fin.readline().strip().strip(',')
    count += 1

    if not line:
      break


    try:
      d = loads(line)
    except ValueError as e:
      print 'Could not decode line %i: %s'%( count, str(e) )
      continue

    query_id = d['query_id']
    if query_id != last_query_id:
      fout.write('%s,%i\n'%(query_id, pos))
    last_query_id = query_id


@wrap_file_function('r')
def load_index(f):
  """
  Load an alignment index file.
  """
  recs = (l.split(',') for l in f)
  d = { r[0] : int(r[1]) for r in recs }
  return d

@wrap_file_function('r')
def get_alignments(fp, query_id, index):
  fp.seek(index[query_id])
  alignments = []
  while True:

    line = fp.readline().strip().strip(',')
    try:
      d = loads(line)
    except ValueError as e:
      print 'Could not decode line: %s'%( str(e) )
      break

    if d['query_id'] == query_id:
      alignments.append(d)
    else:
      break

  return alignments

@wrap_file_function('r')
def get_first_alignment(fp, query_id, index):
  """
  Get the first alignment for query_id in alignments json file fp.
  """
  fp.seek(index[query_id])
  while True:

    line = fp.readline().strip().strip(',')
    try:
      d = loads(line)
      return d
    except ValueError as e:
      print 'Could not decode line: %s'%( str(e) )
      break

  return None


def compute_rescaled_score(aln):
  """
  Compute the total rescaled score value.
  """
  return sum(aln['rescaled_score'].values())

def parse_alignment(a):
  """
  Summarize an alignment with a flat dictionary
  """

  d = OrderedDict()

  score = a['score']
  rescaled_score = a['rescaled_score']
  total_score = sum(score.values())
  total_score_rescaled = sum(rescaled_score.values())


  keys = ['query_id',
          'ref_id',
          'ref_is_forward',
          'total_score',
          'query_num_frags',
          'query_num_sites']

  def update(keys, data = d):
    for k in keys:
      d[k] = a[k]
  
  update(keys)

  keys = ['query_miss_rate',
          'ref_miss_rate',
          'total_miss_rate',
          'query_misses',
          'ref_misses',
          'max_chunk_sizing_score',
          'query_scaling_factor',
          'hit_log_likelihood',
          'miss_log_likelihood',
          'chunk_log_likelihood',
          'total_log_likelihood']
  
  update(keys)

  d['query_miss_score'] = score['query_miss_score']
  d['ref_miss_score'] = score['ref_miss_score']
  d['sizing_score'] = score['sizing_score']
  d['sizing_score_rescaled'] = rescaled_score['sizing_score']
  d['total_score_rescaled'] = total_score_rescaled

  return d

@wrap_file_function('r')
def parse_alignments(f):
  """
  Read an entire alignments file. Flatten the alignment records.
  Return a list of alignments.
  """
  try:
    alignments = []
    for i, l in enumerate(f):
        try:
            data = json.loads(l.strip())
            alignments.append(data)
        except ValueError as e:
            print 'Could not parse line %i:\n\t%s'%(i, l)
            continue

  except ValueError:
    stderr("Could not parse json from file %s\n"%f)
    return []
  alignments = [parse_alignment(a) for a in alignments]
  return alignments

@wrap_file_function('r')
def read_all_alignments(f):
  """
  Read an entire alignments file. Do not flatten the alignment records.
  Return a list of alignments.
  """
  return [r for r in json_iter(f)]

def parse_alignments_iter(f):
  """
  Parse and flatten alignment records
  """
  for i, l in enumerate(f):

      try:
          data = json.loads(l.strip())
          yield parse_alignment(data)
      except ValueError as e:
          print 'Could not parse line %i'%i
          continue


def json_iter(f):
  """
  Read alignment records as dicts.
  """
  for i, l in enumerate(f):

      try:
          data = json.loads(l.strip())
          yield data
      except ValueError as e:
          print 'Could not parse line %i'%i
          continue




