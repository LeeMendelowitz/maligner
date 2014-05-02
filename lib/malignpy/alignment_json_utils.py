"""
Utilities for working with alignment json files.
"""

import json
from json import loads

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

    line = fin.readline().strip().strip(',')
    count += 1

    if not line:
      break

    pos = fin.tell()

    try:
      d = loads(line)
    except Exception as e:
      print 'Could not decode line %i: %s'%( count, str(e) )
      continue

    query_id = d['query_id']
    if query_id != last_query_id:
      fout.write('%s,%i\n'%(query_id, pos))
    last_query_id = query_id


@wrap_file_function('r')
def load_index(f):
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
    except Exception as e:
      print 'Could not decode line: %s'%( str(e) )
      break

    if d['query_id'] == query_id:
      alignments.append(d)
    else:
      break

  return alignments


def compute_rescaled_score(aln):
  return sum(aln['rescaled_score'].values())


