# Utilities to parse alignment files from Schwartz Lab SOMA
# *and* to write alignments in Schwartz lab SOMA format
#
# We use functions and classes defined here both to convert
# from maligner to SOMA and from SOMA to maligner.
#
# This module could be better organized!
# Perhaps separate modules for parsing and writing?
#################################################

import lxml.etree as ET
from lxml.etree import Element, ElementTree
from copy import copy, deepcopy
from ..common import wrap_file_function
from ..core.maligner_dp_alignments import (Alignment as MAlignment, Chunk as MChunk)


##################################################################
# Helper classes and functions for parsing generic XML nodes
def smart_bool(s):
  """Convert an str to a bool, flexibly based on first character"""
  if not s:
    return False
  s = s.upper()
  if s[0] == 'T':
    return True
  return False

class SimpleField(object):
  def __init__(self, converter = str):
    self.converter = converter

  def __call__(self, n, data):
    data[n.tag] = self.converter(n.text)

class SmartBool(SimpleField):
  def __init__(self):
    SimpleField.__init__(self, smart_bool)

class ListNodeField(object):

  def __init__(self, converter = str, data_key = None):
    self.converter = converter
    self.data_key = data_key

  def __call__(self, n, data):
    data_key = self.data_key if self.data_key else n.tag
    cur_data = data.get(data_key, [])
    cur_data.append(self.converter(n))
    data[data_key] = cur_data

def MakeTextElement(tag, text = None):
  elem = Element(tag)
  if text is not None:
    elem.text = text
  return elem

class Blob(object):
  pass

################################################################
# Helper functions for parsing alignment nodes

def _parse_reference_map_field(n, data):
  data[n.tag] = n.find('name').text

def _parse_aligned_map(n, data):

  orientation = n.find('orientation').text
  orientation = "F" if orientation == "N" else "R"

  data.update({
    n.tag: n.find('name').text,
    'orientation': orientation
  })

def _parse_f(n):
  d = dict((c.tag, int(c.text)) for c in n)
  return (d['i'], d['l'], d['r'])

################################################################
# Classes representing parsed nodes 

class NodeClass(object):

  # Register parsers for child nodes of this node by
  # setting _node_parsers in subclass of NodeClass
  _node_parsers = {}

  @classmethod
  def from_node(cls, n, parsers = None):
    """Build a new instance of the class 
    from the node n. If supplied, use the parsers 
    provided. Otherwise, use the class's _node_parsers.
    """
    ret = cls()
    ret.node = n
    ret._from_node(n, parsers)
    return ret

  def _from_node(self, n, parsers = None):
    """Set data extracted from node n.
    If supplied, use the parsers provided. Otherwise
    use the instance's _node_parsers.
    """
    data = {}

    if parsers is None:
      parsers = self._node_parsers

    # Parsing all defined parsers
    for c in n:
      if c.tag in parsers:
        self._node_parsers[c.tag](c, data)
      else:
        data[c.tag] = c.text

    for k, v in data.iteritems():
      setattr(self, k, v)


class RestrictionMap(NodeClass):
  """
  Class to initialze from a node from SOMA XML file
  """

  _node_parsers = {
    'map_block' : SimpleField(lambda s: [float(b) for b in s.split()]),
    'circular' : SmartBool()
  }

  @classmethod
  def from_node(cls, n):

    self = cls()
    self.node = n
    self._from_node(n) # Use the _node_parsers

    # Add additional attributes
    self.frags = [int(1000*f) for f in self.map_block]

    return self

  @property
  def nfrags(self):
    return len(self.frags)

  @property
  def length(self):
    """Length in bp"""
    return sum(self.frags)

  def write_as_maligner(self, f):
    s = '{0}\t{1}\t{2}\t{3}'.format(self.name, self.length, self.nfrags,
      '\t'.join(str(frag) for frag in self.frags))
    f.write(s + '\n')


class Alignment(NodeClass):

  # These fields are populated by cls.from_node():
  uuid = None
  reference_map = None
  aligned_map = None
  soma_score = None
  count = None
  f = []

  _node_parsers = {
    'uuid' : SimpleField(),
    'reference_map': _parse_reference_map_field,
    'aligned_map':  _parse_aligned_map,
    'soma_score' : SimpleField(float),
    'count' : SimpleField(int),
    'f' : ListNodeField(_parse_f)
  }

  def build_chunks(self):
    # Convert SOMA chunks indices to maligned Chunk indices
    self.chunks = [MatchedChunk(*i) for i in self.f]
    self._chunks_merged = merge_chunks_list(self.chunks)

    self.chunks[0].is_boundary = True
    self.chunks[-1].is_boundary = True

    self._chunks_merged[0].is_boundary = True
    self._chunks_merged[-1].is_boundary = True

  # Add map data to the alignment so we can compute other
  # alignment characteristics.
  def add_map_data(self, query_map, ref_map):
    qm, rm = query_map, ref_map

    qfrags = qm.map_block

    # Orient query fragments in the reverse direction to
    # assist with determining chunk lengths
    if self.orientation == 'R':
      qfrags = qfrags[::-1]

    # Assign query lengths and ref lengths to the chunks,
    # convert to kb
    for c in self.chunks:
      c.ql = int(1000*sum(qfrags[c.q_s:c.q_e]))
      c.rl = int(1000*sum(rm.map_block[c.r_s:c.r_e]))

    for c in self._chunks_merged:
      c.ql = int(1000*sum(qfrags[c.q_s:c.q_e]))
      c.rl = int(1000*sum(rm.map_block[c.r_s:c.r_e]))

    self.query_map = self.aligned_map

    self.num_query_frags = len(qfrags)

    cpairs = zip(self._chunks_merged[:-1], self._chunks_merged[1:])
    self.has_query_gap = not all(cl.q_e == cr.q_s for cl,cr in cpairs)
    self.has_ref_gap = not all(cl.r_e  == cr.r_s for cl,cr in cpairs)

    qs = self._chunks_merged[0].q_s
    rs = self._chunks_merged[0].r_s
    qe = self._chunks_merged[-1].q_e
    re = self._chunks_merged[-1].r_e
    self.num_matched_chunks = len(self._chunks_merged)
    self.num_frags_aligned = qe - qs
    self.query_start = qs
    self.query_end = qe
    self.ref_start = rs
    self.ref_end = re
    self.num_query_frags_aligned = qe - qs
    self.trim_left = qs
    self.trim_right = self.num_query_frags - qe
    self.query_length = sum(qfrags)

    interior_chunks = [c for c in self._chunks_merged if not c.is_boundary]

    self.query_length_aligned = sum(c.ql for c in self._chunks_merged)
    self.query_length_aligned_interior = sum(c.ql for c in interior_chunks)
    self.ref_length_aligned = sum(c.rl for c in self._chunks_merged)
    self.ref_length_aligned_interior = sum(c.rl for c in interior_chunks)

    self.query_scaling_factor = float(self.ref_length_aligned_interior)/float(self.query_length_aligned_interior)

    self.frac_length_aligned = float(self.query_length_aligned)/self.query_length
    self.frac_frags_aligned = float(qe - qs)/self.num_query_frags

    self.query_misses = sum(c.query_misses for c in self._chunks_merged)
    self.ref_misses = sum(c.ref_misses for c in self._chunks_merged)
    self.matched_sites = len(self._chunks_merged)+1
    self.query_miss_rate = float(self.query_misses)/(self.query_misses + self.matched_sites)
    self.ref_miss_rate = float(self.ref_misses)/(self.ref_misses + self.matched_sites)

    # Orient chunks for maligner
    self.chunks_maligner = deepcopy(self._chunks_merged)
    if self.orientation == "R":
      nqf = self.num_query_frags
      for cm, cs in zip(self.chunks_maligner, self._chunks_merged):
        cm.q_s = nqf - cs.q_e
        cm.q_e = nqf - cs.q_s

    self.chunks_maligner_rescaled = deepcopy(self.chunks_maligner)

    # Rescale the chunks for maligner, convert to integer
    for c in self.chunks_maligner_rescaled:
      c.rl = int(c.rl)
      c.ql = int(c.ql * self.query_scaling_factor)

  def compute_maligner_sizing_error(self, sd_rate, min_sd ):

    for c in self.chunks_maligner_rescaled:
      csd = max(sd_rate * c.rl, min_sd)
      delta = c.ql - c.rl
      chi2 = (delta/csd)**2
      c.chi2 = chi2
      c.sd = csd

    for c in self.chunks_maligner:
      csd = max(sd_rate * c.rl, min_sd)
      delta = c.ql - c.rl
      chi2 = (delta/csd)**2
      c.chi2 = chi2
      c.sd = csd

    self.sizing_score = sum(c.chi2 if not c.is_boundary else 0 for c in self.chunks_maligner)
    self.sizing_score_rescaled = sum(c.chi2 if not c.is_boundary else 0 for c in self.chunks_maligner_rescaled)

  def compute_scores(self, query_miss_penalty, ref_miss_penalty, sd_rate = 0.05, min_sd = 1000.0):
    miss_score = self.query_misses*query_miss_penalty + self.ref_misses*ref_miss_penalty
    self.compute_maligner_sizing_error(sd_rate, min_sd)
    self.total_score = self.sizing_score + miss_score
    self.total_rescaled_score = self.sizing_score_rescaled + miss_score

  @property
  def chunk_string(self):
    # Write maligner chunk string.
    def c_to_str(c):
      fields = [c.q_s, c.q_e, c.ql, c.r_s, c.r_e, c.rl]
      return ','.join([str(f) for f in fields])
    return ';'.join(c_to_str(c) for c in self.chunks_maligner_rescaled)

##################################################################
# Methods for converting to a SOMA Alignments File

def make_experiment():
  pass


def restriction_map_to_node(name, length, frags, is_consensus = True, enzymes="None", circular=False):
  """Initialize from name, frags (bp)"""

  frags = list(frags)

  sub_elements = [

    MakeTextElement("name", text = name),
    MakeTextElement("type", text = "consensus" if is_consensus else "opmap"),
    MakeTextElement("enzymes", text = enzymes),
    MakeTextElement("circular", text =str(circular).lower()),
    MakeTextElement("num_frags", text = str(len(frags))),
    MakeTextElement("map_block", text = " ".join("%.3f"%(frag/1000.0) for frag in frags))

  ]

  element = MakeTextElement("restriction_map")
  element.extend(sub_elements)

  return element

def alignment_to_node(uuid, aln, query_frags, ref_frags, soma_score = 0.0):
  """Return an ETree Element for writing to Alignment File.

  aln should be an instance of maligner_dp_alignments.Alignment

  """

  # <map_alignment>
  #   <uuid>3ca38422-b6e7-43e8-a925-ef36d475a2d7</uuid>
  #   <reference_map>
  #     <name>chr1</name>
  #   </reference_map>
  #   <aligned_map>
  #     <name>2248849_0_44</name>
  #     <orientation>R</orientation>
  #   </aligned_map>
  #   <soma_score>3.3207884</soma_score>
  #   <count>18</count>
  #   <f><i>1</i><l>26398</l><r>26398</r></f>
  #   <f><i>2</i><l>26399</l><r>26399</r></f>
  #   <f><i>3</i><l>26400</l><r>26401</r></f>
  #   <f><i>4</i><l>26402</l><r>26402</r></f>
  #   <f><i>5</i><l>26403</l><r>26403</r></f>
  #   <f><i>6</i><l>26404</l><r>26405</r></f>
  #   <f><i>7</i><l>26406</l><r>26406</r></f>
  #   <f><i>8</i><l>26406</l><r>26406</r></f>
  #   <f><i>9</i><l>26407</l><r>26407</r></f>
  #   <f><i>10</i><l>26408</l><r>26408</r></f>
  #   <f><i>11</i><l>26409</l><r>26409</r></f>
  #   <f><i>12</i><l>26410</l><r>26410</r></f>
  #   <f><i>13</i><l>26411</l><r>26411</r></f>
  #   <f><i>14</i><l>26412</l><r>26413</r></f>
  #   <f><i>15</i><l>26414</l><r>26414</r></f>
  #   <f><i>16</i><l>26415</l><r>26415</r></f>
  #   <f><i>17</i><l>26416</l><r>26419</r></f>
  #   <f><i>18</i><l>26420</l><r>26421</r></f>
  # </map_alignment>


  ref_element = MakeTextElement("reference_map")
  ref_element.append(MakeTextElement("name", text = aln.ref_map))

  aligned_map_element = MakeTextElement("aligned_map")
  aligned_map_element.append(MakeTextElement("name", text = aln.query_map))
  aligned_map_element.append(MakeTextElement("orientation", text = "N" if aln.is_forward=="F" else "R"))

  num_chunks = len(aln.chunks)
  num_query_frags = aln.query_end - aln.query_start

  soma_frags = convert_alignment_chunks(aln, query_frags, ref_frags)
  soma_frag_elements = [f.to_element() for f in soma_frags]

  sub_elements = [

    MakeTextElement("uuid", text=uuid),
    ref_element,
    aligned_map_element,
    MakeTextElement("soma_score", text=str(soma_score)),
    MakeTextElement("count", text = str(num_query_frags))

  ]

  sub_elements.extend(soma_frag_elements)

  element = MakeTextElement("map_alignment")
  element.extend(sub_elements)

  return element


###################################################################################
# Functions for working with XML files and ET root nodes.

class ReturnValue(object):
  pass

def get_root(f):
  tree = ET.parse(open(f))
  return tree.getroot()

def get_maps(tr):
  nodes = tr.iter('restriction_map')
  return [RestrictionMap.from_node(n) for n in nodes]

def get_consensus_maps(r):
  maps = (RestrictionMap(m) for m in r.iter('restriction_map'))
  return [m for m in maps if m.type == 'consensus']

@wrap_file_function('r')
def get_consensus_maps_from_file(aln_file):
  return get_consensus_maps(get_root(aln_file))

# Get consensus maps from multiple files
def get_all_consensus_maps(file_list):
  maps = []
  for f in file_list:
    cmaps = get_consensus_maps(get_root(f))
    maps.extend(cmaps)
  return maps

def get_query_maps(tr):
  """Return query maps from a SOMA alignments file"""
  all_maps = get_maps(tr)
  query_maps = [m for m in all_maps if m.type != 'consensus']
  return query_maps

def get_all_maps(tr):
  """Return all maps from a SOMA alignments file"""
  maps = get_maps(tr)
  ret = ReturnValue()
  ret.reference_maps = dict((m.name, m) for m in maps if m.type == 'consensus')
  ret.query_maps = dict((m.name, m) for m in maps if m.type != 'consensus')
  return ret

def make_alignment_from_node(n, maps):
  """Make an alignment from the xml node. 
  This means constructing an Alignment object,
  linking it with the restriction maps
  to determine various alignment attributes.
  """
  aln = Alignment.from_node(n)

  qm = maps.query_maps.get(aln.aligned_map, None)
  rm = maps.reference_maps.get(aln.reference_map, None)

  aln.build_chunks()

  aln.add_map_data(qm, rm)

  return aln

def iter_alignments(tr, maps):
  """tr: root node of the xml cElementTree
     maps: dictionary of reference and query maps, returned by get_all_maps
  """

  alns = tr.iter('map_alignment')

  return (make_alignment_from_node(n, maps) for n in alns)

def read_file(f):
  """
  Read an entire xml and return the root node. Beware of memory usage if
  """
  return ET.parse(f)







##########################################################################
# Classes for dealing with "chunks" (alignment blocks/intervals)
# There are differences between how these are represented in maligned and in Schwartz format.

class MergeException(Exception):
  pass
  
class MatchedChunk(object):
  """Helper class to assist with converting indices of Schwartz
  map_alignment to the maligner format
  """
  def __init__(self, i, l, r):
    self.i = i
    self.l = l
    self.r = r #inclusive!

    # Intervals of chunks in query and reference, with
    # ending indices exclusive (as in python)
    self.q_s = i
    self.q_e = i + 1
    self.r_s = self.l
    self.r_e = self.r + 1

    self.ql = None
    self.rl = None

    self.is_boundary = False

  def overlaps_in_reference(self, other):
    return not ((self.r < other.l) or (other.r < self.l))

  def should_merge(self, other):
    """Return True if this MatchedChunk should be merged with other on the right"""
    if not self.overlaps_in_reference(other):
      return False      
    if other.q_s == self.q_e:
      return True
    return False

  def merge(self, other):
    """Merge this MatchedChunk with a right MatchedChunk"""

    if not self.should_merge(other):
      raise MergeException("Chunks should not be merged.")

    # Assert that the Chunks are adjacent in the query?
    # Not necessarily true if there is a deletion from query]
    self.q_e = other.q_e
    self.r_s = min(self.r_s, other.r_s)
    self.r_e = max(self.r_e, other.r_e)

  def __repr__(self):
    return '({0},{1},{ql})-({2},{3},{rl})'.format(self.q_s, self.q_e,
      self.r_s, self.r_e, ql=self.ql, rl = self.rl)

  @property
  def query_misses(self):
    return self.q_e - self.q_s - 1

  @property
  def ref_misses(self):
    return self.r_e - self.r_s - 1

  @property
  def rel_error(self):
    d = self.ql - self.rl
    return d/self.rl


def merge_chunks_list(cl):
  """Merge overlapping chunks in effort to convert from SOMA format
  to maligner format. This makes a copy of the chunk list and then modifies."""

  if not cl:
    return []

  cl = [copy(c) for c in cl]

  cur_chunks = cl

  while True:
    made_change = False
    new_chunks = []
    cur_chunk = cur_chunks[0]
    for i in range(1, len(cur_chunks)):
      next_chunk = cur_chunks[i]
      if cur_chunk.should_merge(next_chunk):
        cur_chunk.merge(next_chunk)
        made_change = True
      else:
        new_chunks.append(cur_chunk)
        cur_chunk = next_chunk
        next_chunk = None
    if cur_chunk:
      new_chunks.append(cur_chunk)

    cur_chunks = new_chunks

    if not made_change:
      break

  return cur_chunks

class SomaFrag(object):

  def __init__(self, i, l, r):
    self.i = i
    self.l = l
    self.r = r
    assert(l <= r)

  def is_sane(self):
    return (self.i >= 0 and self.l <= self.r)


  def to_element(self):
    f_element = Element('f')
    f_element.append(MakeTextElement('i', text = str(self.i)))
    f_element.append(MakeTextElement('l', text = str(self.l)))
    f_element.append(MakeTextElement('r', text = str(self.r)))
    return f_element

def maligner_chunk_to_soma_frags(chunk, query_frags, ref_frags):

  frags = _maligner_chunk_to_soma_frags(chunk.qs, chunk.qe,
      chunk.rs, chunk.re, query_frags, ref_frags)

  assert(_check_soma_frags_sane(frags))

  return frags

def _maligner_chunk_to_soma_frags(qs, qe, rs, re, query_frags, ref_frags):
  """Convert the indices of a maligner chunk to one or more soma frags
  This is tricky because there are multiple cases.

  Throughput [qs,qe,rs,re] is a maligner chunk give by indices of query start, end & ref start, end.
  <i, l, r> is a soma chunk showing that fragment i is matched to reference map fragments l inclusive through r inclusive
  A. one query frag to one ref. frag
  [1,2,10,11] -> <1, 10, 10> 

  B. one query frag to muliple ref frag.
  [1,2,10,12] -> <1, 10, 11>

  C multiple query frag to multiple ref frag.
    In this case there are two possible translations:
    i) [1, 3, 10, 12] -> <1, 10, 10>, <2, 10, 11> OR 
    ii) [1, 3, 10, 12] -> <1, 10, 11> <2, 11, 11>
  
  We will pick based on where the restriction cuts are.

  Case i):
    1        2
  |-----|----------------| query
  |-------------|--------| reference
      10           11

  Case ii):
    1               2
  |--------------|--------| query
  |------|----------------| reference
      10           11

  """

  nq = qe-qs
  nr = re-rs

  assert(qe <= len(query_frags))
  assert(re <= len(ref_frags))

  assert(nq > 0)
  assert(nr > 0)

  if nq == nr == 1:

    # one to one
    return [SomaFrag(qs, rs, rs)]
  
  elif nq == 1 and nr >1:

    # one to many
    return [SomaFrag(qs, rs, re-1)]

  elif nq > 1 and nr ==1:

    return [SomaFrag(qs+i, rs, rs) for i in range(nq)]

  elif nq > 1 and nr > 1:

    # many to many. Pretty tricky.
    # go from left to right and decide on the ref. indices based on fragment size.

    q = query_frags[qs]
    r = ref_frags[rs]

    if q <= r:
      return [SomaFrag(qs, rs, rs)] + _maligner_chunk_to_soma_frags(qs+1, qe, rs, re, query_frags, ref_frags)
    else:

      # Keep adding reference chunks until reference grows too large:
      rr = rs
      last_r = re - 1
      while  rr < last_r and r < q:
        rr += 1
        r += ref_frags[rr]
      return [SomaFrag(qs, rs, rr)] + _maligner_chunk_to_soma_frags(qs+1, qe, rr, re, query_frags, ref_frags)  
  else:
    raise RuntimeError("how'd you get here?")

def _check_soma_frags_sane(frags):

  is_sane = all(f.is_sane() for f in frags)

  if len(frags) == 1:
    return is_sane

  for fl, fr in zip(frags[:-1], frags[1:]):
    is_sane = is_sane and (fr.i == fl.i + 1)
    is_sane = is_sane and (fl.l <= fr.l)
    is_sane = is_sane and ((fl.r <= fr.r))

  if not is_sane:
    import pdb; pdb.set_trace()

  return is_sane

def convert_alignment_chunks(aln, query_frags, ref_frags):
  """Convert maligner chunks to soma frags"""

  # If the alignment is reverse, we need to flip the query_frags AND
  # flip the indices in the chunks. This is so we match soma's
  # output of reverse alignments.

  chunks = aln.chunks
  num_query_frags = len(query_frags)
  if aln.is_forward == "R":

    query_frags = query_frags[::-1]

    chunks = [copy(c) for c in chunks] # Make a copy
    for c in chunks:
      c.flip_query_coords(num_query_frags) # Flip query coords

  # Check that things are sane:
  cqe_prev = chunks[0].qe
  for c in chunks[1:]:
    assert(c.qs == cqe_prev)
    cqe_prev = c.qe

  frag_lists = [maligner_chunk_to_soma_frags(c, query_frags, ref_frags) for c in chunks]
  soma_frags = [f for fl in frag_lists for f in fl]
  assert(_check_soma_frags_sane(soma_frags))
  return soma_frags



class AlignmentWriter(object):
  """Write alignment records in the maligner format"""

  # Maligner output to Alignment attr
  maligner_fields_mapping = [
    ('query_map', 'query_map'),
    ('ref_map', 'reference_map'),
    ('is_forward', 'orientation'),
    ('num_matched_chunks', 'num_matched_chunks' ),
    ('query_misses', 'query_misses'),
    ('ref_misses', 'ref_misses'),
    ('query_miss_rate', 'query_miss_rate' ),
    ('ref_miss_rate', 'ref_miss_rate' ),
    ('total_score', 'total_score'),
    ('total_rescaled_score', 'total_rescaled_score'),
    ('sizing_score', "sizing_score"),
    ('sizing_score_rescaled', "sizing_score_rescaled" ),
    ('query_scaling_factor', 'query_scaling_factor'),
    ('chunk_string', "chunk_string")
  ]

  maligner_fields = [m for m,a in maligner_fields_mapping]

  def __init__(self):
    pass

  def _formatter(self, val):
    if isinstance(val, float):
      return '%.6f'%val
    return str(val)

  def _gen_maligner_fields(self, aln):
    for m, a in self.maligner_fields_mapping:
      if a is None:
        yield "NA"
        continue
      yield getattr(aln, a, "NA")

  def to_maligner_str(self, aln):
    """Output a string for maligner output file"""
    f = self._formatter
    fields = [f(v) for v in self._gen_maligner_fields(aln)]
    return '\t'.join(fields)

  def write_maligner_alignment(self, f, aln):
    """Output a line to the maligner output file"""
    f.write(self.to_maligner_str(aln) + '\n')

  def write_maligner_header(self, f):
    """Write maligner header"""
    f.write('\t'.join(self.maligner_fields) + '\n')





