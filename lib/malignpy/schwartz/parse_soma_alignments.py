# Parse alignment files from Schwartz Lab SOMA

import xml.etree.cElementTree as ET
from copy import copy, deepcoy

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

################################################################
# Helper functions for parsing alignment nodes

def _parse_reference_map_field(n, data):
	data[n.tag] = n.find('name').text

def _parse_aligned_map(n, data):
	data.update({n.tag: n.find('name').text,
	       'orientation': n.find('orientation').text})

def _parse_f(n):
	d = dict((c.tag, int(c.text)) for c in n)
	return (d['i'], d['l'], d['r'])

################################################################
# Classes representing parsed nodes	

class NodeClass(object):

	# Register parsers for child nodes of this node by
	# setting _node_parsers in subclass of NodeClass
	_node_parsers = {}

	def __init__(self, n):
		self.node = n

	@classmethod
	def from_node(cls, n):
		aln = cls(n)

		data = {}

		# Parsing all defined parsers
		for c in n:
			if c.tag in cls._node_parsers:
				cls._node_parsers[c.tag](c, data)
			else:
				data[c.tag] = c.text

		for k, v in data.iteritems():
			setattr(aln, k, v)

		return aln

class RestrictionMap(NodeClass):

	_node_parsers = {
		'map_block' : SimpleField(lambda s: [float(b) for b in s.split()]),
		'circular' : SmartBool()
	}

	@classmethod
	def from_node(cls, n):
		self = NodeClass.from_node(n)
		self.frags = [int(1000*f) for f in self.map_block]

  @property
  def nfrags(self):
    return len(self.frags)

  @property
  def length(self):
  	return sum(self.frags)

  def write_as_maligner(self, f):
    s = '{0}\t{1}\t{2}\t{3}'.format(self.name, self.length, self.nfrags,
      '\t'.join(str(frag) for frag in self.frags))
    f.write(s + '\n')


class Alignment(NodeClass):

	def build_chunks(self):
		self.chunks = [MatchedChunk(*i) for i in self.f]
		self.chunks_merged = merge_chunks_list(self.chunks)

	# Add map data to the alignment so we can compute other
	# alignment characteristics.
	def add_map_data(self, query_map, ref_map):
		qm, rm = query_map, ref_map

		qfrags = qm.map_block

		# Orient query fragments in the reverse direction to
		# assist with determining chunk lengths
		if self.orientation == 'R':
			qfrags = qfrags[::-1]

		# Assign query lengths and ref lengths to the chunks
		for c in self.chunks:
			c.ql = sum(qfrags[c.q_s:c.q_e])
			c.rl = sum(rm.map_block[c.r_s:c.r_e])

		for c in self.chunks_merged:
			c.ql = sum(qfrags[c.q_s:c.q_e])
			c.rl = sum(rm.map_block[c.r_s:c.r_e])

		self.query_map = self.aligned_map

		self.num_query_frags = len(qfrags)

		cpairs = zip(self.chunks_merged[:-1], self.chunks_merged[1:])
		self.has_query_gap = not all(cl.q_e == cr.q_s for cl,cr in cpairs)
		self.has_ref_gap = not all(cl.r_e  == cr.r_s for cl,cr in cpairs)

		qs = self.chunks_merged[0].q_s
		rs = self.chunks_merged[0].r_s
		qe = self.chunks_merged[-1].q_e
		re = self.chunks_merged[-1].r_e
		self.num_frags_aligned = qe - qs
		self.query_start = qs
		self.query_end = qe
		self.ref_start = rs
		self.ref_end = re
		self.num_query_frags_aligned = qe - qs
		self.trim_left = qs
		self.trim_right = self.num_query_frags - qe
		self.query_length = sum(qfrags)
		self.query_length_aligned = sum(c.ql for c in self.chunks_merged)
		self.ref_length_aligned = sum(c.rl for c in self.chunks_merged)
		self.frac_length_aligned = self.query_length_aligned/self.query_length
		self.frac_frags_aligned = float(qe - qs)/self.num_query_frags
		self.query_scaling_factor = self.query_length_aligned/self.ref_length_aligned

		self.query_misses = sum(c.query_misses for c in self.chunks_merged)
		self.ref_misses = sum(c.ref_misses for c in self.chunks_merged)
		self.matched_sites = len(self.chunks_merged)+1
		self.query_miss_rate = float(self.query_misses)/(self.query_misses + self.matched_sites)
		self.ref_miss_rate = float(self.ref_misses)/(self.ref_misses + self.matched_sites)

		# Orient chunks for maligner
		self.chunks_maligner = self.chunks_merged

		if self.orientation == "R":
			self.chunks_maligner = deepcopy(self.chunks_merged)
			nqf = self.num_query_frags
			for cm, cs in zip(self.chunks_maligner, self.chunks_merged):
				cm.q_s = nqf - cs.q_e
				cm.q_e = nqf - cs.q_s

	@property
	def chunk_string(self):
		# Write maligner chunk string.
		def c_to_str(c):
			fields = [c.q_s, c.q_e, c.ql, c.r_s, c.r_e, c.rl]
			return ','.join([str(f) for f in fields])
		return ';'.join(c_to_str(c) for c in self.chunks_maligner)

	_node_parsers = {
		'uuid' : SimpleField(),
		'reference_map': _parse_reference_map_field,
		'aligned_map':  _parse_aligned_map,
		'soma_score' : SimpleField(float),
		'count' : SimpleField(int),
		'f' : ListNodeField(_parse_f)
	}



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


class AlignmentWriter(object):
	"""Write alignment records in the maligner format"""
	
	# Maligner output to Alignment attr
	maligner_fields_mapping = [
		('query_map', )
		('ref_map', ),
		('is_forward', ),
		('num_matched_chunks', ),
		('query_misses', ),
		('ref_misses', ),
		('query_miss_rate', ),
		('ref_miss_rate', ),
		('total_score', ),
		('total_rescaled_score', ),
		('sizing_score', ),
		('sizing_score_rescaled', ),
		('query_scaling_factor', ),
		('chunk_string', )
	]


	]

	def __init__(self):
		pass

	def to_maligner(self, aln):



