############################
# File: MalignerMap.py
# Author: Lee Mendelowitz (lmendelo@umiacs.umd.edu)
#
# Description:
# Define the SOMA Map class for reading, writing, and creating restriction map
# objects.
#
#
# The SOMA map file is tab delimited with one record per line:
# [MapId] [Length (bp)] [number of fragments] [frag 0 length (bp)] [frag 1 length (bp)] ...
import sys

from ..common import wrap_file_function

class MalignerMap(object):

    def __init__(self, *args, **kwargs):
        if 'line' in kwargs:
            self.makeFromLine(kwargs['line'])
        else:
            self.makeFromAttributes(**kwargs)
        self.checkMap()

    # Create a MalignerMap from a line in maps file
    def makeFromLine(self, line):
        """
        Create a MalignerMap from a line in a SOMA Maps file.
        """
        fields = line.strip().split()
        self.mapId = fields[0]
        self.length = int(fields[1])
        self.numFrags = int(fields[2])
        self.frags = [int(f) for f in fields[3:]]

    # Create from mapId and frags attribute
    def makeFromAttributes(self, **kwargs):
        self.frags = list(kwargs['frags'])
        self.mapId = kwargs['mapId']
        self.numFrags = len(self.frags)
        self.length = sum(self.frags)
       
        # Add any other attributes from kwargs 
        for attr in kwargs.iterkeys():
            if attr not in self.__dict__:
                self.__dict__[attr] = kwargs[attr]

    # write the MalignerMap to file
    def write(self, handle):
        f = handle
        fields = [self.mapId,
                  str(self.length),
                  str(len(self.frags))]
        fields.extend(str(frag) for frag in self.frags)
        outS = '\t'.join(fields) + '\n'
        f.write(outS)

    def get_data(self):
        """
        Return dict representation of a MalignerMap
        """
        fields = ['mapId', 'length', 'numFrags', 'frags']
        d = { k:getattr(self,k) for k in fields }
        return d

    # Check the consistency of the object
    def checkMap(self):
        if len(self.frags) != self.numFrags:
            raise Exception('MalignerMap attributes are inconsistent!')

    @property
    def name(self):
        return self.mapId

    

def smooth(inputMap, minFrag=500, has_boundary = True):
    if has_boundary:
        return smooth_with_boundary(inputMap, minFrag)
    else:
        return smooth_without_boundary(inputMap, minFrag)

#########################################################
# Smooth a map by removing small restriction fragments.
def smooth_without_boundary(inputMap, minFrag=500):
    frags = []
    mergeCount = 0
    fragGen = iter(inputMap.frags)
    curFrag = []
    for frag in fragGen:
        curFrag = [frag]
        curLength = frag
        while curLength < minFrag:
            try:
                nextFrag = fragGen.next()
                curFrag.append(nextFrag)
                curLength += nextFrag
                mergeCount += 1
            except StopIteration:
                break
        assert(sum(curFrag) == curLength)
        frags.append(curLength)
        curFrag = []
   
    sys.stderr.write('smooth merged %i fragments out of %i.\n'%(mergeCount, len(inputMap.frags)))
    outputMap = MalignerMap(frags=frags, mapId = inputMap.mapId)
    assert(sum(outputMap.frags) == sum(inputMap.frags))
    assert(outputMap.length == inputMap.length)
    return outputMap

#########################################################
# Smooth a map by removing small restriction fragments.
# Leave boundary fragments as is
def smooth_with_boundary(inputMap, minFrag=500):
    frags = []
    mergeCount = 0

    if not inputMap.frags:
        return MalignerMap(frags = [], mapId = inputMap.mapId)

    first_frag = inputMap.frags[0]
    last_frag = inputMap.frags[-1] if len(inputMap.frags) > 1 else None

    fragGen = iter(inputMap.frags[1:-1])

    frags.append(first_frag)

    # Smooth interior fragments
    curFrag = []
    for frag in fragGen:
        curFrag = [frag]
        curLength = frag
        while curLength < minFrag:
            try:
                nextFrag = fragGen.next()
                curFrag.append(nextFrag)
                curLength += nextFrag
                mergeCount += 1
            except StopIteration:
                break
        assert(sum(curFrag) == curLength)
        frags.append(curLength)
        curFrag = []

    # Append last fragment
    if last_frag:
        frags.append(last_frag)
   
    sys.stderr.write('smooth merged %i fragments out of %i.\n'%(mergeCount, len(inputMap.frags)))
    outputMap = MalignerMap(frags=frags, mapId = inputMap.mapId)
    assert(sum(outputMap.frags) == sum(inputMap.frags))
    assert(outputMap.length == inputMap.length)
    return outputMap
