############################
# File: SOMAMap.py
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

class SOMAMap(object):

    def __init__(self, *args, **kwargs):
        if 'line' in kwargs:
            self.makeFromLine(kwargs['line'])
        else:
            self.makeFromAttributes(**kwargs)
        self.checkMap()

    # Create a SOMAMap from a line in maps file
    def makeFromLine(self, line):
        """
        Create a SOMAMap from a line in a SOMA Maps file.
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

    # write the SOMAMap to file
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
        Return dict representation of a SOMAMap
        """
        fields = ['mapId', 'length', 'numFrags', 'frags']
        d = { k:getattr(self,k) for k in fields }
        return d

    # Check the consistency of the object
    def checkMap(self):
        if len(self.frags) != self.numFrags:
            raise Exception('SOMAMap attributes are inconsistent!')


#########################################################
# Smooth a map by removing small restriction fragments.
def smooth(inputMap, minFrag=500):
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
   
    sys.stdout.write('smooth merged %i fragments out of %i.\n'%(mergeCount, len(inputMap.frags)))
    outputMap = SOMAMap(frags=frags, mapId = inputMap.mapId)
    assert(sum(outputMap.frags) == sum(inputMap.frags))
    assert(outputMap.length == inputMap.length)
    return outputMap
