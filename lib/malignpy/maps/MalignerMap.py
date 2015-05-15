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
import numpy as np
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

    

def smooth(inputMap, minFrag=500, has_boundary = False):
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

#########################################################
# Smooth an input map in a symmetric manner, such that the smoothed version of
# the reverse map is the same. We do this by comparing the first frag to the last frag in size,
# and flipping if the first frag is less than the last frag.
def smooth_symmetric(inputMap, minFrag=500):

    if len(inputMap.frags) < 2:
        return inputMap

    first_frag = inputMap.frags[0]
    last_frag =  inputMap.frags[-1]
    if first_frag < last_frag:
        output_map = MalignerMap(frags = inputMap.frags[::-1],
                mapId = inputMap.mapId)
        output_map = smooth_without_boundary(output_map, minFrag = minFrag)
        # Reverse the frags again
        output_map.frags = output_map.frags[::-1]
        return output_map
 
    output_map = smooth_without_boundary(inputMap, minFrag = minFrag)
    return output_map

#########################################################
# Smooth an input map in an iterative fashion, starting with
# with the smallest fragment. We smooth by merging with the smaller
# neighboring fragment.
#
# This strategy will do a better job of producing consistent smoothed maps
# for a reference and forward and reverse subqueries of the reference.
def smooth_left_right(inputMap, minFrag=500, boundary_frags = True):

    from numpy import array, concatenate as concat, sum

    if len(inputMap.frags) < 2:
        return inputMap

    frags = np.array(inputMap.frags)

    # Remove boundary fragments - ignore them while smoothing
    if boundary_frags:

        bleft = frags[0]
        bright = frags[-1]
        frags = frags[1:-1]


    while True:

        num_frags = frags.shape[0]
        if num_frags < 2:
            break

        i = np.argmin(frags)
        f = frags[i]
        last_ind = num_frags - 1

        if f > minFrag:
            break
        if (i == 0):
            # This is the first frag, merge right
            frags = concat(([frags[0] + frags[1]], frags[2:]))
        elif (i == last_ind):
            # This is the laste frag, merge left
            frags = concat((frags[:-2], [frags[-1] + frags[-2]]))
        else:
            # This is an interior frag, merge with min
            left_frag = frags[i-1]
            right_frag = frags[i+1]
            if left_frag < right_frag:
                frags = concat((frags[0:i-1], [frags[i-1] + frags[i]], frags[i+1:]))
            else:
                frags = concat((frags[0:i], [frags[i] + frags[i+1]], frags[i+2:]))
    
    # Reattach boundary fragments
    if boundary_frags:
        frags = concat(([bleft], frags, [bright]))

    output_map = MalignerMap(frags = frags, mapId = inputMap.mapId)
    return output_map


