import sys
import re

from .SOMAMap import SOMAMap
from ..common import wrap_file_function

__all__ = (
    "read_maps",
    "gen_maps",
    "write_maps",
    "convert_maps_file_to_mongo"
)

###################################################
# Convert an optical map from the Schwartz lab format
# to the SOMA format
# opticalMapFile: optical map file in the Schwartz lab format
def convert_optical_maps(opticalMapFile, outputPfx):
    opMapFileOut = '%s.opt'%outputPfx

    msg = '\n'+'*'*50 + \
          '\nReading Optical Map File %s\n'%opticalMapFile + \
          '*'*50 + '\n'
    sys.stderr.write(msg)

    opMapList = readMapDataSchwartz(opticalMapFile)
    enzymeSet = set(om.enzyme for om in opMapList)
    if len(enzymeSet) > 1:
        raise RuntimeError('Different enzymes used in the input optical map set!')
    enzyme = opMapList[0].enzyme

    msg = '\n'+'*'*50 +\
          '\nConverting Optical Map to SOMA Format\n' +\
          '*'*50 + '\n'
    sys.stderr.write(msg)

    # Optical maps for chromosomes 
    # Remove all white space from restriction map names
    for opMap in opMapList:
        opMap.mapId = ''.join(opMap.mapId.split())
    writeMaps(opMapList, opMapFileOut)
    result = { 'enzyme' : enzyme,
               'opMapList' : opMapList,
               'opticalMapFile' : opMapFileOut}
    return result

#############################################################################################
# Read an optical map in the Schwartz format
# Return a list of OpticalMapData instances
def readMapDataSchwartz(filename, verbose = False):
    omaps = []
    fin = open(filename)
    mapGen = mapDataSchwartzGen(fin)

    if verbose:
        omaps = []
        for map in mapGen:
            omaps.append(map)
            print 'Read map from chromosome %s with enzyme %s and %i fragments'%(map.mapId, map.enzyme, len(map.frags))
    else:
        omaps = [m for m in mapGen]
    fin.close()
    return omaps

def mapDataSchwartzGen(handle):
    fin = handle
    for line in fin:
        # Skip if line is empty
        if not line.strip():
            continue
        line2 = fin.next()
        # Read two consecutive lines
        # Line 1: Chromosome name (this is line)
        # Line 2: Enzyme, 'B', List of fragment lengths (this is line2)
        mapId = line.strip()
        fields = line2.strip().split()
        enzyme = fields[0]
        fragLengths = [int(1000.0*float(field)) for field in fields[2:]] # Convert to bp
        omap = SOMAMap(mapId = mapId, frags = fragLengths, enzyme = enzyme)
        yield omap
    

#############################################################################################
# Reads (in a crude fashion) an optical map in xml format
# Return a list of OpticalMapData instances
def read_map_data_xml(fileName):
    fin = open(fileName)
    opticalMapList = []

    # Define regular expression for a line with map data
    # (frag_number) (size) (std. dev) (pfx) (sfx)
    # Note: pfx and sfx are optional!
    fragmentPattern = re.compile(r'S="(\d+)" STDDEV="([\d\.]+)"') # Finds fragment length and standard deviation
    startMapPattern = re.compile(r'<RESTRICTION_MAP ID="([^"]+)" ENZYME="([^"]+)"') # Start of a restriction map
    endMapPattern = re.compile(r'</RESTRICTION_MAP>') # End of a restriction map

    fragLengths = []
    fragSDs = []
    enzyme = ''
    mapName = ''
    for line in fin:
        startMapMatch = startMapPattern.search(line)
        endMapMatch = endMapPattern.search(line)
        fragmentMatch = fragmentPattern.search(line)
        # Assert that only at most one of these patterns has a match
        assert sum(map(int, [m != None for m in [startMapMatch, endMapMatch, fragmentMatch]])) <= 1

        if startMapMatch:
            # Assert that the previous map was ended properly
            assert(len(fragLength)==0)
            assert(len(fragSD)==0)
            assert(len(enzyme)==0)
            assert(len(mapName)==0)
            mapName = startMapMatch.group(1)
            enzyme = startMapMatch.group(2)
        elif fragmentMatch:
            fragLength = float(fragmentMatch.group(1))/1000 # in kbp
            fragSD = float(fragmentMatch.group(2)) # in kbp
            if fragSD==0.000:
                print 'Warning: Ignoring fragment with SD=0.000 in optical map %s'%mapName
            else:
                fragLengths.append(fragLength)
                fragSDs.append(fragSD)
        elif endMapMatch:
            print 'Creating optical map with name %s enzyme %s and %i fragments'%(mapName, enzyme, len(fragLengths))
            opticalMapList.append(SOMAMap(mapId = mapName, frags = fragLengths, fragSD = fragSD, enzyme = enzyme))
            fragLengths = []
            fragSDs = []
            enzyme = ''
            mapName = ''
    fin.close()

    # Assert that the previous map was ended properly
    assert(len(fragLengths)==0)
    assert(len(fragSDs)==0)
    assert(len(enzyme)==0)
    assert(len(mapName)==0)
    return opticalMapList

#########################################################
# Read map file in the SOMA map format
@wrap_file_function('r')
def read_maps(fin):
    """
    Read all SOMA maps in the file and return as a dict.

    Usage:
      maps = read_maps('file.maps')

      # or
      f = open('file.maps')
      maps = readMaps(f)
      f.close()
    """
    maps = [SOMAMap(line=l) for l in fin]
    mapDict = dict((m.mapId, m) for m in maps)
    return mapDict

#############################################
def gen_maps(f):
    """
    Generate maps from the SOMAMap file specified by f.
    f can be a filename str or a file handle.
    """

    should_close = False
    if isinstance(f, str):
        f = open(f)
        should_close = True

    lines = (l for l in f if l)
    for l in lines:
        yield SOMAMap(line=l)

    if should_close:
        f.close()

#########################################################
# Write maps to a file in the SOMA map format
def write_maps(mapList, handle):
    fout = handle
    closeFile = False
    if type(handle) is str:
        fout = open(handle, 'w')
        closeFile = True

    for map in mapList:
        map.write(fout)

    if closeFile:
        fout.close()

#########################################################
@wrap_file_function('r', 'w')
def convert_maps_file_to_mongo(maps_in, json_out):
    map_gen = gen_maps()
