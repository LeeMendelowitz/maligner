#!/usr/bin/env python
"""
Convert maps from the Schwartz lab format to the maligner format.
"""
import argparse, sys, os
from malignpy.maps.SOMAMap import SOMAMap
from malignpy.maps.utils import read_map_data_schwartz_format
from malignpy.common import wrap_file_function

@wrap_file_function('r', 'w')
def run(input_maps_file, fout):
  map_gen = read_map_data_schwartz_format(input_maps_file)
  for map in map_gen:
    map.write(fout)


if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="""Convert maps files from the Shwartz lap maps format to the maligner maps format.""")
  parser.add_argument('maps_file', metavar='MAPS_FILE', type=str,
                   help='Input maps file in the Schwartz lab format.')
  parser.add_argument('-o', '--output', metavar='OUTPUT FILE',
    help = "Output file. (Default: STDOUT)", default = sys.stdout)
  args = parser.parse_args()
  run(args.maps_file, args.output)
