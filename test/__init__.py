import sys, os
from os.path import split, abspath

file_dir = split(abspath(__file__))[0]
root_dir = split(file_dir)[0]
lib_dir = os.path.join(root_dir, 'lib')
sys.path.append(lib_dir)
