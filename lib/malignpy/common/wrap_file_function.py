"""
Utilities for working with alignment json files.
"""
from functools import wraps

class wrap_file_function(object):
  """
  Wrap a function which takes a file or a str as it's first argument.
  If a str is provided, replace the first argument of the wrapped function
  with a file handle, and close the file afterwards

  Example:

  @wrap_file_function('w')
  def write_hi(f):
    f.write('hi!\n')

  # This will write to already open file handle.
  f = open('f1.txt', 'w')
  write_hi(f)
  f.close()

  # This will open file f2.txt with mode 'w', write to it, and close the file.
  write_hi('f2.txt')
  """

  def __init__(self, *args):
    self.modes = args if args else ('r',)

  def __call__(self, func):

    @wraps(func)
    def wrapped(*args, **kwargs):

      close = [] # Files that should be closed
      files = [] # File handles that should be passed to func
      num_files = len(self.modes)

      for i, mode in enumerate(self.modes):

        fp = args[i]

        if isinstance(fp, str):
          fp = open(fp, mode)
          close.append(fp)

        files.append(fp)

      try:

        # Replace the files in args when calling func
        args = files + list(args[num_files:])

        # Make function call and return value
        return func(*args, **kwargs)

      finally:

        for fp in close: 
            fp.close()

    return wrapped



if __name__ == "__main__":
 
  # Demonstration of wrapping a function which writes to a file.
  print '-'*50
  @wrap_file_function('w')
  def write_hi(f):
    f.write('hi!\n')
 
  f = open('temp.txt', 'w')
  write_hi(f)
  write_hi(f)
  write_hi(f)
  f.close()
 
  write_hi('temp2.txt')
  
 
  # Demonstration of wrapping a function which reads from a file.
  print '-'*50
  @wrap_file_function()
  def read_file(f):
      print f.read()
 
  f = open('temp.txt')
  print 'Reading file temp.txt from handle f:'
  read_file(f)
  print 'Reading file temp2.txt'
  read_file('temp2.txt')

  # Demonstration of wrapping a function takes multiple files
  print '-'*50
  @wrap_file_function('r', 'r')
  def read_files(f1, f2):
      print 'reading f1: '
      print f1.read()
      print 'reading f2: '
      print f2.read()

  @wrap_file_function('r', 'w')
  def read_write(f1, f2):
    f2.write(f1.read())
 
  read_files(open('temp.txt'), open('temp2.txt'))
  read_files('temp.txt', 'temp2.txt')
  read_write('temp.txt', 'temp.copy.txt')

  import sys
  print 'writing temp.txt to stdout:'
  read_write('temp.txt', sys.stdout)

  print 'Contents of temp.copy.txt:'
  read_file('temp.copy.txt')

  @wrap_file_function('w')
  def throw_exception(f):
    raise RuntimeError('BLAH!')

  #print '-'*50
  #throw_exception('exception.txt')





