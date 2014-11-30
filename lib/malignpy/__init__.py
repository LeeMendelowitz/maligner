import os
if not 'NO_MALIGNPY' in os.environ:
  from .malignpy_wrapper import *
from .core import *
import maps
