import os

# Only import the malignpy_wrapper, which wraps C++ code using the BoostPython,
# if this environment variable is defined
if 'USE_MALIGNPY_PYTHON' in os.environ:
  from .malignpy_wrapper import *
from .core import *
import maps
