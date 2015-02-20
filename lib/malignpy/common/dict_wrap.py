

class DictWrap(object):
  """
  Wrap a dictionary as a class.
  """
  def __init__(self, d):
    for k,v in d.iteritems():
      setattr(self, k, v)