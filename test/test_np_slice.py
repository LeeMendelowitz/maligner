import numpy as np
import timeit

nr = 20
nc = 40000
b = np.arange(nr*nc).reshape(nr,nc)

def get_special_diag(a):
  m,n = a.shape
  nr, nc = m, n - m + 1
  a_flat = a.flatten()
  inds = np.tile(np.arange(nc), nr) + np.arange(nr).repeat(nc) + n*np.arange(nr).repeat(nc)
  ret = a_flat[inds].reshape(nr, nc)
  return ret

def get_special_diag2(a):
  m,n = a.shape
  nr, nc = m, n - m + 1
  r = np.arange(nr).repeat(nc)
  c = np.tile(np.arange(nc), nr) + np.arange(nr).repeat(nc)
  ret = a[r,c].reshape(nr, nc)
  return ret

def get_special_diag3(a):
  m,n = a.shape
  nr, nc = m, n - m + 1
  ret = np.zeros((nr, nc), dtype = a.dtype)
  for i in xrange(nc):
    ret[:,i] = a.diagonal(i)
  return ret


setup="""
import numpy as np
nr = 20
nc = 40000
b = np.arange(nr*nc).reshape(nr,nc)

def get_special_diag(a):
  m,n = a.shape
  nr, nc = m, n - m + 1
  a_flat = a.flatten()
  inds = np.tile(np.arange(nc), nr) + np.arange(nr).repeat(nc) + n*np.arange(nr).repeat(nc)
  ret = a_flat[inds].reshape(nr, nc)
  return ret
def get_special_diag2(a):

  m,n = a.shape
  nr, nc = m, n - m + 1
  r = np.arange(nr).repeat(nc)
  c = np.tile(np.arange(nc), nr) + np.arange(nr).repeat(nc)
  ret = a[r,c].reshape(nr, nc)
  return ret

def get_special_diag3(a):
  m, n = a.shape
  nr, nc = m, n - m + 1
  ret = np.zeros((nr, nc))
  for i in xrange(nc):
    ret[:,i] = a.diagonal(i)
  return ret

"""

stmt1 = "ret = get_special_diag(b)"
stmt2 = "ret = get_special_diag2(b)"
stmt3 = "ret = get_special_diag3(b)"

def run():
  t1 = timeit.timeit(stmt1, setup, number = 100)
  t2 = timeit.timeit(stmt2, setup, number = 100)
  t3 = timeit.timeit(stmt3, setup, number = 100)
  return (t1, t2, t3)
  