from math import sqrt

def quadratic_formula(a,b,c):
  return (-b+(b**2 -4*a*c)**.5)/2/a, (-b-(b**2 -4*a*c)**.5)/2/a

def p1(p2, ene):
  p2 *= p2
  p1_2 = quadratic_formula(1, 2*p2 - 2*ene - ene**2, p2**2 - 2*p2*ene - p2*ene**2)[0]
  return p1_2**.5

def ene(p1, p2):
  e = (p1**2 + p2**2)/((p1**2 + p2**2 + 1)**.5 + 1)
  return e

# Testing function showing initial energy as a function of initial momentum
def slice_init(t, p1_max, p2_max, dp):
  p1_n = int(p1_max/dp)
  p2_n = 2 * int(p2_max/dp)+1
  slice = numpy.zeros((p1_n, p2_n))
  p_init = numpy.zeros(3)
  for i in range(p1_n):
    p_init[0] = -(i+1)*dp
    for j in range(p2_n):
      p_init[1] = -p2_max + j*dp
      ene_init = p_init.dot(p_init) / (1. + ( 1+p_init.dot(p_init) )**.5 )
      slice[i,j] = ene_init
  return slice

# super testing functions to show initial p1, p2
# (will indexing ever not confuse me? it seems unlikely)
def p1_init(t, p1_max, p2_max, dp):
  p1_n = int(p1_max/dp)
  p2_n = 2 * int(p2_max/dp)+1
  slice = numpy.zeros((p1_n, p2_n))
  p_init = numpy.zeros(3)
  for i in range(p1_n):
    p_init[0] = -(i+1)*dp
    for j in range(p2_n):
      p_init[1] = -p2_max + j*dp
      ene_init = p_init.dot(p_init) / (1. + ( 1+p_init.dot(p_init) )**.5 )
      slice[i,j] = p_init[0]
  return slice
def p2_init(t, p1_max, p2_max, dp):
  p1_n = int(p1_max/dp)
  p2_n = 2 * int(p2_max/dp)+1
  slice = numpy.zeros((p1_n, p2_n))
  p_init = numpy.zeros(3)
  for i in range(p1_n):
    p_init[0] = -(i+1)*dp
    for j in range(p2_n):
      p_init[1] = -p2_max + j*dp
      ene_init = p_init.dot(p_init) / (1. + ( 1+p_init.dot(p_init) )**.5 )
      slice[i,j] = p_init[1]
  return slice

if __name__ == "__main__":
  import h5py, os, numpy
  import matplotlib.pyplot as plt
  
  hdf5_file = h5py.File('p1p2t.h5','r')
  p1_max = hdf5_file.attrs.get('p1_max')
  p2_max = hdf5_file.attrs.get('p2_max')
  #
  mask_dataset = hdf5_file.get('escape_mask')
  mask_keys = mask_dataset.keys()
  mask = mask_dataset.get(mask_keys[0])[...]
  for key in mask_keys:
    mask = numpy.maximum(mask, mask_dataset.get(key)[...])
  X = numpy.zeros(mask.shape)
  Y = numpy.array((X,X,X+1,mask))
  Y = Y.swapaxes(0,1)
  Y = Y.swapaxes(1,2)
  #
  energy_dataset = hdf5_file.get('reetry_energy')
  ene_keys = energy_dataset.keys()
  ene_max_ds = energy_dataset.get(ene_keys[0])[...]
  for key in ene_keys:
    ene_max_ds = numpy.maximum(ene_max_ds, energy_dataset.get(key)[...])
  plt.figure(figsize=(6,4), dpi = 100)
  plt.axes([0.15, 0.15, 0.75, 0.75])
  plt.xlim([-p2_max,p2_max])
  plt.ylim([-p1_max,0.0])
  plt.imshow(ene_max_ds, cmap=plt.cm.hot, vmin=0., vmax=ene_max_ds.max(),
   extent=[-p2_max, p2_max, -p1_max,0.])
  plt.colorbar()
  plt.imshow(Y, extent=[-p2_max, p2_max, -p1_max,0.])
  plt.savefig('max_delta_e.jpg')
  plt.close()
  hdf5_file.close()
  #
  p1_max_ds = ene_max_ds.copy()
  dp = p1_max / p1_max_ds.shape[0]
  it = numpy.nditer(p1_max_ds, flags=['multi_index'], op_flags=['readwrite'])
  while not it.finished:
    p1_init = (it.multi_index[0]+1)*dp
    p2 = it.multi_index[1]*dp - p2_max
    delta_ene = it[0]
    if not mask[it.multi_index[0], it.multi_index[1]]:
      total_ene = delta_ene + ene(p1_init, p2)
      p1_final = p1(p2, total_ene)
      it[0] = p1_final
    else:
      it[0] = 0.
    it.iternext()
  plt.figure(figsize=(6,4), dpi = 100)
  plt.axes([0.15, 0.15, 0.75, 0.75])
  plt.xlim([-p2_max,p2_max])
  plt.ylim([-p1_max,0.0])
  plt.imshow(p1_max_ds, cmap=plt.cm.hot, vmin=0., vmax=p1_max_ds.max(),
   extent=[-p2_max, p2_max, -p1_max,0.])
  plt.colorbar()
  plt.imshow(Y, extent=[-p2_max, p2_max, -p1_max,0.])
  plt.savefig('max_final_p1.jpg')
  plt.close()
  #
  it = numpy.nditer(ene_max_ds, flags=['multi_index'], op_flags=['readwrite'])
  while not it.finished:
    p1_init = (it.multi_index[0]+1)*dp
    p2 = it.multi_index[1]*dp - p2_max
    if not mask[it.multi_index[0], it.multi_index[1]]:
      it[0] += ene(p1_init, p2)
    else:
      it[0] = 0.
    it.iternext()
  plt.figure(figsize=(6,4), dpi = 100)
  plt.axes([0.15, 0.15, 0.75, 0.75])
  plt.xlim([-p2_max,p2_max])
  plt.ylim([-p1_max,0.0])
  plt.imshow(ene_max_ds, cmap=plt.cm.hot, vmin=0., vmax=ene_max_ds.max(),
   extent=[-p2_max, p2_max, -p1_max,0.])
  plt.colorbar()
  plt.imshow(Y, extent=[-p2_max, p2_max, -p1_max,0.])
  plt.savefig('max_final_e.jpg')
  plt.close()
