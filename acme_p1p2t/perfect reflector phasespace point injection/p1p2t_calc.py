import numpy, h5py, math

def boris_push(x, p, t, dt, fields):
  u_minus = p + fields.E(x, t) * dt/2
  tau = fields.B(x, t) * dt/2 * (1.0 + p.dot(p))**(-0.5)
  sigma = 2 * tau / (1 + tau.dot(tau))
  u_plus = numpy.cross(u_minus, sigma)
  u_plus += u_minus + tau * u_minus.dot(sigma)
  p = u_plus + fields.E(x, t) * dt/2
  x += dt * p * (1.0 + p.dot(p))**(-0.5)
  return x, p

class fields(object):
  def __init__(self, a0):
    self.a0 = a0
  def E(self, x, t):
    E = numpy.zeros(3)
    if( x[0] < 0.):
      E[1] = self.a0 * (math.sin(x[0]-t) - math.sin(-x[0]-t))
    return E
  def B(self, x, t):
    B = numpy.zeros(3)
    if( x[0] < 0.):
      B[2] = self.a0 * (math.sin(x[0]-t) + math.sin(-x[0]-t))
    return B

def find_reentry_energy(p, t):
  t0 = t
  #a0, dt is hardwired here, clearly not the best way to do this
  dt = math.pi/100.
  x = numpy.zeros(3)
  laser_field = fields(6.0)
  while ( t < t0 + 2*math.pi):
    x, p = boris_push(x, p, t, dt, laser_field)
    t += dt
    if( x[0] > 0.):
      energy = p.dot(p)/(1. + (p.dot(p) + 1)**.5)
      return energy, False
  return -1.0, True

def find_reentry_time(p, t):
  t0 = t
  #a0, dt is hardwired here, clearly not the best way to do this
  dt = math.pi/100.
  x = numpy.zeros(3)
  laser_field = fields(6.0)
  while ( t < t0 + 4*math.pi):
    x, p = boris_push(x, p, t, dt, laser_field)
    t += dt
    if( x[0] > 0.):
      return t-t0, False
  return t-t0, True
  
def find_reentry_time_energy(p, t):
  t0 = t
  #a0, dt is hardwired here, clearly not the best way to do this
  dt = math.pi/100.
  x = numpy.zeros(3)
  laser_field = fields(6.0)
  while ( t < t0 + 2*math.pi):
    x, p = boris_push(x, p, t, dt, laser_field)
    t += dt
    if( x[0] > 0.):
      energy = p.dot(p)/(1. + (p.dot(p) + 1)**.5)
      return t-t0, energy, False
  return t-t0, -1.0, True

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

import multiprocessing
def find_te_tuple_wrapper(tuppy):
  return find_reentry_time_energy(tuppy[0], tuppy[1])

def par_tem_slice(t, p1_max, p2_max, dp):
  p1_n = int(p1_max/dp)
  p2_n = 2 * int(p2_max/dp)+1
  energy_slice = numpy.zeros((p1_n, p2_n))
  time_slice = numpy.zeros((p1_n, p2_n))
  mask = numpy.zeros((p1_n, p2_n))
  #multiprocessing.Pool.map does not seem to like lambda functions
  #find_energy_wrapper = lambda p2: find_reentry_energy(p2, t)
  #and it doesn't like locally defined functions? mrr
  #def find_energy_wrapper(p):
  #  return find_reentry_energy(p, t)
  pool = multiprocessing.Pool(4)
  for i in range(p1_n):
    p = [numpy.array((-(i+1)*dp, -p2_max + dp*j, 0.,)) for j in range(p2_n)]
    the_result = pool.map(find_te_tuple_wrapper, zip(p, [t for dummy in range(p2_n)]))
    for j in range(p2_n):
      ene_init = p[j].dot(p[j])/( 1 + ( 1 + p[j].dot(p[j]) )**.5)
      time_slice[i,j] = the_result[j][0]
      energy_slice[i,j] = the_result[j][1] - ene_init
      if (the_result[j][2]):
        mask[i,j] = 1.0
  pool.close()
  return time_slice, energy_slice, mask

if __name__ == "__main__":
  import matplotlib.pyplot as plt
  #from pylab import cm, show
  p1_max, p2_max, dp = 12., 6., 0.05
  import os
  current_directory = os.getcwd()
  num_slices = 100
  hdf5_file = h5py.File(current_directory + '/p1p2t.h5', 'w')
  hdf5_file.attrs.create('p1_max', p1_max)
  hdf5_file.attrs.create('p2_max', p2_max)
  reentry_energy_group = hdf5_file.create_group('reetry_energy')
  reentry_time_group = hdf5_file.create_group('reetry_time')
  escape_mask_group = hdf5_file.create_group('escape_mask')
  for i in range(num_slices):
    t = 2*math.pi*i/num_slices
    time_slice, ene_slice, mask = par_tem_slice(t, p1_max, p2_max, dp)
    #
    this_iter_time_dataset = reentry_time_group.create_dataset( "%08d" % i, data = time_slice)
    this_iter_time_dataset.attrs.create('iteration',i)
    this_iter_time_dataset.attrs.create('time', t)
    #
    this_iter_ene_dataset = reentry_energy_group.create_dataset( "%08d" % i, data = ene_slice)
    this_iter_ene_dataset.attrs.create('iteration',i)
    this_iter_ene_dataset.attrs.create('time', t)
    #
    this_iter_mask_dataset = escape_mask_group.create_dataset( "%08d" % i, data = mask)
    this_iter_mask_dataset.attrs.create('iteration',i)
    this_iter_mask_dataset.attrs.create('time', t)
    #
    print "Finished writing the " + str(i) + "th slice"

  hdf5_file.close()