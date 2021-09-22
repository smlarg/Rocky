import numpy, math
from random import random
r = lambda : random() - 0.5

def boris_push(x, p, t, dt, fields):
  u_minus = p + fields.E(x, t) * dt/2
  tau = fields.B(x, t) * dt/2 * (1.0 + p.dot(p))**(-0.5)
  sigma = 2 * tau / (1 + tau.dot(tau))
  u_plus = numpy.cross(u_minus, sigma)
  u_plus += u_minus + tau * u_minus.dot(sigma)
  p = u_plus + fields.E(x, t) * dt/2
  x += dt * p * (1.0 + p.dot(p))**(-0.5)
  return x, p

def boris_push2(x, p, dt, E, B):
  """
  x, p, E, B arrays of shape(n,3)
  a.dot(b.T)[:,0:1] returns an array of shape(n,1), so for instance
  v = p * ( 1 + p.dot(p.T)[:,0:1])**(-0.5) is a shape(n,3) array as well
  it seems like there must be a simpler way to do that, it seems like a logical thing to
  do, but I don't see it right now
  """
  u_minus = p + E * dt/2
  tau = B * dt/2 * (1.0 + p.dot(p.T)[:,0:1])**(-0.5)
  sigma = 2 * tau / (1 + tau.dot(tau))
  u_plus = numpy.cross(u_minus, sigma)
  u_plus += u_minus + tau * u_minus.dot(sigma.T)[:,0:1]
  p = u_plus + E * dt/2
  x += dt * p * (1.0 + p.dot(p.T)[:,0:1])**(-0.5)
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

def dist_func(p):
  if p.dot(p) > 3.5**2:
    return 0.
  return math.exp(-p.dot(p)/0.548)

#if __name__ == "__main__":
if True:
  p_range = 3.5
  xmin = 0.0
  xmax = xmin + 20.0
  dp = 0.05
  dx = 0.025
  
  block_size = 2048
  
  np = int( 2 * p_range/ dp)
  nx = int( (xmax - xmin) / dx)
  
  #p_array = numpy.zeros((np * np * nx, 3), order = 'FORTRAN')
  p_array = numpy.zeros((np * np * nx, 3))
  #x_array = numpy.zeros((np * np * nx, 3), order = 'FORTRAN')
  x_array = numpy.zeros((np * np * nx, 3))
  q_array = numpy.zeros((np * np * nx))
  
  index = 0
  for i in range(nx):
    #x = numpy.array(((i+r())*dx + xmin,0.,0.))
    for j in range(np):
      for k in range(np):
        p = numpy.array( (j+r(),k+r(),0.5*np) )
        x = numpy.array(((i+r()+0.5)*dx + xmin,0.,0.))
        p *= dp
        p -= p_range
        q = dist_func(p)
        if q > 0.:
          p_array[index] = p
          x_array[index] = x
          q_array[index] = q
          index += 1
  
  #real_particle_count = index
  #particle_count = real_particle_count + -real_particle_count%block_size
  particle_count = index
  p_array = p_array[:particle_count]
  x_array = x_array[:particle_count]
  q_array = q_array[:particle_count]
  
  t = 0.
  t_max = 10.
  dt = 0.01
  iter = 0
  
  import time
  time0 = time.clock()
  
  import rocky_py
  #E = numpy.zeros((block_size,3))
  #B = numpy.zeros((block_size,3))
  
  while t < t_max:
    print "at time = " + str(t)
#    for i in xrange(particle_count/block_size):
#      lb, ub = i*block_size, (i+1)*block_size
#      E, B = rocky_py.get_block_of_fields(x_array[lb:ub],t)
#      x, p = x_array[lb:ub], p_array[lb:ub]
#      x, p = rocky_py.push_block(x, p, dt, E, B)
#      x_array[lb:ub], p_array[lb:ub] = x, p
    x_array, p_array = rocky_py.pushall(x_array, p_array, t, dt)
    t += dt
    iter += 1
    time1 = time0
    time0 = time.clock()
    print "timestep length = " + str(time0 - time1) + "s"
  
  import h5py, os
  current_directory = os.getcwd()
  hdf5_file = h5py.File(current_directory + '/particle_data.h5', 'w')
  xDS = hdf5_file.create_dataset('x', data = x_array)
  pDS = hdf5_file.create_dataset('p', data = p_array)
  qDS = hdf5_file.create_dataset('q', data = q_array)
  #
  hdf5_file.attrs.create('TIME',numpy.array((t,)))
  hdf5_file.attrs.create('ITER',numpy.array((iter,)))
  hdf5_file.attrs.create('DT',numpy.array((dt,)))
  #
  hdf5_file.close()
  	