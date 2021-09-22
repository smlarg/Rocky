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

def slice(t, p1_max, p2_max, dp):
  p1_n = int(p1_max/dp)
  p2_n = 2 * int(p2_max/dp)+1
  slice = numpy.zeros((p1_n, p2_n))
  p_init = numpy.zeros(3)
  for i in range(p1_n):
    p_init[0] = -(i+1)*dp
    for j in range(p2_n):
      p_init[1] = -p2_max + j*dp
      ene_init = p_init.dot(p_init) / (1. + ( 1+p_init.dot(p_init) )**.5 )
      ene_final, escapes = find_reentry_energy(p_init, t)
      if (not escapes):
        slice[i,j] = ene_final - ene_init
      else:
        continue
  return slice

import multiprocessing
def find_energy_tuple_wrapper(tuppy):
  return find_reentry_energy(tuppy[0], tuppy[1])

def par_slice(t, p1_max, p2_max, dp):
  p1_n = int(p1_max/dp)
  p2_n = 2 * int(p2_max/dp)+1
  slice = numpy.zeros((p1_n, p2_n))
  mask = numpy.zeros((p1_n, p2_n))
  #multiprocessing.Pool.map does not seem to like lambda functions
  #find_energy_wrapper = lambda p2: find_reentry_energy(p2, t)
  #and it doesn't like locally defined functions? mrr
  #def find_energy_wrapper(p):
  #  return find_reentry_energy(p, t)
  pool = multiprocessing.Pool(4)
  for i in range(p1_n):
    p = [numpy.array((-(i+1)*dp, -p2_max + dp*j, 0.,)) for j in range(p2_n)]
    the_result = pool.map(find_energy_tuple_wrapper, zip(p, [t for dummy in range(p2_n)]))
    for j in range(p2_n):
      ene_init = p[j].dot(p[j])/( 1 + ( 1 + p[j].dot(p[j]) )**.5)
      if (not the_result[j][1]):
        slice[i,j] = the_result[j][0] - ene_init
      else:
        mask[i,j] = 1.0
  pool.close()
  return slice, mask

if __name__ == "__main__":
  import matplotlib.pyplot as plt
  #from pylab import cm, show
  p1_max, p2_max, dp = 12., 6., 0.1
  import os
  pathToFrames = os.getcwd()
  no_o_frames = 20
  for i in range(no_o_frames):
    t = math.pi*i/no_o_frames
    frame, mask0 = par_slice(t, p1_max, p2_max, dp)
    #fig, ax = plt.subplots()
    #im = ax.imshow(slice0, cmap=cm.RdBu, vmin=0., vmax=12., extent=[-p2_max, p2_max, -p1_max,0.])
    #cb = fig.colorbar(im, ax=ax)
    #ax.set_title("t = " + str(t))
    #show()
    #
    X = numpy.zeros(mask0.shape)
    Y = numpy.array((X,X,X+1,mask0))
    Y = Y.swapaxes(0,1)
    Y = Y.swapaxes(1,2)
    #
    plt.figure(figsize=(6,4), dpi = 100)
    plt.axes([0.15, 0.15, 0.75, 0.75])
    plt.xlim([-p2_max,p2_max])
    plt.ylim([-p1_max,0.0])
    plt.title('time = ' + str(t) )
    plt.imshow(frame, cmap=plt.cm.hot, vmin=0., vmax=frame.max(), extent=[-p2_max, p2_max, -p1_max,0.])
    plt.colorbar()
    plt.imshow(Y, extent=[-p2_max, p2_max, -p1_max,0.])
    index = "%03d" % i
    plt.savefig( pathToFrames + "/" + "p1p2t_frame" + index + ".jpg")
    plt.close()
    print "Finished writing the " + str(i) + "th frame"
