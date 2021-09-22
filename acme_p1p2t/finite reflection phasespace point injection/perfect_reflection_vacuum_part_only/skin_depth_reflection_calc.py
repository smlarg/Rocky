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

class reflected_laser_field(object):
  """
  A perfectly reflected laser (in terms of energy absorption)
  from a finite plasma conductor, with the reflection plane at x1 = 0
  and the laser incident from negative x1.
  Note that this is not really correct anyway - it ignores the vXB force,
  and we're mostly using this for a0 > 1, so we can't ignore that.
  But ... we did it anyway.
  """
  def __init__(self, a0, omega_0, omega_p):
    self.a0 = a0
    # usually one or the other of these is 1
    self.omega_p = omega_p
    self.omega_0 = omega_0
    self.k_0 = self.omega_0 # to, just possibly, avoid some confusion
    if omega_p <= omega_0:
      raise RuntimeError("An error")
    self.n_p = omega_p**2 # we normalize away the 4pi too? seems odd...
    self.n_0 = omega_0**2
    self.k_c = (self.omega_p**2 - self.omega_0**2)**.5 # k_conductor
                                         # not k_plasma because that could be w_p/c, which it's not (quite)
    self.phi = math.atan( - self.k_0 / self.k_c )
  def E(self, x, t):
    E = numpy.zeros(3)
    if( x[0] < 0.):
      E[1] = 2 * self.a0 * math.sin(self.k_0 * x[0] + self.phi) * math.sin(self.omega_0 * t)
    else:
      E[1] = 2 * self.a0 * (-self.omega_0/self.omega_p) * math.exp( -self.k_c * x[0]) * math.sin(self.omega_0 * t)
    return E
  def B(self, x, t):
    B = numpy.zeros(3)
    if( x[0] < 0.):
      B[2] = 2 * self.a0 * math.cos(self.k_0 * x[0] + self.phi) * math.cos(self.omega_0 * t)
    else:
      B[2] = -self.k_c * 2 * self.a0 * (-1.0/self.omega_p) * math.exp( -self.k_c * x[0]) * math.cos(self.omega_0 * t)
    return B
  def A(self, x, t):
    #return (1./self.omega_0) * self.E(x,t) * math.cos(self.omega_0 * t)/ math.sin(self.omega_0 * t)
    # nope, division by zero
    A = numpy.zeros(3)
    if( x[0] < 0.):
      A[1] = 2 * self.a0 *  (1.0/self.omega_0) * math.sin(self.k_0 * x[0] + self.phi) * math.cos(self.omega_0 * t)
    else:
      A[1] = 2 * self.a0 * (-1.0/self.omega_p) * math.exp( -self.k_c * x[0]) *          math.cos(self.omega_0 * t)
    return A    

def find_p_final(p_initial, t_initial, dt, fields):
  """
  returns the 'final' momentum, and whether that momentum was inside the plasma or not
  ('inside' being 10 skin depths within 1.5 laser cycles)
  'corrects' the p_initial to conserve canonical momentum
  (this only really makes sense in the low-a0 limit, but let's see what we get...)
  """
  t = t_initial
  x = numpy.zeros(3)
  p = numpy.array(p_initial)
  #p -= fields.A(x, t) # I hope that's minus...........
  enters_vacuum = False
  while t < t_initial + 3*math.pi :
    x, p = boris_push(x, p, t, dt, fields)
    t += dt
    if( x[0] > 10./fields.k_c):
      return p, True, enters_vacuum
    if x[0] < 0.:
      enters_vacuum = True
  return p, False, enters_vacuum

def ene_angle_slice(t, ene_max, ene_min, n_ene, n_theta, fields, result_array):
  """
  moves through a (log spaced) array of energies, and finds the minimum energy
  (from among those energies) which can give a final energy and angle,
  for a given time of injection but cycling over all angles.
  Now that I write what this does down it doesn't seem at all like the right way
  to decompose the dimensionality of the problem, but I did it.
  also this wont work with the rest of what's written now, because the minimum comparisons are convoluted
  """
  slice = result_array.data.copy()
  slice = 0.*slice + 2.*ene_max
  log_ene_min = math.log(ene_min)
  log_ene_max = math.log(ene_max)
  log_ene_range = log_ene_max - log_ene_min
  for i in range(n_ene):
    log_ene = log_ene_min + log_ene_range * i/(n_ene-1)
    ene = math.exp(log_ene)
    p = (ene**2 + 2*ene)**.5
    for j in range(n_theta):
      theta = 2*math.pi*j/n_theta
      p1 = p * math.sin(theta)
      p2 = p * math.cos(theta)
      p_final, part_returns = find_p_final( (p1,p2,0), t, 0.1, fields)
      if part_returns and p_final[0] > 0.:
        ene_final = p_final.dot(p_final) / ( (p_final.dot(p_final) + 1)**.5 + 1)
        if ene_final >= result_array.ene_min:
          theta_fin = math.atan( p_final[1]/ p_final[0] )
          theta_fin *= 180./math.pi
          n_theta_fin = int(theta_fin)
          n_ene_fin = int((math.log(ene_final)-result_array.log_ene_min)/result_array.d_log_ene)
          slice[n_ene_fin, n_theta_fin] = min( slice[n_ene_fin, n_theta_fin], ene)
  return slice

def final_energy_and_angle(t, dt, ene, n_theta, fields, result_array):
  """
  for a given energy and injection time, finds 'all' the final angles and energies of
  particles which start at the surface, by cycling over 'all' initial angles (the only free variable)
  and deposits them to an array as a binary value
  (1 for there exists an initial angle end here in phases space, 0 for there doesn't)
  this would be used to make a contour map (with the exterior loop cycling over time)
  with each contour in the ene-theta space
  being an initial energy (this is my best guess how g.e. kemp made figure 3)
  (though it's not clear why time is cycled exterior to this function, but angle is cycled inside it)
  """
  slice = result_array.data.copy()
  slice *= 0.
  mag_p = (ene**2 + 2*ene)**.5 #magnitude of p
  for i in range(n_theta):
    theta_initial = 2 * math.pi * i/n_theta
    p = numpy.zeros(3)
    p[0] = mag_p * math.sin(theta_initial)
    p[1] = mag_p * math.cos(theta_initial)
    p_final, if_returns, if_enters_vacuum = find_p_final( p, t, dt, fields)
    if if_returns and if_enters_vacuum and p_final[0] > 0:
      ene_final = p_final.dot(p_final) / ( (p_final.dot(p_final) + 1)**.5 + 1)
      if ene_final >= result_array.ene_min:
        theta_fin = math.atan( p_final[1]/ p_final[0] )
        theta_fin *= 180./math.pi
        n_theta_fin = int(theta_fin) + 90
        n_ene_fin = int((math.log10(ene_final)-result_array.log_ene_min)/result_array.d_log_ene)
        slice[n_ene_fin, n_theta_fin] = 1.
  return slice
    

class kemp_ge_fig3_array(object):
  def __init__(self, n_ene):
    self.ene_min = 100./511
    self.ene_max = self.ene_min * 100
    self.log_ene_min = math.log10(self.ene_min)
    self.log_ene_max = math.log10(self.ene_max)
    self.log_ene_range = self.log_ene_max - self.log_ene_min
    self.n_ene = n_ene
    self.data = numpy.zeros( (n_ene, 180))
    self.d_log_ene = self.log_ene_range / n_ene


def kemp_ge_fig_3_test():
  result_array = kemp_ge_fig3_array(100)
  nt = 24 # number of slices in time
  laser_field = reflected_laser_field( 1.5, 1., 45.**.5)
  ene_max = 1.*10**.5/511.
  ene_min = ene_max/10.
  n_ene = 8
  for j in range(n_ene):
    ene = ene_min * 10.**(j*1.0/(n_ene-1))
    for i in range(nt):
      t = 2*math.pi*i/nt
      slice = final_energy_and_angle(t, 0.1, ene, 48, laser_field, result_array)
      result_array.data = numpy.maximum(result_array.data, slice)
  return result_array


def test():
  x = kemp_ge_fig_3_test()
  import matplotlib.pyplot as plt
  plt.figure(figsize=(6,4), dpi = 100)
  plt.axes([0.15, 0.15, 0.75, 0.75])
  plt.xlim([-90,90])
  plt.ylim([-1.0,1.0])
  plt.imshow(x.data, extent=[-90, 90, -1.0, 1.0], aspect = 'auto', origin = 'lower')
  plt.show()

if __name__ == "__main__":
  import os, h5py
  current_directory = os.getcwd()
  max_energy = (100.*10**.5)/511.
  min_energy = (0.1/10**.5)/511.
  slices_per_decade = 4
  decades = math.log10(max_energy/min_energy)
  total_slices = slices_per_decade * decades
  # we won't actually get to the max energy, but that makes more sense I think
  # it's very unclear what energies kemp actually used for the figure
  logged_spacing = (math.log10(max_energy) - math.log10(min_energy))/(slices_per_decade * decades)
  hdf5_file = h5py.File(current_directory + '/geKempfig3_ene_slices.h5', 'w')
  ene_group = hdf5_file.create_group('energy_slices')
  laser_field = reflected_laser_field( 1.5, 1., 1.e12)
  nt = 96
  #loop over all slices selected, saving each one to an hdf5 dataset
  for i in range(int(total_slices)):
    ene = 10.**(math.log10(min_energy) + i * logged_spacing)
    result_array = kemp_ge_fig3_array(100)
    #loop over one full laser cycle in time
    for j in range(nt):
      t =  2*math.pi * j/nt
      slice = final_energy_and_angle(t, 0.1, ene, 96, laser_field, result_array)
      result_array.data = numpy.maximum(result_array.data, slice)
    this_ene_dataset = ene_group.create_dataset( "%08d" % i, data=result_array.data)
    this_ene_dataset.attrs.create('energy', ene)
    print "Finished writing the slice for energy " + str(ene)
  hdf5_file.close()
    
    
  