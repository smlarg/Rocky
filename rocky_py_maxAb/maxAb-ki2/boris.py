import numpy, math, time, rocky_py
from math import *
from random import random
r = lambda : random() - 0.5
import boris_diag as bd

def init_plasma(p_max, dp, x_min, x_max, dx, Temp):
  '''
  Initialize a 1D2V Maxwellian plasma
  '''
  def dist_func(p, p_max, T):
    if p.dot(p) > p_max**2:
      return 0.
    return math.exp(-p.dot(p)/T)
  
  nx = int( (x_max - x_min) / dx)  
  np = int( 2 * p_max/ dp)
  
  x = numpy.zeros((np * np * nx, 3))
  p = numpy.zeros((np * np * nx, 3))
  q = numpy.zeros((np * np * nx))
  
  index = 0
  for i in range(nx):
    for j in range(np):
      for k in range(np):
        x_i = numpy.array(((i+r()+0.5)*dx + x_min,0.,0.))
        p_i = numpy.array( (j+r(),k+r(),0.5*np) )
        p_i *= dp
        p_i -= p_max
        q_i = dist_func(p_i, p_max, Temp)
        if q_i > 0.:
          x[index] = x_i
          p[index] = p_i
          q[index] = q_i
          index += 1
  
  x = x[:index]
  p = p[:index]
  q = q[:index]
  
  return x, p, q

if __name__ == "__main__":
  
  x, p, q = init_plasma(p_max = 3.5, dp = 0.05,
   x_min = 0.0, x_max = 25.2, dx = 0.05, Temp = 0.548)
  
  def vector_potential(x, t):
    if x < 0:
      return -1.65835921350013*sin(1.0*t + 1.0*x) + 6.0*cos(1.0*t - 1.0*x) - 3.31671842700025*cos(1.0*t + 1.0*x)
    else:
      return (-1.65835921350013*sin(1.0*t - 2.23606797749979*x) + 2.68328157299975*cos(1.0*t - 2.23606797749979*x))*exp(-2.0*x)
  
  for i in xrange(x.shape[0]):
    p[i][1] -= vector_potential(x[i][0], 0)
  
  t_max = 25.2
  dt = 0.01
  t = 0.
  iter = 0
  
  diag_n = 25
  p1x1_axes_params = [bd.ps_axis_param('p', dim = 0, min = -5., max = 15., n = 200),
   bd.ps_axis_param('x', 0, -3.15, 12.6-3.15, 200)]
  p1p2_axes_params = [bd.ps_axis_param('p', dim = 0, min = -5., max = 15., n = 200),
   bd.ps_axis_param('p', 1, -12.5, 12.5, 250)]
  def p1p2_w_func(x, *args):
    return x[0] < 12.6-3.15
  def p1p2_fake_func(x, p, *args):
    theta = math.atan2(p[0], p[1])
    r = math.sqrt(p[0]**2 + p[1]**2)
    return math.sin(r+theta)
  
  time0 = time.clock()
  
  while (t < t_max):
    print "at time = " + str(t)
    x, p = rocky_py.pushall(x, p, t, dt)
    if not iter%diag_n:
      print 'running diagnostics'
      part_dict = {'x':x, 'p':p, 'q':q}
      time_dict = {'TIME':t, 'ITER':iter, 'DT':dt}
      bd.write_2d_grid_file( bd.get_2d_ps_fort( part_dict, p1x1_axes_params ),
       time_dict, p1x1_axes_params )
      bd.write_2d_grid_file( bd.get_2d_ps_fort( part_dict, p1p2_axes_params, p1p2_w_func ),
       time_dict, p1p2_axes_params )
    time1 = time0
    time0 = time.clock()
    print "timestep length = " + str(time0 - time1) + "s"
    t += dt 
    iter += 1
