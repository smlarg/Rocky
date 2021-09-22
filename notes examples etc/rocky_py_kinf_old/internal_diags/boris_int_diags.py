import numpy, math, time, rocky_py
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
        p_i = numpy.array( (j+r(),k+r(),0.5*np) )
        x_i = numpy.array(((i+r()+0.5)*dx + x_min,0.,0.))
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
  
  t_max = 25.2
  dt = 0.01
  t = 0.
  iter = 0
  
  diag_n = 251
  p1x1_axes_params = [bd.ps_axis_param('p', dim = 0, min = -5., max = 15., n = 200),
   bd.ps_axis_param('x', 0, -3.15, 12.6-3.15, 200)]
  p1p2_axes_params = [bd.ps_axis_param('p', dim = 0, min = -5., max = 15., n = 200),
   bd.ps_axis_param('p', 1, -10., 10., 200)]
  def p1p2_w_func(x, *args):
    return x[0] < 12.6-3.15
  
  time0 = time.clock()
  
  while (t < t_max):
    print "at time = " + str(t)
    x, p = rocky_py.pushall(x, p, t, dt)
    time1 = time0
    time0 = time.clock()
    print "timestep length = " + str(time0 - time1) + "s"
    if not iter%diag_n:
      part_dict = {'x':x, 'p':p, 'q':q}
      time_dict = {'TIME':t, 'ITER':iter, 'DT':dt}
      bd.write_2d_grid_file( bd.get_2d_ps( part_dict, p1x1_axes_params ),
       time_dict, p1x1_axes_params )
      bd.write_2d_grid_file( bd.get_2d_ps( part_dict, p1p2_axes_params, p1p2_w_func ),
       time_dict, p1p2_axes_params )
    t += dt 
    iter += 1
