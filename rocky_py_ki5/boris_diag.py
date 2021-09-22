import numpy, h5py, os, rocky_py

class ps_axis_param(object):
  def __init__(self, type = 'x', dim = 0, min = 0., max = 1., n = 10):
    self.type = type
    self.dim = dim
    self.min = min
    self.max = max
    self.n = n
  def dx(self):
    return (self.max - self.min)/self.n
  def units(self):
    return {
     'x' : 'c / \\omega_p',
     'p' : 'm_e c',
     }.get(self.type)
  def name(self):
    return self.type + str(self.dim + 1)
  def long_name(self):
    return self.type + '_' + str(self.dim + 1)

def get_2d_ps(part_data, axes_params, weighting_function = False):
  '''
  part_data = {'x':x, 'p':p, 'q':q}
  axes_params = [ps_axis_param0, ps_axis_param1]
  
  weighting_function should expect input x, p, q as arrays
  (wait, or is q a scalar? probably)
  '''
  num_par = part_data['x'].shape[0]
  if weighting_function:
    w = weighting_function
  else:
    def w(*args, **kwargs):
      return 1.
  ax0 = axes_params[0]
  ax1 = axes_params[1]
  shape = []
  for axis in axes_params:
    shape.append(axis.n)
  ps = numpy.zeros(shape)
  data0 = part_data[ax0.type][:,ax0.dim]
  data1 = part_data[ax1.type][:,ax1.dim]
  rdx0 = 1.0/ax0.dx()
  rdx1 = 1.0/ax1.dx()
  for i in xrange(num_par):
    index0 = int((data0[i] - ax0.min)*rdx0)
    index1 = int((data1[i] - ax1.min)*rdx1)
    ww = w(part_data['x'][i], part_data['p'][i], part_data['q'][i])
    frac0 = (data0[i] - ax0.min)%ax0.dx()
    for j in (index0, index0 + 1):
      if j in xrange( ax0.n ):
        frac1 = (data1[i] - ax1.min)%ax1.dx()
        for k in (index1, index1 + 1):
          if k in xrange(ax1.n):
            ps[j,k] += part_data['q'][i]*frac0*frac1*ww
          frac1 = 1. - frac1
      frac0 = 1. - frac0
  return ps

def get_2d_ps_fort(part_data, axes_params, weighting_function = False):
  num_par = part_data['x'].shape[0]
  w = numpy.ones(num_par)
  if weighting_function:
    for i in xrange(num_par):
      w[i] = weighting_function(part_data['x'][i], part_data['p'][i], part_data['q'][i])
  ax0 = axes_params[0]
  ax1 = axes_params[1]
  shape = []
  data0 = part_data[ax0.type][:,ax0.dim]
  data1 = part_data[ax1.type][:,ax1.dim]
  return rocky_py.twod_ps(ax0.n, ax1.n,
   ax0.min, ax0.max, ax1.min, ax1.max,
   data0, data1, w, part_data['q'] )

class visxd_grid_attributes:
  dimensions = 0
  root_group_attrs = {
   'NAME':'species_1',
   'TYPE':'grid',
   'TIME': -1.0,
   'ITER': -1,
   'DT' : -1.0,
   'TIME UNITS':'1 / \\omega_p',
   'XMIN':numpy.array([-1., -1., -1.]),
   'XMAX':numpy.array([ 1.,  1.,  1.]),
   'PERIODIC':numpy.array([0, 0, 0]),
   'MOVE C':numpy.array([0, 0, 0])
   }
  dataset_attrs = {
   'UNITS':'',
   'LONG_NAME':''
   }
  axis_group_attrs = {
   }
  axis1_attrs = {
   'TYPE':'linear',
   'UNITS':'',
   'NAME':'',
   'LONG_NAME':''
   }
  axis2_attrs = {
   'TYPE':'linear',
   'UNITS':'',
   'NAME':'',
   'LONG_NAME':''
   }
  axis3_attrs = {
   'TYPE':'linear',
   'UNITS':'',
   'NAME':'',
   'LONG_NAME':''
   }

def write_2d_grid_file( data, time_dict, axes ):
  
  attributes = visxd_grid_attributes
  att = attributes
  
  attributes.dimensions = 2
  attributes.root_group_attrs['TIME'] = time_dict['TIME']
  attributes.root_group_attrs['ITER'] = time_dict['ITER']
  attributes.root_group_attrs['DT'] = time_dict['DT']
  
  p0 = axes[0]
  p1 = axes[1]

  #reverse (and offset by one, of course)
  #(of course)
  att.axis1_attrs['UNITS'] = p1.units()
  att.axis1_attrs['NAME'] = p1.name()
  att.axis1_attrs['LONG_NAME'] = p1.long_name()
  
  att.axis2_attrs['UNITS'] = p0.units()
  att.axis2_attrs['NAME'] = p0.name()
  att.axis2_attrs['LONG_NAME'] = p0.long_name()
  
  att.dataset_attrs['LONG_NAME'] = p0.long_name() + p1.long_name()
  
  save_dir = ''
  save_dir += os.getcwd() + '/'
  save_dir += 'MS/PHA' + '/'
  #re...reverse
  ps_name = attributes.axis2_attrs['NAME'] + attributes.axis1_attrs['NAME']
  save_dir += ps_name + '/'
  save_dir += attributes.root_group_attrs['NAME'] + '/'
  try:
    os.makedirs( save_dir )
  except OSError:
    pass
  file_name = ps_name + '_' + "%06d" % attributes.root_group_attrs['ITER'] + '.h5'
  data_file = h5py.File( save_dir + file_name, 'w' )
  for pair in attributes.root_group_attrs.items():
    x, y = pair
    data_file.attrs.create(x,y)
  ps_ds = data_file.create_dataset( ps_name, data = data )
  for pair in attributes.dataset_attrs.items():
    x, y = pair
    ps_ds.attrs.create(x,y)
  axis_group = data_file.create_group('AXIS')
  axis1_ds = axis_group.create_dataset('AXIS1', data = numpy.array((axes[1].min, axes[1].max)))
  for pair in attributes.axis1_attrs.items():
    x, y = pair
    axis1_ds.attrs.create(x,y)
  axis2_ds = axis_group.create_dataset('AXIS2', data = numpy.array((axes[0].min, axes[0].max)))
  for pair in attributes.axis2_attrs.items():
    x, y = pair
    axis2_ds.attrs.create(x,y)
  data_file.close()

if __name__ == '__main__':
  
  hdf5_data_file = h5py.File('particle_data.h5', 'r')
  x = hdf5_data_file.get('x')[:]
  p = hdf5_data_file.get('p')[:]
  q = hdf5_data_file.get('q')[:]
  #
  t = hdf5_data_file.attrs.get('TIME')
  iter = hdf5_data_file.attrs.get('ITER')
  dt = hdf5_data_file.attrs.get('DT')
  hdf5_data_file.close()
  
  ps_axis_param0 = ps_axis_param('p', dim = 0, min = -5., max = 15., n = 200)
  p0 = ps_axis_param0
  ps_axis_param1 = ps_axis_param('x', 0, -3.15, 12.6-3.15, 200)
  p1 = ps_axis_param1
  axes_params = [p0, p1]
  
  part_dict = {'x':x, 'p':p, 'q':q}
  def w_func(x, p, q):
    import math
    theta = math.atan2(x[0],p[0])
    r = math.sqrt(x[0]**2 + p[0]**2)
    return math.sin(10*(theta + r))
    
  ps_ds = get_2d_ps( part_dict, axes_params, w_func )
  
  time_dict = {'TIME':t, 'ITER':iter, 'DT':dt}
  write_2d_grid_file( ps_ds, time_dict, axes_params )
  