import numpy, h5py, os

def p1x1_ps(x, p, q, num_par, p1_range, x1_range):
  p1x1 = numpy.zeros((p1_range[2],x1_range[2]))
  dp1 = (p1_range[1] - p1_range[0]) / p1_range[2]
  dx1 = (x1_range[1] - x1_range[0]) / x1_range[2]
  for i in xrange(num_par):
    p1_int = int((p[i,0] - p1_range[0])/dp1)
    x1_int = int((x[i,0] - x1_range[0])/dx1)
    #
    p1_frac = (p[i,0] - p1_range[0])%dp1
    for j in (p1_int, p1_int + 1):
      if j in xrange(p1_range[2]) :
        x1_frac = (x[i,0] - x1_range[0])%dx1
        for k in (x1_int, x1_int + 1):
          if k in xrange(x1_range[2]):
            p1x1[j,k] += q[i] * p1_frac * x1_frac
          x1_frac = 1. - x1_frac
      p1_frac = 1. - p1_frac
  return p1x1


hdf5_data_file = h5py.File('particle_data.h5', 'r')
x = hdf5_data_file.get('x')[:]
p = hdf5_data_file.get('p')[:]
q = hdf5_data_file.get('q')[:]
#
t = hdf5_data_file.attrs.get('TIME')
iter = hdf5_data_file.attrs.get('ITER')
dt = hdf5_data_file.attrs.get('DT')
hdf5_data_file.close()

num_par = x.shape[0]

p1_range = ( -5., 15., 200)
x1_range = (-3.15, 6.30, 200)
p1x1 = p1x1_ps(x, p, q, num_par, p1_range, x1_range)

current_directory = os.getcwd()
#create file
hdf5_p1x1_file = h5py.File(current_directory + '/p1x1.h5', 'w')
#add attributes
hdf5_p1x1_file.attrs.create('NAME', data = numpy.array('py-spec p1x1',dtype = '|S256'))
hdf5_p1x1_file.attrs.create('TYPE', data = numpy.array('grid',dtype = '|S4'))
hdf5_p1x1_file.attrs.create('TIME', data = t)
hdf5_p1x1_file.attrs.create('ITER', data = iter)
hdf5_p1x1_file.attrs.create('DT', data = dt)
hdf5_p1x1_file.attrs.create('TIME UNITS', numpy.array(['1 / \\omega_p'], dtype='|S256'))
# XMIN, XMAX, PERIODIC are of the actual simulation, and are not used anywhere I don't think
hdf5_p1x1_file.attrs.create('XMIN', numpy.array([-1., -1., -1.]))
hdf5_p1x1_file.attrs.create('XMAX', numpy.array([ 1.,  1.,  1.]))
hdf5_p1x1_file.attrs.create('PERIODIC', numpy.array([0, 0, 0], dtype='int32'))
# I don't think MOVE C is used either
hdf5_p1x1_file.attrs.create('MOVE C', numpy.array([0, 0, 0], dtype='int32'))

#create actual dataset
p1x1DS = hdf5_p1x1_file.create_dataset('p1x1', data = p1x1)
#add attributes
p1x1DS.attrs.create('UNITS', numpy.array(['a.u.'], dtype='|S80'))
p1x1DS.attrs.create('LONG_NAME', numpy.array(['p_1x_1'], dtype='|S80'))

#create group for axes (no attributes)
axis_group = hdf5_p1x1_file.create_group('AXIS')

#create first axis dataset
axis1 = numpy.array((x1_range[0], x1_range[1]))
axis1_DS = axis_group.create_dataset('AXIS1', data = axis1)
#add attributes
axis1_DS.attrs.create('TYPE', numpy.array(['linear'], dtype='|S6'))
axis1_DS.attrs.create('UNITS', numpy.array(['c / \\omega_p'], dtype='|S256'))
axis1_DS.attrs.create('NAME', numpy.array(['x1'], dtype='|S256'))
axis1_DS.attrs.create('LONG_NAME', numpy.array(['x_1'], dtype='|S256'))

#create second axis dataset
axis2 = numpy.array((p1_range[0], p1_range[1]))
axis2_DS = axis_group.create_dataset('AXIS2', data = axis2)
#add attributes
axis2_DS.attrs.create('TYPE', numpy.array(['linear'], dtype='|S6'))
axis2_DS.attrs.create('UNITS', numpy.array(['m_e c'], dtype='|S256'))
axis2_DS.attrs.create('NAME', numpy.array(['p1'], dtype='|S256'))
axis2_DS.attrs.create('LONG_NAME', numpy.array(['p_1'], dtype='|S256'))

hdf5_p1x1_file.close()


