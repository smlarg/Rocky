import h5py, os, numpy
import matplotlib.pyplot as plt
hdf5_file = h5py.File('geKempfig3_ene_slices.h5','r')
slices_group = hdf5_file.get('energy_slices')
keys = slices_group.keys()
keys.sort()
summed_data_set = slices_group.get(keys[0])[:]
summed_data_set *= 0.
total_data_set = summed_data_set.copy()
total_data_set += 0.
for i in range(4):
  summed_data_set *= 0.
  for j in range(4*i, 4*i + 4):
    summed_data_set = numpy.maximum( summed_data_set, slices_group.get(keys[j])[:])
  summed_data_set *= (4-i)
  total_data_set = numpy.maximum( total_data_set, summed_data_set)

total_data_set *= -1
total_data_set += 4
plt.figure(figsize=(6,4), dpi = 100)
plt.axes([0.15, 0.15, 0.75, 0.75])
plt.xlim([-90,90])
plt.ylim([-1.0,1.0])
plt.xlabel(r'Degrees', size = 16)
plt.ylabel(r'Final Energy (log)', size = 16)
plt.title('No Skin Layer')
plt.imshow(total_data_set, extent=[-90, 90, -1.0, 1.0], aspect = 'auto', origin = 'lower',
  vmin = 0., vmax = 4.)
cb = plt.colorbar()
cb.set_label('Initial Temp (log)')
plt.savefig('unskinned.eps')
plt.close()
