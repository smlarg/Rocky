if __name__ == "__main__":
  import h5py, os
  import matplotlib.pyplot as plt
  
  hdf5_file = h5py.File('p1p2t.h5','r')
  p1_max = hdf5_file.attrs.get('p1_max')
  p2_max = hdf5_file.attrs.get('p2_max')
  #
  energy_dataset = hdf5_file.get('reetry_energy')
  ene_keys = energy_dataset.keys()
  ene_keys.sort()
  energy_slices = []
  times = []
  iters = []
  for key in ene_keys:
    energy_slices.append(energy_dataset.get(key)[...])
    times.append(energy_dataset.get(key).attrs.get('time'))
    iters.append(energy_dataset.get(key).attrs.get('iteration'))
  zipped_up = zip(energy_slices, times, iters)
  try:
    os.mkdir('energy_frames')
  except:
    pass
  maximum = 0.
  for ene in energy_slices:
    maximum = max(maximum, ene.max())
  for pants in zipped_up:
    ds = pants[0]
    t = pants[1]
    i = pants[2]
    plt.figure(figsize=(6,4), dpi = 100)
    plt.axes([0.15, 0.15, 0.75, 0.75])
    plt.xlim([-p2_max,p2_max])
    plt.ylim([-p1_max,0.0])
    plt.title('time = ' + str(t) )
    plt.imshow(ds, cmap=plt.cm.hot, vmin=0., vmax=maximum, extent=[-p2_max, p2_max, -p1_max,0.])
    plt.colorbar()
    #plt.imshow(Y, extent=[-p2_max, p2_max, -p1_max,0.])
    index = "%03d" % i
    plt.savefig( 'energy_frames' + "/" + "ene_frame" + index + ".jpg")
    plt.close()
  #
  time_dataset = hdf5_file.get('reetry_time')
  time_keys = time_dataset.keys()
  time_keys.sort()
  time_slices = []
  times = []
  iters = []
  for key in time_keys:
    time_slices.append(time_dataset.get(key)[...])
    times.append(time_dataset.get(key).attrs.get('time'))
    iters.append(time_dataset.get(key).attrs.get('iteration'))
  zipped_up = zip(time_slices, times, iters)
  try:
    os.mkdir('time_frames')
  except:
    pass
  maximum = 0.
  for time in time_slices:
    maximum = max(maximum, time.max())
  for pants in zipped_up:
    ds = pants[0]
    t = pants[1]
    i = pants[2]
    plt.figure(figsize=(6,4), dpi = 100)
    plt.axes([0.15, 0.15, 0.75, 0.75])
    plt.xlim([-p2_max,p2_max])
    plt.ylim([-p1_max,0.0])
    plt.title('time = ' + str(t) )
    plt.imshow(ds, cmap=plt.cm.hot, vmin=0., vmax=maximum, extent=[-p2_max, p2_max, -p1_max,0.])
    plt.colorbar()
    #plt.imshow(Y, extent=[-p2_max, p2_max, -p1_max,0.])
    index = "%03d" % i
    plt.savefig( 'time_frames' + "/" + "time_frame" + index + ".jpg")
    plt.close()
  #
  #
  hdf5_file.close()
