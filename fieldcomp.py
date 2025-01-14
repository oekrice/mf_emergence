from scipy.io import netcdf_file
import numpy as np
import matplotlib.pyplot as plt
#Script for comparing two magnetic fields

runs = [0,9]

snap = 80

comps = []
for ri, run in enumerate(runs):
    fname = '/nobackup/trcn27/mf3d0/%03d/%04d.nc' % (run, snap)

    data = netcdf_file(fname, 'r', mmap=False)

    bz = np.swapaxes(data.variables['bz'][:],0,2)


    data.close()

    comps.append(bz[:,:,0])

comps = np.array(comps)

diff = comps[1] - comps[0]

fig, axs = plt.subplots(1, 3, figsize = (10,3.5))
axs[0].pcolormesh(comps[0])
axs[1].pcolormesh(comps[1])
axs[2].pcolormesh(diff)

print('error', np.max(np.abs(diff)))


plt.show()
