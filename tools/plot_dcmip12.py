#!/usr/bin/env python3

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import sys

file_path = sys.argv[1]

var_name = 'q1'
z_name = 'z'
itime = 24
interp_z = np.linspace(0, 12000, 61)
levels = np.linspace(0, 1.1, 12)

f = xr.open_dataset(file_path, decode_cf=False)

lat = f.lat
z = f[z_name].mean('lon')

tmp = f[var_name].mean(f[var_name].dims[3])
var1 = np.ndarray((interp_z.size,lat.size))
var2 = np.ndarray((interp_z.size,lat.size))
for i in range(lat.size):
	var1[:,i] = np.interp(interp_z, z[12,::-1,i], tmp[12,::-1,i])
	var2[:,i] = np.interp(interp_z, z[24,::-1,i], tmp[24,::-1,i])
z = np.repeat(interp_z, lat.size).reshape((interp_z.size, lat.size))

fig, ax = plt.subplots(2, 1, figsize=(10, 8))

fig.subplots_adjust(right=0.8, wspace=0.02, hspace=0.1)

ax[0] = plt.subplot(211)
ax[0].set_xticks([])
ax[1] = plt.subplot(212)
x = np.tile(lat, (z.shape[0],1))

ax[0].contourf(x, z, var1, cmap='rainbow', levels=levels)
im = ax[1].contourf(x, z, var2, cmap='rainbow', levels=levels)

fig.colorbar(im, cax=fig.add_axes([0.83, 0.1, 0.02, 0.8]))

plt.show()
