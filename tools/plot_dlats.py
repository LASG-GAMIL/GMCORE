#!/usr/bin/env python3

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
import sys

if len(sys.argv) > 1:
	f = Dataset(sys.argv[1], 'r')
	num_full_lat = f.variables['lat'].size
	num_half_lat = f.variables['ilat'].size
	full_lat = np.radians(f.variables['lat'])
	half_lat = np.radians(f.variables['ilat'])
	dlat = np.diff(full_lat)
else:
	num_full_lat = 340
	num_half_lat = num_full_lat - 1
	dlat0        = np.radians(0.5)
	pole_mul     = 2
	pole_decay   = 100

	full_lat = np.ndarray([num_full_lat])
	half_lat = np.ndarray([num_half_lat])
	dlat = np.ndarray([num_half_lat])

	dlat1 = np.pi / num_half_lat
	for j in range(num_half_lat):
		half_lat[j] = -np.pi / 2 + (j + 0.5) * dlat1

	for j in range(num_half_lat):
		dlat[j] = dlat0 * (1 + (pole_mul - 1) * np.exp(-pole_decay * (np.abs(half_lat[j]) - np.pi / 2)**2))
	dlat = dlat / np.sum(dlat) * np.pi

	half_lat[0] = -np.pi / 2 + dlat[0] / 2
	for j in range(1, num_half_lat):
		half_lat[j] = half_lat[j-1] + (dlat[j-1] + dlat[j]) / 2

fig = plt.figure(figsize=(10, 5))
ax = fig.add_subplot(1, 1, 1)
ax.set_title(f'NLAT = {num_full_lat}')
plt.plot(np.degrees(half_lat), np.degrees(dlat))
if 'dlat0' in locals():
	dlat[:] = dlat0
	plt.plot(np.degrees(half_lat), np.degrees(dlat))
plt.show()
plt.close()
