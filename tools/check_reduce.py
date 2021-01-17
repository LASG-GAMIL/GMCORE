#!/usr/bin/env python3

import numpy as np
import f90nml
import sys
import matplotlib.pyplot as plt

namelist = f90nml.read(sys.argv[1])

num_full_lat = namelist['gmcore_control']['num_lat']
num_half_lat = num_full_lat - 1
dlon         = np.radians(2 * np.pi / namelist['gmcore_control']['num_lon'])
dlat0        = dlon
pole_mul     = namelist['gmcore_control']['coarse_pole_mul']
pole_decay   = namelist['gmcore_control']['coarse_pole_decay']

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

reduce_factors = np.ndarray([num_half_lat])
reduce_factors.fill(1)
for j in range(num_half_lat):
	if half_lat[j] <= 0:
		jr = j
	else:
		jr = num_half_lat - 1 - j
	if jr < len(namelist['gmcore_control']['reduce_factors']):
		reduce_factors[j] = namelist['gmcore_control']['reduce_factors'][jr]

le_lat = np.ndarray([num_half_lat])
for j in range(num_half_lat):
	le_lat[j] = np.cos(half_lat[j]) * dlon * reduce_factors[j]

fig = plt.figure(figsize=(10, 5))
ax = fig.add_subplot(1, 1, 1)
plt.plot(np.degrees(half_lat), np.degrees(le_lat), marker='o')
plt.xticks(np.arange(-90, 90, 5))
plt.show()
plt.close()
