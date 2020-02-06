#!/usr/bin/env python3

from netCDF4 import Dataset
import numpy as np
import sys

f = Dataset(sys.argv[1], 'r')

h0 = f.variables['h'][0,:]
h1 = f.variables['h'][-1,:]
cos_lat = np.cos(np.radians(f.variables['lat'][:]))

L1 = np.sum(np.sum(np.abs(h1 - h0), axis=1) * cos_lat) / np.sum(np.sum(np.abs(h0), axis=1) * cos_lat)
L2 = np.sqrt(np.sum(np.sum((h1 - h0)**2, axis=1) * cos_lat)) / np.sqrt(np.sum(np.sum(h0**2, axis=1) * cos_lat))
print(L1, L2)