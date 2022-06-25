#!/usr/bin/env python3

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import sys

file_path = sys.argv[1]

f = xr.open_dataset(file_path, decode_cf=False)

fig = plt.figure(figsize=(10, 8))

ax0 = plt.subplot(121, projection=ccrs.SouthPolarStereo())
ax0.set_extent([-180, 180, -90, -70], ccrs.PlateCarree())
ax0.contourf(f.lon, f.lat, f.q1[9,0,:,:], transform=ccrs.PlateCarree(), cmap='rainbow')

ax1 = plt.subplot(122, projection=ccrs.PlateCarree(central_longitude=180))
ax1.contourf(f.lon, f.lat, f.q1[12,0,:,:], transform=ccrs.PlateCarree(), cmap='rainbow')
gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='k', alpha=0.5, linestyle='--')
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER

plt.show()
