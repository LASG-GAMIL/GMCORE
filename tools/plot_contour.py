#!/usr/bin/env python3

import argparse
import cartopy.crs as crs
from cartopy.io.shapereader import Reader as ShapeReader
from cartopy.feature import ShapelyFeature
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import os
import numpy as np
import pendulum

parser = argparse.ArgumentParser('Plot contour for given variable.')
parser.add_argument('-i', dest='input', help='Input h0 data file', nargs='+', required=True)
parser.add_argument('-o', dest='output', help='Output figure file')
parser.add_argument('-v', dest='var', help='Variable name', default='pv')
parser.add_argument('-t', dest='time_step', help='Time step', required=True, type=int)
parser.add_argument('-l', dest='with_lines', help='Draw contour lines', action='store_true')
args = parser.parse_args()

if not args.output:
  args.output = args.var + '.png'

f = []
for file_path in args.input:
  if not os.path.isfile(file_path):
    print(f'[Error]: File {file_path} does not exist!')
    exit(1)
  f.append(Dataset(file_path, 'r'))

proj = crs.PlateCarree()

ax = plt.axes(projection=proj)

ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)

def parse_time(time_var, time_step):
  dt, base_time_str = time_var.units.split('since')
  base_time = pendulum.from_format(base_time_str.strip(), 'YYYY-MM-DDTHH_mm_ss')
  if dt == 'hours ':
    return base_time.add(hours=time_var[time_step])

def plot(lon, lat, var):
  plt.pcolormesh(lon, lat, var, transform=proj, cmap='rainbow')
  plt.colorbar(orientation='horizontal')

  if args.with_lines:
    xlon, xlat = np.meshgrid(lon, lat)
    plt.contour(xlon, xlat, var, transform=proj)

  plt.title(f'{f[i].variables[args.var].long_name} @ {parse_time(f[0].variables["time"], args.time_step)}')

for i in range(len(f)):
  plot(f[i].variables['lon'][:], f[i].variables['lat'][:], f[i].variables[args.var][args.time_step,:,:])

plt.savefig(args.output)
plt.close()
