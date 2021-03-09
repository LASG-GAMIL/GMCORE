#!/usr/bin/env python3

import argparse
import cartopy.crs as crs
from cartopy.feature import COASTLINE
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import os
import warnings
warnings.filterwarnings('ignore')

parser = argparse.ArgumentParser('Plot contour for given variable.')
parser.add_argument('-i', dest='input', help='Input h0 data file', nargs='+', required=True)
parser.add_argument('-o', dest='output', help='Output figure file')
parser.add_argument('-v', dest='var', help='Variable name', default='pv')
parser.add_argument('-t', dest='time_step', help='Time step', required=True, type=int)
parser.add_argument('-k', dest='level_idx', help='Vertical level index', type=int)
parser.add_argument('-p', dest='plev', help='Pressure level in hPa', type=float)
parser.add_argument('-c', '--colormap', help='Colormap to use', default='Spectral_r')
parser.add_argument('-l', dest='with_lines', help='Draw contour lines', action='store_true')
parser.add_argument('--no-grid-lines', dest='without_grid_lines', help='Draw grid liens', action='store_true')
parser.add_argument('--coast-lines', dest='with_coast_lines', help='Draw coast lines', action='store_true')
args = parser.parse_args()

ds = xr.open_mfdataset(args.input, concat_dim='time', combine='by_coords')

if not args.var in ds:
	print(f'[Error]: Invalid variable name {args.var}!')
	exit(1)

proj = crs.PlateCarree()

ax = plt.axes(projection=proj)
if not args.without_grid_lines:
	gl = ax.gridlines(crs=proj, draw_labels=True, linewidth=1, color='k', alpha=0.5, linestyle='--')
	gl.xlabels_top = False
	gl.ylabels_right = False
	gl.xformatter = LONGITUDE_FORMATTER
	gl.yformatter = LATITUDE_FORMATTER

cbar_kwargs = {
	'orientation': 'horizontal',
#	'shrink': 0.8
}

var = ds[args.var].isel(time=args.time_step)
if 'lev' in var.dims:
	if args.level_idx != None:
		var = var.isel(lev=args.level_idx)
	else:
		print('[Error]: Please select a vertical level with -k option!')
		exit(1)

var.plot(ax=ax, robust=True, transform=proj, cmap=args.colormap, cbar_kwargs=cbar_kwargs)
if args.with_lines:
	var.plot.contour(ax=ax, transform=proj)

if args.with_coast_lines:
	ax.add_feature(COASTLINE.with_scale('50m'), linewidth=0.5)

if args.output:
	plt.savefig(args.output)
else:
	plt.savefig(f'{os.path.basename(args.input[0]).replace(".nc", "")}.{args.var}.png')
