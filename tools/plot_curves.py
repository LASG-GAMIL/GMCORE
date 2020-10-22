#!/usr/bin/env python3

import argparse
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import xarray as xr
import os
import warnings

warnings.filterwarnings('ignore')

parser = argparse.ArgumentParser(description='Plot curves from model data.')
parser.add_argument('-i', '--input', help='Model data', required=True)
parser.add_argument('-o', '--output', help='Output figure file name')
args = parser.parse_args()

if not args.output:
	args.output = os.path.basename(args.input).replace('.nc', '') + '.te_tpe.pdf'

ds = xr.open_dataset(args.input)

te = ds['te'][:]
tpe = ds['tpe'][:]

time = ds['time']

#plt.figure(figsize=(8, 4))

ax1 = plt.axes()
ax2 = ax1.twinx()

color1 = 'red'
color2 = 'blue'

te.plot(ax=ax1, color=color1, linewidth=3)
ax1.tick_params(axis='y', colors=color1)
ax1.set_ylabel('Total energy', color=color1)
tpe.plot(ax=ax2, color=color2, linewidth=3)
ax2.tick_params(axis='y', colors=color2)
ax2.set_ylabel('Total potential enstrophy', color=color2)

plt.savefig(args.output)

print(f'[Notice]: Output {args.output}.')
