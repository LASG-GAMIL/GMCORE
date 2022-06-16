#!/usr/bin/env python3

import argparse
import cdsapi
import pendulum
import os

parser = argparse.ArgumentParser(description='Get ERA5 reanalysis data for creating initial data.')
parser.add_argument('-t', dest='time', help='Time (YYYYMMDDHH)', required=True)
parser.add_argument('--res', help='Resoultion (1p0, 0p25)', default='0p25')
args = parser.parse_args()

args.time = pendulum.from_format(args.time, 'YYYYMMDDHH')

plev_file = f'era5_{args.time.format("YYYYMMDDHH")}_plev.nc'
sfc_file  = f'era5_{args.time.format("YYYYMMDDHH")}_sfc.nc'

grids = { '1p0': '1.0/1.0', '0p25': '0.25/0.25' }

c = cdsapi.Client()

if not os.path.isfile(plev_file):
	c.retrieve(
		'reanalysis-era5-pressure-levels',
		{
			'product_type': 'reanalysis',
			'grid': grids[args.res],
			'format': 'netcdf',
			'variable': [
				'temperature',
				'u_component_of_wind',
				'v_component_of_wind',
				'specific_humidity'
			],
			'pressure_level': [
				'1'  , '2'  , '3'  ,
				'5'  , '7'  , '10' ,
				'20' , '30' , '50' ,
				'70' , '100', '125',
				'150', '175', '200',
				'225', '250', '300',
				'350', '400', '450',
				'500', '550', '600',
				'650', '700', '750',
				'775', '800', '825',
				'850', '875', '900',
				'925', '950', '975',
				'1000',
			],
			'year': args.time.format('YYYY'),
			'month': args.time.format('MM'),
			'day': args.time.format('DD'),
			'time': args.time.format('HH:00'),
		},
		plev_file
	)

if not os.path.isfile(sfc_file):
	c.retrieve(
		'reanalysis-era5-single-levels',
		{
			'product_type': 'reanalysis',
			'grid': grids[args.res],
			'format': 'netcdf',
			'variable': [
				'mean_sea_level_pressure',
			],
			'year': args.time.format('YYYY'),
			'month': args.time.format('MM'),
			'day': args.time.format('DD'),
			'time': args.time.format('HH:00'),
		},
		sfc_file
	)

out_file = f'era5_{args.time.format("YYYYMMDDHH")}.nc'
if not os.path.isfile(out_file):
	os.system(f'cdo merge {plev_file} {sfc_file} {out_file}')
	os.system(f'rm -f {plev_file} {sfc_file}')
