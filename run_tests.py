#!/usr/bin/env python3

import argparse
import os
import subprocess
import sys

parser = argparse.ArgumentParser('Run GMCORE tests.')
parser.add_argument('--slurm', help='Use SLURM job manager', action='store_true')
parser.add_argument('-q', '--queue', help='Job queue')
parser.add_argument('-n', '--np', help='Processes to use for running tests', type=int, default=2)
parser.add_argument('-p', '--ntasks-per-node', type=int, default=20)
parser.add_argument('-m', '--node-list')
parser.add_argument('-w', '--work-root', help='Where to run tests', required=True)
parser.add_argument('-c', '--cases', help='Which cases to run', nargs='+', default=[])
args = parser.parse_args()

gmcore_root = os.path.dirname(os.path.realpath(__file__))

if args.node_list:
	args.node_list = os.path.abspath(args.node_list)

if not os.path.isdir(args.work_root):
	os.makedirs(args.work_root)

os.chdir(args.work_root)
testbed_root = os.getcwd() + '/GMCORE-TESTBED/'

def run(cmd):
	print(f'==> {cmd}')
	res = subprocess.run(cmd, shell=True, check=True)

def mpiexec(exe, namelist, args):
	if args.slurm:
		if not args.queue:
			print('[Error]: No job queue is provided!')
			exit(1)
		run(f'srun -p {args.queue} -n {args.np} --ntasks-per-node {args.ntasks_per_node} --mpi=pmi2 --exclusive {gmcore_root}/build/{exe} {namelist}')
	elif args.node_list:
		run(f'mpiexec -rr -f {args.node_list} -np {args.np} {gmcore_root}/build/{exe} {namelist}')
	else:
		run(f'mpiexec -np {args.np} {gmcore_root}/build/{exe} {namelist}')

def run_case(exe, case):
	if len(args.cases) == 0 or case in args.cases:
		os.chdir(testbed_root + case)
		mpiexec(exe, 'namelist', args)

if not os.path.isdir('GMCORE-TESTBED'):
	run('git clone https://gitee.com/dongli85/GMCORE-TESTBED')

run_case('gmcore_swm_driver.exe', 'swm.rh.180x90')
run_case('gmcore_swm_driver.exe', 'swm.rh.360x180')
run_case('gmcore_swm_driver.exe', 'swm.mz.180x90')
run_case('gmcore_swm_driver.exe', 'swm.mz.360x180')
run_case('gmcore_driver.exe', 'rh.180x90')
run_case('gmcore_driver.exe', 'rh.360x180')
run_case('gmcore_driver.exe', 'mz.180x90')
run_case('gmcore_driver.exe', 'mz.360x180')
run_case('gmcore_driver.exe', 'bw.180x90')
run_case('gmcore_driver.exe', 'bw.360x180')
