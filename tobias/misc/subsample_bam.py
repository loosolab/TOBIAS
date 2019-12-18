#!/usr/bin/env python

"""
SubsampleBam: Samples from bam in percentage-steps with replicates per step  

@author: Mette Bentsen
@contact: mette.bentsen (at) mpi-bn.mpg.de
@license: MIT

"""

import os
import sys
import multiprocessing as mp
import subprocess
import argparse

from tobias.parsers import add_subsample_arguments
from tobias.utils.utilities import *
from tobias.utils.logger import *

#-------------------------------------------------------------------#
def run_commandline(command):

	print("{0} RUNNING: \"{1}\"".format(datetime.now(), command))
	print(command.split())
	p = subprocess.call(command.split())

	return(1)

#-------------------------------------------------------------------#
def run_subsampling(args):

	check_required(args, ["bam"])

	args.prefix = os.path.splitext(os.path.basename(args.bam))[0] if args.prefix == None else args.prefix
	args.outdir = os.path.abspath(args.outdir) if args.outdir != None else os.path.abspath(os.getcwd())

	#---------------------------------------------------#

	logger = TobiasLogger("SubsampleBam", args.verbosity)
	logger.begin()

	#### Getting ready for running 
	cmd_calls = []
	for frac in range(args.start,args.end,args.step):
		for rand in range(1,args.no_rand+1):
			outfile = os.path.join(args.outdir, args.prefix + "_{0}_r{1}.bam".format(frac, rand))
			call = "samtools view -s {0} -bh -o {1} {2}".format((frac)/float(100)+rand, outfile, args.bam)
			cmd_calls.append(call)

	#### Run tasks ###
	all_tasks = len(cmd_calls)
	complete_tasks = 0

	if args.cores > 1:

		prev_complete_tasks = 0

		#Send all jobs to pool
		pool = mp.Pool(processes=args.cores)
		task_list = [pool.apply_async(run_commandline, cmd) for cmd in cmd_calls]
		pool.close()
		monitor_progress(task_list, logger)
		pool.join()

	else:
		for command in cmd_calls:
			res = run_commandline(command)
			complete_tasks += 1
			print(complete_tasks / float(all_tasks) * 100.0)

	logger.end()

#--------------------------------------------------------------------------------------------------------#
if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser = add_subsample_arguments(parser)
	args = parser.parse_args()

	if len(sys.argv[1:]) == 0:
		parser.print_help()
		sys.exit()

	run_subsampling(args)
