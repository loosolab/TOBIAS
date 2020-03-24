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

	#Split commands on ;
	for single_command in command.split(";"):
		try:
			print("{0} RUNNING: \"{1}\"".format(datetime.now(), single_command))
			p = subprocess.call(single_command.split())
		except:
			sys.exit("ERROR running: {0}".format(single_command))
	return(1)

#-------------------------------------------------------------------#
def run_subsampling(args):

	check_required(args, ["bam"])
	check_files([args.bam])

	args.bam = os.path.abspath(args.bam)
	args.prefix = os.path.splitext(os.path.basename(args.bam))[0] if args.prefix == None else args.prefix
	args.outdir = os.path.abspath(args.outdir) if args.outdir != None else os.path.abspath(os.getcwd())

	make_directory(args.outdir)

	#---------------------------------------------------#

	logger = TobiasLogger("SubsampleBam", args.verbosity)
	logger.begin()

	#---------------------------------------------------#

	#Check if samtools is available
	try:
		p = subprocess.call("samtools", stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL) 	#stdout to /dev/null
	except:
		logger.error("samtools is not available. Please install samtools to be able to use SubsampleBam.")
		sys.exit(1)

	#### Getting ready for running 
	cmd_calls = []
	for frac in range(args.start, args.end+1, args.step):	#+1 ensures that args.end is included
		for rand in range(1, args.no_rand+1):
			outfile = os.path.join(args.outdir, args.prefix + "_{0}_r{1}.bam".format(frac, rand))
			call = "samtools view -s {0} -bh -o {1} {2};samtools index {1}".format((frac)/float(100)+rand, outfile, args.bam)
			
			#Check if outfile already exists
			if args.force or not os.path.exists(outfile):
				logger.debug("Adding call to create {0}".format(outfile))
				cmd_calls.append(call)
			else:
				logger.warning("Output {0} already exists - skipping! Use --force to force creation of subsamples.".format(outfile))	

	#### Run tasks ###
	n_tasks = len(cmd_calls)
	logger.debug("Created {0} tasks".format(n_tasks))

	complete_tasks = 0
	if args.cores > 1:
		logger.debug("Running tasks using multiprocessing ({0} cores)".format(args.cores))
		prev_complete_tasks = 0

		#Send all jobs to pool
		pool = mp.Pool(processes=args.cores)
		task_list = [pool.apply_async(run_commandline, (cmd,)) for cmd in cmd_calls]
		pool.close()
		monitor_progress(task_list, logger)
		pool.join()

	else:
		for command in cmd_calls:
			res = run_commandline(command)
			complete_tasks += 1
			logger.info("- Progress: {0:.1f}%".format(round(complete_tasks / float(n_tasks) * 100.0, 3)))

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
