#!/usr/bin/env python

"""
ScoreBigwig: Calculate footprint tracks from cutsite bigwig

@author: Mette Bentsen
@contact: mette.bentsen (at) mpi-bn.mpg.de
@license: MIT

"""

import os
import sys
import argparse
import numpy as np
import math
import textwrap
import logging

import pyBigWig
import multiprocessing as mp
from scipy import stats

#Internal functions and classes
from tobias.parsers import add_scorebigwig_arguments
from tobias.utils.utilities import *
from tobias.utils.regions import OneRegion, RegionList
from tobias.utils.sequences import *
from tobias.utils.signals import *
from tobias.utils.logger import TobiasLogger

#--------------------------------------------------------------------------------#
def calculate_scores(regions, args):

	logger = TobiasLogger("", args.verbosity, args.log_q)

	pybw_signal = pyBigWig.open(args.signal) 	#cutsites signal
	pybw_header = pybw_signal.chroms()			
	chrom_lengths = {chrom: int(pybw_header[chrom]) for chrom in pybw_header}

	#Set flank to enable scoring in ends of regions
	flank = args.region_flank

	#Go through each region
	for i, region in enumerate(regions):

		logger.debug("Calculating scores for region: {0}".format(region))

		#Extend region with necessary flank
		region.extend_reg(flank)
		reg_key = (region.chrom, region.start+flank, region.end-flank)	#output region

		#Get bigwig signal in region
		signal = region.get_signal(pybw_signal, logger=logger)
		signal = np.nan_to_num(signal).astype("float64")

		#-------- Prepare signal for score calculation -------#
		if args.absolute:
			signal = np.abs(signal)

		if args.min_limit != None:
			signal[signal < args.min_limit] = args.min_limit
		if args.max_limit != None:
			signal[signal > args.max_limit] = args.max_limit

		#------------------ Calculate scores ----------------#
		if args.score == "sum":
			scores = fast_rolling_math(signal, args.window, "sum")

		elif args.score == "mean":
			scores = fast_rolling_math(signal, args.window, "mean")

		elif args.score == "footprint":
			scores = tobias_footprint_array(signal, args.flank_min, args.flank_max, args.fp_min, args.fp_max)		#numpy array

		elif args.score == "FOS":
			scores = FOS_score(signal, args.flank_min, args.flank_max, args.fp_min, args.fp_max)
			#scores = -scores

		elif args.score == "none":
			scores = signal
		
		else:
			sys.exit("Scoring {0} not found".format(args.score))
		
		#----------------- Post-process scores --------------#
		
		#Smooth signal with args.smooth bp
		if args.smooth > 1:
			scores = fast_rolling_math(scores, args.smooth, "mean")

		#Remove ends to prevent overlap with other regions
		if flank > 0:
			scores = scores[flank:-flank]

		args.writer_qs["scores"].put(("scores", reg_key, scores))

	return(1)

#------------------------------------------------------------------------------------------#

def run_scorebigwig(args):
	
	check_required(args, ["signal", "output", "regions"])
	check_files([args.signal, args.regions], "r")
	check_files([args.output], "w")

	#---------------------------------------------------------------------------------------#
	# Create logger and write info to log
	#---------------------------------------------------------------------------------------#

	logger = TobiasLogger("ScoreBigwig", args.verbosity)
	logger.begin()

	parser = add_scorebigwig_arguments(argparse.ArgumentParser())
	logger.arguments_overview(parser, args)
	logger.output_files([args.output])

	logger.debug("Setting up listener for log")
	logger.start_logger_queue()
	args.log_q = logger.queue

	#---------------------------------------------------------------------------------------#
	#----------------------- I/O - get regions/bigwig ready --------------------------------#
	#---------------------------------------------------------------------------------------#

	logger.info("Processing input files")

	logger.info("- Opening input cutsite bigwig")
	pybw_signal = pyBigWig.open(args.signal)
	pybw_header = pybw_signal.chroms()
	chrom_info = {chrom:int(pybw_header[chrom]) for chrom in pybw_header}
	logger.debug("Chromosome lengths from input bigwig: {0}".format(chrom_info))

	#Decide regions 
	logger.info("- Getting output regions ready")
	if args.regions:
		regions = RegionList().from_bed(args.regions)

		#Check whether regions are available in input bigwig
		not_in_bigwig = list(set(regions.get_chroms()) - set(chrom_info.keys()))
		if len(not_in_bigwig) > 0:
			logger.warning("Contigs {0} were found in input --regions, but were not found in input --signal. These regions cannot be scored and will therefore be excluded from output.".format(not_in_bigwig))
			regions = regions.remove_chroms(not_in_bigwig)
		
		regions.apply_method(OneRegion.extend_reg, args.extend)
		regions.merge()
		regions.apply_method(OneRegion.check_boundary, chrom_info, "cut")

	else:
		regions = RegionList().from_chrom_lengths(chrom_info)

	#Set flank to enable scoring in ends of regions
	if args.score == "sum":
		args.region_flank = int(args.window/2.0)
	elif args.score == "footprint" or args.score == "FOS":
		args.region_flank = int(args.flank_max)
	else:
		args.region_flank = 0

	#Go through each region
	for i, region in enumerate(regions):
		region.extend_reg(args.region_flank)
		region = region.check_boundary(chrom_info, "cut")
		region.extend_reg(-args.region_flank)

	#Information for output bigwig
	reference_chroms = sorted(list(chrom_info.keys()))
	header = [(chrom, chrom_info[chrom]) for chrom in reference_chroms]
	regions.loc_sort(reference_chroms)

	#---------------------------------------------------------------------------------------#
	#------------------------ Calculating footprints and writing out -----------------------#
	#---------------------------------------------------------------------------------------#

	logger.info("Calculating footprints in regions...")
	regions_chunks = regions.chunks(args.split)

	#Setup pools
	args.cores = check_cores(args.cores, logger)
	writer_cores = 1	
	worker_cores = max(1, args.cores - writer_cores)
	logger.debug("Worker cores: {0}".format(worker_cores))
	logger.debug("Writer cores: {0}".format(writer_cores))

	worker_pool = mp.Pool(processes=worker_cores)
	writer_pool = mp.Pool(processes=writer_cores)
	manager = mp.Manager()

	#Start bigwig file writers
	q = manager.Queue()
	writer_pool.apply_async(bigwig_writer, args=(q, {"scores":args.output}, header, regions, args))
	writer_pool.close() #no more jobs applied to writer_pool
	writer_qs = {"scores": q}

	args.writer_qs = writer_qs

	#Start calculating scores
	pool = mp.Pool(processes=args.cores)
	task_list = [pool.apply_async(calculate_scores, args=[chunk, args]) for chunk in regions_chunks]
	no_tasks = len(task_list)
	pool.close()
	monitor_progress(task_list, logger)
	results = [task.get() for task in task_list]

	#Stop all queues for writing
	logger.debug("Stop all queues by inserting None")
	for q in writer_qs.values():
		q.put((None, None, None))

	#Done computing
	writer_pool.join() 
	worker_pool.terminate()
	worker_pool.join()
	
	logger.stop_logger_queue()

	#Finished scoring
	logger.end()


#--------------------------------------------------------------------------------------------------------#
if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser = add_scorebigwig_arguments(parser)
	args = parser.parse_args()

	if len(sys.argv[1:]) == 0:
		parser.print_help()
		sys.exit()
		
	run_footprinting(args)
