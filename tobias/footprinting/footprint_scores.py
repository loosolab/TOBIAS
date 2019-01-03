#!/usr/bin/env python

"""
FootprintScores: Calculate footprint tracks from cutsite bigwig

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
from datetime import datetime
from scipy import stats

#Internal functions and classes
from tobias.utils.utilities import *
from tobias.utils.regions import *
from tobias.utils.sequences import *
from tobias.utils.signals import *
from tobias.utils.ngs import *

#-----------------------------------------------------------------#
def add_footprint_arguments(parser):

	parser.formatter_class = lambda prog: argparse.RawDescriptionHelpFormatter(prog, max_help_position=40, width=90)
	description = "FootprintScores calculates footprint scores from ATAC-seq cutsites (.bigwig format) calculated using the ATACorrect tool.\n\n"
	description += "Usage: FootprintScores --signal <cutsites.bw> --output <output.bw>\n\n"
	description += "Output:\n- <footprint_scores.bigwig>"
	parser.description = format_help_description("FootprintScores", description)
	
	parser._action_groups.pop()	#pop -h

	#Required arguments
	required = parser.add_argument_group('Required arguments')
	required.add_argument('-s', '--signal', metavar="<bigwig>", help="A .bw file of ATAC-seq cutsite signal")
	required.add_argument('-o', '--output', metavar="<bigwig>", help="Full path to output bigwig")			
	required.add_argument('-r', '--regions', metavar="<bed>", help="Genomic regions to run footprinting in")

	optargs = parser.add_argument_group('Optional arguments')
	optargs.add_argument('--extend', metavar="<int>", type=int, help="Extend input regions with bp (default: 100)", default=100)
	optargs.add_argument('--score', metavar="<score>", choices=["tobias", "FOS", "sum"], help="Type of scoring to perform on cutsites (tobias/FOS/sum) (default: tobias)", default="tobias")
	
	optargs.add_argument('--fp_min', metavar="<int>", type=int, help="Minimum footprint width (default: 15)", default=15)
	optargs.add_argument('--fp_max', metavar="<int>", type=int, help="Maximum footprint width (default: 50)", default=50)
	optargs.add_argument('--flank_min', metavar="<int>", type=int, help="Minimum range of flanking regions (default: 10)", default=10)
	optargs.add_argument('--flank_max', metavar="<int>", type=int, help="Maximum range of flanking regions (default: 30)", default=30)
	optargs.add_argument('--sum_window', metavar="<int>", type=int, help="For --score == sum; Window for calculation of sum (default: 200)", default=200)
	optargs.add_argument('--smooth', metavar="<int>", type=int, help="Smooth output signal by mean in <bp> windows (default: 5)", default=5)
	optargs.add_argument('--min_limit', metavar="<float>", type=float, help="Limit input bigwig score range (default: no lower limit)") 			#default none
	optargs.add_argument('--max_limit', metavar="<float>", type=float, help="Limit input bigwig score range (default: no upper limit)") 	#default none

	runargs = parser.add_argument_group('Run arguments')
	runargs.add_argument('--cores', metavar="<int>", type=int, help="Number of cores to use for computation (default: 1)", default=1)
	runargs.add_argument('--split', metavar="<int>", type=int, help="Split of multiprocessing jobs (default: 100)", default=100)
	runargs.add_argument('--verbosity', metavar="<int>", help="Level of output logging (1 (sparse) / 2 (normal) / 3 (debug)) (default: 2)", choices=[1,2,3], default=2, type=int)
	runargs.add_argument('--log', metavar="<file>", help="Full path of logfile (default: log is printed to stdout)")
	runargs.add_argument('--debug', help=argparse.SUPPRESS, action='store_true')

	return(parser)

#------------------------------------------------------------------#
def calculate_scores(regions, args):

	pybw_signal = pyBigWig.open(args.signal) 	#cutsites signal
	pybw_header = pybw_signal.chroms()			
	chrom_lengths = {chrom:int(pybw_header[chrom]) for chrom in pybw_header}

	output_dict = {}

	#Set flank to enable scoring in ends of regions
	if args.score == "sum":
		flank = int(args.sum_window/2.0)
	elif args.score == "tobias" or args.score == "FOS":
		flank = args.flank_max
	else:
		flank = args.flank_max

	#Go through each region
	for i, region in enumerate(regions):

		#Extend region with necessary flank
		region.extend_reg(flank)
		region.check_boundary(chrom_lengths, "cut")	
		reg_key = (region.chrom, region.start+flank, region.end-flank)	#output region

		#Get bigwig signal in region
		signal = region.get_signal(pybw_signal)
		signal = np.nan_to_num(signal).astype("float64")

		if args.min_limit != None:
			signal[signal < args.min_limit] = args.min_limit
		if args.max_limit != None:
			signal[signal > args.max_limit] = args.max_limit

		#Calculate scores
		if args.score == "sum":
			signal = np.abs(signal)
			scores = fast_rolling_math(signal, args.sum_window, "sum")

		elif args.score == "FOS":
			scores = calc_FOS(signal, args.fp_min, args.fp_max, args.flank_min, args.flank_max)

		elif args.score == "tobias":
			scores = tobias_footprint_array(signal, args.fp_min, args.fp_max, args.flank_min, args.flank_max)		#numpy array
		
		else:
			sys.exit("{0} not found".format(args.score))
		
		#Smooth signal with args.smooth bp
		if args.smooth > 1:
			scores = fast_rolling_math(scores, args.smooth, "mean")

		#Remove ends to prevent overlap with other regions
		if flank > 0:
			scores = scores[flank:-flank]
		output_dict[reg_key] = scores

	return(output_dict)

#------------------------------------------------------------------#

def run_footprinting(args):
	
	begin_time = datetime.now()
	
	check_required(args, ["signal", "output", "regions"])
	check_files([args.signal, args.regions], "r")
	check_files([args.output], "w")

	#---------------------------------------------------------------------------------------#
	# Create logger and write info to log
	#---------------------------------------------------------------------------------------#

	logger = create_logger(args.verbosity, args.log)
	logger.comment("#TOBIAS FootprintScores (run started {0})\n".format(begin_time))
	logger.comment("#Command line call: {0}\n".format(" ".join(sys.argv)))

	parser = add_footprint_arguments(argparse.ArgumentParser())
	logger.comment(arguments_overview(parser, args))

	logger.comment("# ----- Output files -----")
	logger.comment("# - {0}".format(args.output))
	logger.comment("\n\n")


	#---------------------------------------------------------------------------------------#
	#----------------------- I/O - get regions/bigwig ready --------------------------------#
	#---------------------------------------------------------------------------------------#

	logger.critical("Processing input files")
	logger.info("Opening cutsite bigwig")
	pybw_signal = pyBigWig.open(args.signal)
	pybw_header = pybw_signal.chroms()
	chrom_info = {chrom:int(pybw_header[chrom]) for chrom in pybw_header}

	logger.info("Getting output regions ready")

	#Decide regions 
	if args.regions:
		regions = RegionList().from_bed(args.regions)
		regions.apply_method(OneRegion.extend_reg, args.extend)
		regions.merge()
		regions.apply_method(OneRegion.check_boundary, chrom_info)
	else:
		regions = RegionList().from_list([OneRegion([chrom, 0, chrom_info[chrom]]) for chrom in chrom_info])	

	#Getting bigwig files ready
	logger.info("Opening output bigwig")
	reference_chroms = sorted(list(chrom_info.keys()))
	header = [(chrom, chrom_info[chrom]) for chrom in reference_chroms]
	pybw_footprints = pyBigWig.open(args.output, "w")
	pybw_footprints.addHeader(header)

	regions.loc_sort(reference_chroms)


	#---------------------------------------------------------------------------------------#
	#------------------------ Calculating footprints and writing out -----------------------#
	#---------------------------------------------------------------------------------------#

	logger.critical("Calculating footprints in regions...")

	if args.debug == True:
		fp = calculate_scores(regions[:100], args)

	regions_chunks = regions.chunks(args.split)

	#Start correction
	pool = mp.Pool(processes=args.cores)
	task_list = [pool.apply_async(calculate_scores, args=[chunk, args]) for chunk in regions_chunks]
	no_tasks = len(task_list)
	pool.close()

	chrom_order = {reference_chroms[i]:i for i in range(len(reference_chroms))}	

	#Process results as they come in
	write_idx = 0		#index in task_list
	prev_progress = (-1, -1)	#done/write
	while write_idx < len(task_list):
	
		if task_list[write_idx].ready() == True:	#Write result to file

			footprints = task_list[write_idx].get()

			#Write tracks to bigwig file
			for region in sorted(footprints.keys(), key=lambda tup: (chrom_order[tup[0]], tup[1], tup[2])): #Ensure that positions are written to bigwig in correct order
				logger.debug("Processing {0}".format(region))
				signal = footprints[region]

				chrom, reg_start, reg_end = region
				positions = np.arange(reg_start, reg_end)	#genomic 0-based positions

				zeros = np.logical_or(np.isclose(signal, 0), np.isnan(signal))
				signal[zeros] = 0						#adjust for weird numpy floating point
				
				included = signal.nonzero()[0]
				pos = positions[included]
				val = signal[included]

				if len(pos) > 0:
					try:
						pybw_footprints.addEntries(chrom, pos, values=val, span=1)
					except:
						logger.critical("Error writing region {0}.".format(region))
						sys.exit()
				
				#Clear memory
				footprints[region] = None

			#Wait for next region to print
			write_idx += 1 #This is also the same as the number of prints done 

		tasks_done = sum([task.ready() for task in task_list])
		if tasks_done != prev_progress[0] or write_idx != prev_progress[1]:
			logger.info("Calculation progress: {0:.0f}% | Writing progress: {1:.0f}%".format(tasks_done/no_tasks*100, write_idx/no_tasks*100))
		
		prev_progress = (tasks_done, write_idx)

	#Done computing
	pool.join()

	#Done writing to files
	logger.critical("Closing bigwig file (this can take a while)")
	pybw_footprints.close()

	end_time = datetime.now()
	logger.critical("Finished FootprintScores run (total time of {0})".format(end_time-begin_time))


#--------------------------------------------------------------------------------------------------------#

if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser = add_footprint_arguments(parser)
	args = parser.parse_args()

	if len(sys.argv[1:]) == 0:
		parser.print_help()
		sys.exit()
		
	run_footprinting(args)
