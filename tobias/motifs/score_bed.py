#!/usr/bin/env python

"""
Scorebed: Scores a bedfile with signal from bigwig file(s)

@author: Mette Bentsen
@contact: mette.bentsen (at) mpi-bn.mpg.de
@license: MIT

"""

#--------------------------------------------------------------------------------------------------------#
#----------------------------------------- Import libraries ---------------------------------------------# 
#--------------------------------------------------------------------------------------------------------#

import os
import sys
import argparse
import pyBigWig
import numpy as np
from datetime import datetime

#Bio-stuff
import pybedtools as pb

#Utils from TOBIAS
from tobias.utils.regions import *
from tobias.utils.utilities import * 
from tobias.utils.logger import *

#-------------------------------------------------------------------------------------------#
#-------------------------------- Command line arguments -----------------------------------#
#-------------------------------------------------------------------------------------------#	

def add_scorebed_arguments(parser):

	parser.formatter_class = lambda prog: argparse.RawDescriptionHelpFormatter(prog, max_help_position=40, width=90)
	description = "ScoreBed is a utility to score .bed-file regions with values from a .bigwig-file. The output is a .bed-file with the bigwig value(s) as extra column(s). Options --position and --math can be used to adjust scoring scheme."
	parser.description = format_help_description("ScoreBed", description)

	parser._action_groups.pop()	#pop -h
	
	#Required arguments
	required = parser.add_argument_group('Required arguments')
	required.add_argument('--bed', metavar="", help="Sites to score (.bed file)")
	required.add_argument('--bigwigs', metavar="", nargs="*",  help="Scores to assign to regions in .bed (.bw file(s))")
	
	#Optional arguments
	optional = parser.add_argument_group('Optional arguments')
	optional.add_argument('--output', metavar="", help="Path to output .bed-file (default: scored sites are written to stdout") 
	optional.add_argument('--subset', metavar="", help="Subset scoring to .bed regions and set all other sites to --null value (default: all sites in input file will be scored)")
	optional.add_argument('--null', metavar="", help="If --subset is given, which score/label to add to non-scored regions (default: 0)", default="0", type=float)
	optional.add_argument('--position', metavar="", help="Position in sites to score (start/mid/end/full) (default: full)", choices=["mid", "start", "end", "full"], default="full")
	optional.add_argument('--math', metavar="", help="If position == full, choose math to perform on signal (min/max/mean/sum) (default: mean)", choices=["min", "max", "mean", "sum"], default="mean")
	optional = add_logger_args(optional)
	#optional.add_argument('--buffer', metavar="", help="Lines to buffer before writing (default: 10000)", type=int, default=10000)

	return(parser)

#--------------------------------------------------------------------------------#

def get_score_func(args):

	if args.position == "start":
		func = lambda signal: signal[0]

	elif args.position == "mid":
		func = lambda signal: signal[int(len(signal)/2.0)]
	
	elif args.position == "end":
		func = lambda signal: signal[-1]

	elif args.position == "full":
		if args.math == "min":
			func = lambda signal: np.min(signal)
		elif args.math == "max":	
			func = lambda signal: np.max(signal)
		elif args.math == "mean":
			func = lambda signal: np.mean(signal)
		elif args.math == "sum":
			func = lambda signal: np.sum(signal)
	
	return(func)

#--------------------------------------------------------------------------------#
def run_scorebed(args):

	#Verbosity is 0 if output is written to stdout
	if args.output == None:
		args.verbosity = 0
	
	#Start logger
	logger = TobiasLogger("ScoreBed", args.verbosity)
	logger.begin()

	parser = add_scorebed_arguments(argparse.ArgumentParser())
	logger.arguments_overview(parser, args)
	logger.output_files([args.output])

	#Check input 
	check_required(args, ["bed", "bigwigs"])
	check_files([getattr(args, arg) for arg in ["bed", "bigwigs", "subset"]], "r")
	check_files([getattr(args, "output")], "w")

	#Setup output file/stdout
	if args.output != None:
		sys.stdout = open(args.output, 'w')
	no_bigwigs = len(args.bigwigs)



	#----------------------------------------------------------------------#
	#----------------------- Overlap bed if needed ------------------------#
	#----------------------------------------------------------------------#

	if args.subset != None:
		logger.info("Overlapping {0} to regions in {1}".format(args.bed, args.subset))
		
		#Make overlap
		pb_bed = pb.BedTool(args.bed)
		pb_region = pb.BedTool(args.subset)
		overlap = pb_bed.intersect(pb_region, c=True)
		args.bed = overlap.fn

	# Setup score function
	score_func = get_score_func(args)


	#----------------------------------------------------------------------#
	#---------------------- Open bw and run scoring -----------------------#
	#----------------------------------------------------------------------#

	pybw = {bigwig_f: pyBigWig.open(bigwig_f) for bigwig_f in args.bigwigs}		#open filehandles for all bigwigs

	logger.info("Starting scoring...")
	#buff = [""]*args.buffer
	#buff_idx = 0
	count = 0
	with open(args.bed) as f:
		for line in f:
			
			columns = line.strip().split("\t")
			outline = "\t".join(columns) if args.subset == None else "\t".join(columns[:-1])	 #don't save overlap count
			chrom, start, end = columns[0], int(columns[1]), int(columns[2])

			#Decide whether to get signal or not
			overlap = 1 if args.subset == None else int(columns[-1])
			if overlap != 0: 
				for bigwig_f in args.bigwigs:	#preserves order of bigwigs
					
					signal = pybw[bigwig_f].values(chrom, start, end, numpy=True)

					if len(signal) > 0:
						signal = np.nan_to_num(signal)
						score = round(score_func(signal), 5)
					else:
						score = args.null
										
					outline += "\t" + "{0:.5f}".format(score)
			else:
				outline += "\t" + "\t".join([args.null]*no_bigwigs) 
			
			#outline += "\n"
			#buff[buff_idx] = outline
			#buff_idx += 1

			print(outline)

			#if buff_idx == args.buffer: 	#eg. if idx is 100 (already 0-99 lines in buffer) and args.buffer is 100, this will be out of bounds for next line > write out
			#	fout.write("\n".join(buff) + "\n")

			#	#Init
			#	buff = [""]*args.buffer
			#	buff_idx = 0

	for bigwig_f in pybw:
		pybw[bigwig_f].close()
	
	logger.end()

#--------------------------------------------------------------------------------------------------------#
if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser = add_scorebed_arguments(parser)
	args = parser.parse_args()

	if len(sys.argv[1:]) == 0:
		parser.print_help()
		sys.exit()
		
	run_scorebed(args)
