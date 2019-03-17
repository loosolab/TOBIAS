#!/usr/bin/env python

"""
Maxpos: Find the position of the maximum signal for each input region.  

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

#Utils from TOBIAS
from tobias.utils.regions import *
from tobias.utils.utilities import * 
from tobias.utils.logger import *


#-------------------------------------------------------------------------------------------#
#-------------------------------- Command line arguments -----------------------------------#
#-------------------------------------------------------------------------------------------#	

def add_maxpos_arguments(parser):

	parser.formatter_class = lambda prog: argparse.RawDescriptionHelpFormatter(prog, max_help_position=40, width=90)
	description = "MaxPos identifies the position of maximum signal (from bigwig) within a given set of .bed-regions. Used to identify peak of signals such as accessibility or footprint scores.\n\n"
	description += "Usage:\nTOBIAS MaxPos --bed <regions.bed> --bigwig <signal.bw>\n\n"
	description += "The output is a 4-column .bed-file (default output is stdout, but you can use --output to set a specific output file).\n"
	parser.description = format_help_description("MaxPos", description)

	parser._action_groups.pop()	#pop -h
	
	#Required arguments
	required = parser.add_argument_group('Required arguments')
	required.add_argument('--bed', metavar="", help="Regions to search for max position within")
	required.add_argument('--bigwig', metavar="", help="Scores used to identify maximum value")

	#Optional arguments
	optional = parser.add_argument_group('Optional arguments')
	optional.add_argument('--output', metavar="", help="Path to output .bed-file (default: scored sites are written to stdout)") 
	optional.add_argument('--invert', help="Find minimum position instead of maximum position", action='store_true', default=False)

	return(parser)

#--------------------------------------------------------------------------------#
def get_minmax_func(args):

	if args.invert == False:
		func = lambda signal: [i for i in range(len(signal)) if signal[i] == np.max(signal)]
	else:
		func = lambda signal: [i for i in range(len(signal)) if signal[i] == np.min(signal)]

	return(func)


#--------------------------------------------------------------------------------#
def run_maxpos(args):

	check_required(args, ["bed", "bigwig"])
	check_files([args.bed, args.bigwig], "r")
	check_files([args.output], "w")


	#Setup output file/stdout
	if args.output != None:
		sys.stdout = open(args.output, 'w')

	# Setup score function
	minmax_func = get_minmax_func(args)

	#----------------------------------------------------------------------#
	#---------------------- Open bw and run scoring -----------------------#
	#----------------------------------------------------------------------#

	pybw = pyBigWig.open(args.bigwig)
	#todo: test length of bigwig to scope of bed

	#Get maximum position(s) for each region
	with open(args.bed) as f:
		for line in f:

			columns = line.strip().split("\t")
			chrom, start, end = columns[0], int(columns[1]), int(columns[2])

			#Get signal from bigwig
			signal = pybw.values(chrom, start, end, numpy=True)
			signal = np.nan_to_num(signal)
			positions = minmax_func(signal)
	
			#output is list of index positions -> convert to bed
			rel_max_start, rel_max_end = min(positions), max(positions)    #position relative to peak start
			outline = "\t".join([chrom, str(start+rel_max_start-1), str(start+rel_max_end), ".", str(np.max(signal)), "."])	
			print(outline)

	pybw.close() 

#--------------------------------------------------------------------------------------------------------#
if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser = add_maxpos_arguments(parser)
	args = parser.parse_args()

	if len(sys.argv[1:]) == 0:
		parser.print_help()
		sys.exit()
		
	run_maxpos(args)
