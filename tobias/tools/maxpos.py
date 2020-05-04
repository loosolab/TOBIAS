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
from tobias.utils.regions import OneRegion, RegionList
from tobias.utils.utilities import * 
from tobias.utils.logger import *

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
			region = OneRegion([chrom, start, end])
			signal = region.get_signal(pybw, numpy_bool=True)
			positions = minmax_func(signal)
	
			#output is list of index positions -> convert to bed
			rel_max_start, rel_max_end = min(positions), max(positions)    #position relative to peak start
			outline = "\t".join([chrom, str(start+rel_max_start), str(start+rel_max_end+1), ".", str(np.max(signal)), "."])	
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
