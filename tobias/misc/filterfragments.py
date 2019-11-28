#!/usr/bin/env python

"""
Small utility to filter .bam-file fragments based on overlap with .bed-regions
One use-case is to filter out fragments arising from gene-containing plasmids which contaminate ATAC-seq with reads mapping to exons

@author: Mette Bentsen
@contact: mette.bentsen (at) mpi-bn.mpg.de
@license: MIT
"""

import os
import sys
import argparse
import pysam

from tobias.utils.regions import *
from tobias.utils.utilities import *

#----------------------------------------------------------------------------------------#
def add_filterfragments_arguments(parser):

	parser.formatter_class = lambda prog: argparse.RawDescriptionHelpFormatter(prog, max_help_position=40, width=90)
	description = "FilterFragments can filter out fragments from a .bam-file of paired-end reads based on the overlap with regions in the given .bed-file."
	parser.description = format_help_description("FilterFragments", description)

	parser._action_groups.pop()	#pop -h

	IO = parser.add_argument_group('Input / output arguments')
	IO.add_argument('--bam', metavar="", help=".bam-file to filter")
	IO.add_argument('--regions', metavar="", help=".bed-file containing regions to filter fragments from")
	IO.add_argument('--mode', help="Mode 1: Remove fragment if both reads overlap .bed-regions. Mode 2: Remove whole fragment if one read is overlapping .bed-regions (default: 1)", choices=[1,2], default=1)
	IO.add_argument('--output', metavar="", help="Path to the filtered .bam-file (default: <prefix of --bam>_filtered.bam)")

	IO = add_logger_args(IO)

	return(parser)

#----------------------------------------------------------------------------------------#
def run_filterfragments(args):

	#Get output filename
	args.output = os.path.splitext(os.path.basename(args.bam))[0] + "_filtered.bam" if args.output == None else args.output

	#Start logger
	logger = TobiasLogger("FilterFragments", args.verbosity)
	logger.begin()

	parser = add_filterfragments_arguments(argparse.ArgumentParser())
	logger.arguments_overview(parser, args)

	check_required(args, ["bam", "regions"])


	################### Find fragments to filter ####################

	#Read regions
	regions = RegionList().from_bed(args.regions)
	logger.info("Read {0} regions from --regions".format(len(regions)))

	#Open bam and overlap with regions
	logger.info("Fetching reads from regions")
	bam_obj = pysam.AlignmentFile(args.bam, "rb")

	all_reads = {}	#dict for counting fragments within regions
	for region in regions:
		logger.debug(region)
		reads = bam_obj.fetch(region.chrom, region.start, region.end)
		for read in reads:
			all_reads[read.query_name] = all_reads.get(read.query_name, []) + [read]	#query_name is the unique fragment id

	logger.info("Found a total of {0} fragments overlapping regions".format(len(all_reads)))

	#Filter fragments based on mode
	if args.mode == 1:
		excluded_reads = set([name for name in all_reads if len(all_reads[name]) > 1])
	elif args.mode == 2:
		excluded_reads = list(all_reads.keys())
	else: 
		pass

	logger.info("Found {0} fragments to filter (mode {1})".format(len(excluded_reads), args.mode))


	################### Write filtered bam ###################

	logger.info("Writing filtered file")
	obam = pysam.AlignmentFile(args.output, "wb", template=bam_obj, threads=60)
	for b in bam_obj.fetch(until_eof=True):
		if b.query_name not in excluded_reads:
			obam.write(b)

	bam_obj.close()
	obam.close()

	#Index newly created bam
	pysam.index(args.output, args.output + ".bai")

	logger.end()
