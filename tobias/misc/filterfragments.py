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
	IO.add_argument('--threads', metavar="", help="Number of threads used for decompressing/compressing bam (default: 10)", default=10)

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

	#Open bam 
	logger.info("Fetching reads from regions")
	bam_obj = pysam.AlignmentFile(args.bam, "rb")

	#Check match between region chroms and bam obj:
	bam_chroms = bam_obj.references
	region_chroms = regions.get_chroms()
	if not set(region_chroms).issubset(set(bam_chroms)):
		chroms_not_found = set(region_chroms) - set(bam_chroms)
		logger.error("Contigs given in --regions do not align with contigs found in --bam! The bamfile contains: {0}. The following contigs from --regions are not found in --bam: {1}".format(bam_chroms, list(chroms_not_found)))
		sys.exit(1)

	#Get reads per fragment from bam within regions
	all_frags = {}	#dict for counting fragments within regions
	for region in regions:
		logger.debug(region)
		reads = bam_obj.fetch(region.chrom, region.start, region.end)
		for read in reads:
			all_frags[read.query_name] = all_frags.get(read.query_name, []) + [read]	#query_name is the unique fragment id
	logger.info("Found a total of {0} fragments overlapping regions".format(len(all_frags)))

	#Filter fragments based on mode
	if args.mode == 1:
		excluded_fragments = set([name for name in all_frags if len(all_frags[name]) > 1])	#only exclude if both reads are found within regions
	elif args.mode == 2:
		excluded_fragments = list(all_frags.keys())
	else: 
		pass

	logger.info("Found {0} fragments to filter (mode {1})".format(len(excluded_fragments), args.mode))


	################### Write filtered bam ###################

	logger.info("Writing filtered file")
	obam = pysam.AlignmentFile(args.output, "wb", template=bam_obj, threads=args.threads)
	for read in bam_obj.fetch(until_eof=True, multiple_iterators=True):		#open again to reset iterator
		if read.query_name not in excluded_fragments:
			obam.write(read)

	bam_obj.close()
	obam.close()

	#Index newly created bam
	pysam.index(args.output, args.output + ".bai")

	logger.end()
