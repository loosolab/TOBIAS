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

from tobias.parsers import add_filterfragments_arguments
from tobias.utils.regions import RegionList
from tobias.utils.utilities import *
from tobias.utils.ngs import index_bam

from tobias import __version__ as TOBIAS_VERSION

#----------------------------------------------------------------------------------------#
def run_filterfragments(args):

	#Get output filename
	args.output = os.path.splitext(os.path.basename(args.bam))[0] + "_filtered.bam" if args.output == None else args.output

	#Start logger
	logger = TobiasLogger("FilterFragments", args.verbosity)
	logger.begin()

	#Get parser arguments
	parser = add_filterfragments_arguments(argparse.ArgumentParser())
	logger.arguments_overview(parser, args)
	check_required(args, ["bam", "regions"])

	################### Find fragments to filter ####################

	#Read regions
	regions = RegionList().from_bed(args.regions)
	logger.stats("Read {0} regions from --regions".format(len(regions)))

	#Open bam 
	logger.info("Fetching reads from regions")
	bam_obj = pysam.AlignmentFile(args.bam)
	bam_obj = index_bam(bam_obj, logger)

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
		logger.spam("Fetching reads from region: {0}".format(region))
		reads = bam_obj.fetch(region.chrom, region.start, region.end)
		for read in reads:
			all_frags[read.query_name] = all_frags.get(read.query_name, []) + [read]	#query_name is the unique fragment id
	logger.stats("Found a total of {0} fragments overlapping regions".format(len(all_frags)))

	#Filter fragments based on mode
	if args.mode == 1:
		excluded_fragments = set([name for name in all_frags if len(all_frags[name]) > 1])	#only exclude if both reads are found within regions
	elif args.mode == 2:
		excluded_fragments = list(all_frags.keys())
	else: 
		pass

	logger.stats("Selected {0} fragments to remove (mode {1})".format(len(excluded_fragments), args.mode))

	################### Write filtered bam ###################

	logger.debug("Setting up header for filtered file")
	header_dict = bam_obj.header.to_dict()

	#Add FilterFragments to PG
	if not "PG" in header_dict:
		header_dict["PG"] = []

	pg_dict = {"ID":"FilterFragments", "PN":"TOBIAS", "CL": logger.cmd, "VN": TOBIAS_VERSION}
	header_dict["PG"].append(pg_dict)
	logger.debug("New header dict: {0}".format(header_dict))

	logger.info("Writing filtered file")
	obam = pysam.AlignmentFile(args.output, "wb", header=header_dict, threads=args.threads)
	logger.debug("Opened {0}".format(args.output))
	for read in bam_obj.fetch(until_eof=True, multiple_iterators=True):		#open again to reset iterator
		if read.query_name not in excluded_fragments:
			obam.write(read)

	bam_obj.close()
	obam.close()

	#Index newly created bam
	logger.debug("Indexing filtered bam ({0})".format(args.output))
	pysam.index(args.output, args.output + ".bai")

	logger.end()
