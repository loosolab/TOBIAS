#!/usr/bin/env python

"""
ClusterTFBS: Performs TF clustering based on positions of TFBS across genome by calculating overlaps of sites

@author: Mette Bentsen
@contact: mette.bentsen (at) mpi-bn.mpg.de
@license: MIT
"""

import argparse
import tempfile
import sys
import os
import numpy as np
import itertools
from datetime import datetime
import copy
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
import matplotlib.pyplot as plt

from tobias.utils.regions import *
from tobias.utils.utilities import *
from tobias.utils.logger import *



def add_clustering_arguments(parser):

	parser.formatter_class = lambda prog: argparse.RawDescriptionHelpFormatter(prog, max_help_position=40, width=90)
	description = ""
	parser.description = format_help_description("ClusterTFBS", description)

	parser._action_groups.pop()	#pop -h
	
	args = parser.add_argument_group('Arguments')
	args.add_argument('-b', '--bedfiles', metavar="", nargs="*", help="Bedfile(s) containing TFBS sites with name in the 4th column")
	args.add_argument('-o', '--outdir', metavar="", help="Path of output folder (default: current directory)", default=".")
	args.add_argument('-p', '--prefix', metavar="", help="Prefix of output files (default: ClusterTFBS)", default="ClusterTFBS")
	args.add_argument('-c', '--cores', metavar="<int>", type=int, help="Number of cores to use for computation (default: 1)", default=1)
	args.add_argument('--threshold', metavar="<float>", help="Threshold for clustering (default: 0.5)", type=lambda f: restricted_float(f,0,1), default=0.5)
	args.add_argument('--method', metavar="", help="Method for hierarchical clustering (single/complete/average/weighted/centroid/median/ward) (default: complete)", choices=["single", "complete", "average", "weighted", "centroid", "median", "ward"], default="complete")
	args = add_logger_args(args)
	
	return(parser)

def overlap_sites(f):
	""" overlap bed sites from file """
	
	sites = RegionList().from_bed(f)
	overlap = sites.count_overlaps()
	
	return(overlap)

#############################################################################################

def run_clustering(args):

	#Create logger
	logger = TobiasLogger("ClusterTFBS", args.verbosity)
	logger.begin()

	parser = add_clustering_arguments(argparse.ArgumentParser())
	logger.arguments_overview(parser, args)

	#Output files
	#output_files = ""
	#logger.output_files(output_files)

	#Check input
	check_required(args, ["bedfiles"])	#Check input arguments
	make_directory(args.outdir)			#Check output directory; create if needed

	#----------------------------------------------------------#

	#Join all sites
	logger.info("Joining and sorting all sites")

	#tempprefix
	temp_prefix = next(tempfile._get_candidate_names())
	logger.debug("Temp prefix: {0}".format(temp_prefix))

	joined_f = os.path.join(args.outdir, "{0}.tmp".format(temp_prefix))
	cmd = "cat {0} | cut -f 1-4 | sort -k1,1 -k2,2n > {1}".format(" ".join(args.bedfiles), joined_f)
	os.system(cmd)

	logger.info("Splitting sites per chromosome")
	handles = {}
	splitted = []
	with open(joined_f) as f:
		for line in f:

			columns = line.split("\t")
			chrom = columns[0]

			if chrom not in handles:
				logger.debug("Opening handle for chrom {0}".format(chrom))
				split_f = os.path.join(args.outdir, "{0}_{1}.tmp".format(temp_prefix, chrom))
				splitted.append(split_f)
				handles[chrom] = open(split_f, "w")

			handles[chrom].write(line)	

	logger.debug("Closing handles")
	for chrom in handles:
		handles[chrom].close()

	#Overlap sites
	logger.info("Overlapping sites within subsets")
	results = run_parallel(overlap_sites, splitted, [], args.cores, logger)

	#Join results from individual overlaps
	logger.info("Joining results from subsets")
	global_overlap = merge_dicts(results)


	#----------------- Create distance matrix -------------------#

	clustering = RegionCluster(global_overlap)
	logger.stats("Found {0} unique IDS".format(clustering.n))
	
	clustering.cluster(threshold=args.threshold, method=args.method)
	
	#Write out distance mat	
	outfile = os.path.join(args.outdir, args.prefix + "_distance_matrix.txt")
	clustering.write_distance_mat(outfile)

	#Write out clusters
	outfile = open(os.path.join(args.outdir, args.prefix + "_clusters.txt"), "w")
	outfile.write("#cluster_name\tcluster_size\tcluster_members\n")

	i = 1
	for cluster in clustering.clusters:
		members = clustering.clusters[cluster]["member_names"]
		outfile.write("C_{0}\t{1}\t{2}\n".format(i, len(members), ", ".join(members)))

	outfile.close()

	#--------------------------- Plot dendrogram -----------------------------#

	fig, ax = plt.subplots(figsize = (3,max(10,len(IDS)/10)))	#5 rows per one column width
	dendro_dat = dendrogram(clustering.linkage_mat, labels=clustering.names, orientation="right", ax=ax, above_threshold_color="black", link_color_func=lambda k: clustering.node_color[k])

	figure_f = os.path.join(args.outdir, args.prefix + "_dendrogram.pdf")
	fig.savefig(figure_f, bbox_inches='tight')
	logger.info("Plotted dendrogram to: {0}".format(figure_f))

	#---------------------- Remove tempfiles and end -------------------------#

	for fil in [joined_f] + splitted:
		logger.debug("Removing {0}".format(fil))
		try:
			os.remove(fil)
		except:
			logger.warning("Could not remove tempfile {0}".format(fil))

	logger.end()

#------------------------------------------------------------------------------------------#

if __name__ == "__main__":

	parser = argparse.ArgumentParser()
	parser = add_clustering_arguments(parser)
	args = parser.parse_args()

	if len(sys.argv[1:]) == 0:
		parser.print_help()
		sys.exit()

	run_clustering(args)

