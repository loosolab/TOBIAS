#!/usr/bin/env python

"""
CreateNetwork: Creates a TF-TF gene regulation network from annotated transcription factor binding sites

@author: Mette Bentsen
@contact: mette.bentsen (at) mpi-bn.mpg.de
@license: MIT

"""

import os
import sys
import argparse
import pyBigWig
import numpy as np
import glob

#import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd

#Utils from TOBIAS
from tobias.utils.utilities import * 
from tobias.utils.logger import *

#--------------------------------------------------------------------------------#

def add_network_arguments(parser):

	parser.formatter_class = lambda prog: argparse.RawDescriptionHelpFormatter(prog, max_help_position=40, width=90)
	description = "Creates a TF-TF gene regulation network from annotated transcription factor binding sites" 
	parser.description = format_help_description("CreateNetwork", description)

	parser._action_groups.pop()	#pop -h
	
	#Required arguments
	required = parser.add_argument_group('Required arguments')
	required.add_argument('--TFBS', metavar="", help="File(s) containing TFBS (with annotation) to create network from", nargs="*")
	required.add_argument('--origin', metavar="", help="File containing mapping of TF <-> origin gene")

	#Optional arguments
	optional = parser.add_argument_group('Optional arguments')
	optional.add_argument('--start', metavar="", help="Name of node to start in (default: all nodes)")
	optional.add_argument('--max-len', metavar="", help="Maximum number of nodes in paths through graph (default: 4)", type=int, default=4)
	#optional.add_argument('--unique', action='store_true', help="Only include edges once (default: edges can occur multiple times in case of multiple binding sites)")
	
	runargs = parser.add_argument_group("Run arguments")
	runargs.add_argument('--outdir', metavar="", help="Path to output directory (default: tobias_network)", default="tobias_network") 
	runargs = add_logger_args(runargs)

	return(parser)

#--------------------------------------------------------------------------------#
def dfs(adjacency, path, all_paths = [], options={"max_length":3}):

	last_node = path[-1] 
	if last_node in adjacency:
		target_nodes = adjacency[last_node].get("targets", [])
		if len(target_nodes) > 0:

			#Go through all targets and get paths
			for target_node in target_nodes:
				if target_node in path:    #cyclic; add node and close path
					new_path = path + [target_node]
					all_paths += [new_path]

				else:					  #not cyclic; add node and search for new paths
					new_path = path + [target_node]
					
					if len(new_path) == options["max_length"]:
						all_paths += [new_path]			#Save current path
					else:
						all_paths = dfs(adjacency, new_path, all_paths, options)		#Find targets of last node in path
		else:
			all_paths += [path] 
	
	else: 	#node is not in adjacency/has no targets; save current path
		all_paths += [path]

	return all_paths


#--------------------------------------------------------------------------------#
def run_network(args):

	make_directory(args.outdir)
	check_required(args, ["TFBS", "origin"])

	logger = TobiasLogger("CreateNetwork", args.verbosity)
	logger.begin()

	#-------------------------- Origin file translating motif name -> gene origin -----------------------------------#
	#translation file, where one motif can constitute more than one gene (jun::fos) 
	#and more genes can encode transcription factors with same motifs (close family members with same target sequence)
	origin_table = pd.read_csv(args.origin, sep="\t", header=None)
	origin_table.columns = ["Origin_" + str(element) for element in origin_table.columns]

	#------------------------ Transcription factor binding sites with mapping to target genes -----------------------#

	logger.info("Reading all input binding sites")
	
	#todo: read in parallel
	dataframes = []
	for fil in args.TFBS:
		logger.debug("- {0}".format(fil))

		df = pd.read_csv(fil, sep="\t", header=None)
		dataframes.append(df)

	logger.debug("Joining sites from all files")
	sites_table = pd.concat(dataframes, sort=False)
	sites_table.columns = ["Sites_" + str(element) for element in sites_table.columns]
	logger.info("- Total of {0} sites found\n".format(sites_table.shape[0]))

	#---------------------------------------- Match target columns to origin ----------------------------------------#

	logger.info("Matching target genes back to TFs using --origin")

	origin_table_str_columns = list(origin_table.dtypes.index[origin_table.dtypes == "object"])
	sites_table_str_columns = list(sites_table.dtypes.index[sites_table.dtypes == "object"])

	origin_table = origin_table.apply(lambda x: x.astype(str).str.upper())
	sites_table = sites_table.apply(lambda x: x.astype(str).str.upper())
	
	#Establishing matched columns
	logger.debug("Establishing which columns should be used for mapping target -> source")
	
	matching = []
	for sites_column in sites_table_str_columns:
		sites_column_content = set(sites_table[sites_column])
		for origin_column in origin_table_str_columns:
			origin_column_content = set(origin_table[origin_column])

			#Overlap
			overlap = len(origin_column_content & sites_column_content)
			matching.append((sites_column, origin_column, overlap))

	sorted_matching = sorted(matching, key = lambda tup: -tup[-1])
	logger.debug("Match tuples: {0}".format(sorted_matching))

	#Columns for matching
	source_col_tfbs = sites_table.columns[3]	#Source is the 4th column of the bedfile 
	source_col_origin = [match[1] for match in sorted_matching if match[0] == source_col_tfbs][0]  #source id in the origin table
	target_col_tfbs = [match[0] for match in sorted_matching if match[0] != source_col_tfbs][0]	 #target id in from sites 
	target_col_origin = [match[1] for match in sorted_matching if match[0] != source_col_tfbs][0]	 #target id in the origin table
		
	#Intersect of sources and targets
	source_ids_tfbs = set(sites_table[source_col_tfbs])
	logger.debug("Source ids from TFBS: {0}".format(", ".join(list(source_ids_tfbs)[:4]) + " (...)"))

	source_ids_origin = set(origin_table[source_col_origin])
	logger.debug("Matched source ids from origin table: {0}".format(", ".join(list(source_ids_origin)[:4]) + " (...)"))

	target_ids_tfbs = set(sites_table[target_col_tfbs])
	logger.debug("Target ids from TFBS: {0}".format(", ".join(list(target_ids_tfbs)[:4]) + " (...)"))

	target_ids_origin = set(origin_table[target_col_origin])
	logger.debug("Matched target ids from origin table: {0}".format(", ".join(list(target_ids_origin)[:4]) + " (...)"))

	common_ids = source_ids_tfbs & source_ids_origin
	if len(common_ids) != len(source_ids_tfbs):
		missing_ids = source_ids_tfbs - common_ids
		logger.warning("Source ids from TFBS could not be found in the --origin table: {0}".format(missing_ids))

	#Add target motif id to sites
	sites_table_convert = sites_table.merge(origin_table[[source_col_origin, target_col_origin]], left_on=target_col_tfbs, right_on=target_col_origin, how="inner")
	logger.info("Subset sites to {0} with annotation genes encoding TFs".format(sites_table_convert.shape[0]))
	
	#Subset to unique edges if chosen
	#if args.unique == True:
	#	sites_table_convert.drop_duplicates(subset=[source_col_tfbs, source_col_origin], inplace=True)
	#	logger.info("Flag --unique is on: Sites were further subset to {0} unique edges".format(sites_table_convert.shape[0]))

	#------------------------------------ Write out edges / adjacency ----------------------------------------#

	##### Write out edges #####
	edges_f = os.path.join(args.outdir, "edges.txt")
	logger.info("Writing edges list to: {0}".format(edges_f))
	sites_table_convert.to_csv(edges_f, sep="\t", index=False)

	###### Create adjacency list ####
	logger.info("Creating adjacency matrix")
	adjacency = {source: {"targets":[]} for source in common_ids}
	for index, row in sites_table_convert.iterrows():
		source, target = row[source_col_tfbs], row[source_col_origin]
		if target not in adjacency[source]["targets"]:
			adjacency[source]["targets"].append(target)

	adjacency_f = os.path.join(args.outdir, "adjacency.txt")
	logger.info("- Writing matrix to: {0}".format(adjacency_f))
	with open(adjacency_f, "w") as f:
		f.write("Source\tTargets\n")
		for source in sorted(adjacency):
			f.write("{0}\t{1}\n".format(source, ", ".join(adjacency[source]["targets"])))

	#-------------------------------------- Find paths through graph ---------------------------------------#

	#Create possible paths through graph
	logger.info("")
	logger.info("Create possible paths through graph")

	#Starting node can be used to subset path-finding to specific nodes; speeds up computation
	if args.start != None:
		start_nodes = [one_id for one_id in common_ids if args.start.upper() in one_id]
		logger.info("Starting nodes are: {0}".format(start_nodes))
	else:
		start_nodes = common_ids
		logger.info("Finding paths starting from all nodes. This behavior can be changed using --start.")
	
	#Find paths starting at source nodes
	for source_id in start_nodes:
		logger.info("- Finding paths starting from {0}".format(source_id))

		#Recursive function to find paths; initiate with no paths found
		paths = dfs(adjacency = adjacency, path = [source_id], all_paths=[], options={"max_length": args.max_len})

		paths_f = os.path.join(args.outdir, "{0}_paths.txt".format(source_id))
		logger.debug("-- Writing paths to: " + paths_f)
		paths_out = open(paths_f, "w")
		paths_out.write("Regulatory_path\tn_nodes\n")

		path_edges = []		#Collect edges while writing paths
		for path in paths:
			if len(path) > 1: #only write paths longer than 1 node

				#String formatting of path
				str_paths = ""
				for i, node in enumerate(path[:-1]):
					str_paths += "{0} --> ".format(node)
				str_paths += path[-1]

				n_nodes = len(path)
				paths_out.write("{0}\t{1}\n".format(str_paths, n_nodes))

				#Make pairwise edges across path
				path_edges += ["\t".join([path[i], path[j], str(j)]) for i,j in zip(range(0, len(path)-1), range(1,len(path)))]

		paths_out.close()
		
		#Write out the edges included in paths
		source_paths_f = os.path.join(args.outdir, "{0}_path_edges.txt".format(source_id))
		with open(source_paths_f, "w") as f:
			f.write("Source\tTarget\tLevel\n")
			for path in sorted(set(path_edges), key=lambda tup: (tup[-1], tup[0])):
				f.write(path + "\n")
	
	#Finish CreateNetwork
	logger.end()
