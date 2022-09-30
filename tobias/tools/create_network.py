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
from tobias.parsers import add_network_arguments
from tobias.utils.utilities import * 
from tobias.utils.logger import TobiasLogger

#--------------------------------------------------------------------------------#
def dfs(adjacency, path, all_paths = [], options={"max_length":3}):
	""" Recursive function collecting all_paths through graph """

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
	check_required(args, ["TFBS", "origin"])	#check if anything is given for parameters
	check_files([args.TFBS, args.origin])

	logger = TobiasLogger("CreateNetwork", args.verbosity)
	logger.begin()

	#-------------------------- Origin file translating motif name -> gene origin -----------------------------------#
	#translation file, where one motif can constitute more than one gene (jun::fos) 
	#and more genes can encode transcription factors with same motifs (close family members with same target sequence)
	origin_table = pd.read_csv(args.origin, sep="\t", header=None)
	origin_table.columns = ["Origin_" + str(element) for element in origin_table.columns]
	origin_table.fillna("", inplace=True) #replace NaN with empty string

	#------------------------ Transcription factor binding sites with mapping to target genes -----------------------#

	logger.info("Reading all input binding sites")
	
	#todo: read in parallel
	dataframes = []
	for fil in args.TFBS:
		logger.debug("- {0}".format(fil))

		try:
			df = pd.read_csv(fil, sep="\t", header=None)
		except pd.errors.EmptyDataError:
			logger.warning(f"- File '{fil}' is empty: no sites could be read")
		except Exception as e:
			logger.warning(f"- Error occurred while reading '{fil}'")
			raise e

		dataframes.append(df)

	logger.debug("Joining sites from all files")
	sites_table = pd.concat(dataframes, sort=False)
	sites_table.columns = ["Sites_" + str(element) for element in sites_table.columns]
	sites_table.fillna("", inplace=True) #replace NaN with empty string
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

	#Match source columns
	source_col_tfbs = sites_table.columns[3]	#Source id (TF name) is the 4th column of the sites bedfile 
	source_col_origin = [match[1] for match in sorted_matching if match[0] == source_col_tfbs][0]  	#source id column (TF name) in the origin table
	
	#Remove tfbs->origin columns from list
	sorted_matching = [tup for tup in sorted_matching if not ((tup[0] == source_col_tfbs) or (tup[1] == source_col_origin))]
	
	#Find match-pair for targets
	target_col_tfbs = [match[0] for match in sorted_matching if match[0] != source_col_tfbs][0]	 	#target id (gene id) column from sites
	target_col_origin = [match[1] for match in sorted_matching if (match[0] != source_col_tfbs)][0]	#target id (gene id) in the origin table
	
	logger.debug("From TFBS table: source_col='{0}' | target_col='{1}'".format(source_col_tfbs, target_col_tfbs))
	logger.debug("From origin table: source_col='{0}' | target_col='{1}'".format(source_col_origin, target_col_origin))
	
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
		logger.warning("The following source ids (4th column) from '--TFBS' could not be found in the '--origin' table: {0}".format(missing_ids))

		#Subset sites_table to those with source within common_ids
		n_rows = sites_table.shape[0]
		sites_table = sites_table[sites_table[source_col_tfbs].isin(common_ids)]
		logger.warning("Subset {0} sites to {1} sites with a valid source id in the '--origin' table".format(n_rows, sites_table.shape[0]))

	#Remove sites without targets (NAN)
	n_rows = sites_table.shape[0]
	sites_table = sites_table[sites_table[target_col_tfbs] != ""]	#targets not empty
	logger.info("Subset {0} sites to {1} sites with any target given".format(n_rows, sites_table.shape[0]))

	#Subset sites_table to targets within target_ids_origin (through origin table) - e.g. only keep targets which are TFs
	n_rows = sites_table.shape[0]
	valid_targets = origin_table[target_col_origin]
	sites_table = sites_table[sites_table[target_col_tfbs].isin(valid_targets)]
	logger.info("Subset {0} sites to {1} sites with matching target id in '--origin'".format(n_rows, sites_table.shape[0]))

	#Merge sites with origin table (to get target motif ids)
	n_rows = sites_table.shape[0]
	sites_table_convert = sites_table.merge(origin_table[[source_col_origin, target_col_origin]], left_on=target_col_tfbs, right_on=target_col_origin, how="inner")
	logger.info("Merged sites/targets with '--origin' table: Continuing with {0} TF-target connections".format(sites_table_convert.shape[0]))
	if n_rows < sites_table_convert.shape[0]:
		msg = "NOTE: The number of TF-target connections is higher than the number of unique sites, which occurs when the '--origin' table contains "
		msg += "target ids assigned to multiple motifs. In this case, CreateNetwork will treat each motif as an independent TF in the graph. "
		logger.info(msg)
	
	#Subset to unique edges if chosen
	#if args.unique == True:
	#	sites_table_convert.drop_duplicates(subset=[source_col_tfbs, source_col_origin], inplace=True)
	#	logger.info("Flag --unique is on: Sites were further subset to {0} unique edges".format(sites_table_convert.shape[0]))

	#------------------------------------ Write out edges / adjacency ----------------------------------------#
	logger.info("")	#create space in logger output

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
