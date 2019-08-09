#!/usr/bin/env python

"""
CreateNetwork: Creates a TF-TF gene regulation network from transcription factor binding sites

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
#--------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------#

def add_network_arguments(parser):

	parser.formatter_class = lambda prog: argparse.RawDescriptionHelpFormatter(prog, max_help_position=40, width=90)
	description = "" 
	parser.description = format_help_description("CreateNetwork", description)

	parser._action_groups.pop()	#pop -h
	
	#Required arguments
	required = parser.add_argument_group('Required arguments')
	required.add_argument('--TFBS', metavar="", help="File(s) containing TFBS to create network from", nargs="*")
	required.add_argument('--origin', metavar="", help="File containing gene origins of TF <-> gene")

	#Optional arguments
	optional = parser.add_argument_group('Optional arguments')
	optional.add_argument('--start', metavar="", help="Name of node to start in")
	optional.add_argument('--max-len', default=3)
	#optional.add_argument('--additional', metavar="", help="Additional information on genes to add; for example RNA-seq")
	#optional.add_argument('--subset', metavar="", help="File containing subset of names to filter on")
	#optional.add_argument('--TFBS_columns', metavar="", nargs="*", help="Source TF -> target gene columns", default=["TFBS_name", "gene_id"])
	#optional.add_argument('--origin_columns', metavar="", nargs="*", help="Source TF, source gene", default=["TOBIAS_motif_id", "ensembl_gene_id"])

	optional.add_argument('--output', metavar="", help="Path to output directory (default: tobias_network)", default="tobias_network") 

	return(parser)

#--------------------------------------------------------------------------------#
def dfs(adjacency, path, all_paths = [], options={"max_length":4}):


	last_node = path[-1] 
	if last_node in adjacency:
		target_nodes = adjacency[last_node].get("targets", [])
		if len(target_nodes) > 0:

			#Go through all targets and get paths
			for target_node in target_nodes:
				if target_node in path:    #cyclic
					new_path = path + [target_node]
					all_paths += [new_path]
				else:
					new_path = path + [target_node]
					
					if len(new_path) == options["max_length"]:
						all_paths += [new_path]			#Save current path
					else:
						all_paths = dfs(adjacency, new_path, all_paths, options)		#Find targets of last node in path
		else:
			all_paths += [path] 

	else:
		all_paths += [path]

	return all_paths


#--------------------------------------------------------------------------------#
def run_network(args):

	make_directory(args.output)
	check_required(args, ["TFBS", "origin"])

	#-------------------------- Origin file translating motif name -> gene origin -----------------------------------#
	#translation file, where one motif can constitute more than one gene (jun::fos) 
	#and more genes can encode transcription factors with same motifs (close family members with same target sequence)
	origin_table = pd.read_csv(args.origin, sep="\t")
	print(origin_table)

	#------------------------ Transcription factor binding sites with mapping to target genes -----------------------#

	print("Reading all binding sites to one file")
	#todo: read in parallel
	dataframes = []
	for fil in args.TFBS:
		print("- {0}".format(fil))

		df = pd.read_csv(fil, sep="\t", header=None)
		dataframes.append(df)

	print("Joining dataframes")
	sites_table = pd.concat(dataframes, sort=False)
	print("Total of {0} sites found".format(sites_table.shape[0]))


	#------------------------------------- Match target columns to origin ------------------------------------#
	print("Matching target genes back to TFs")
	origin_table_str_columns = list(origin_table.dtypes.index[origin_table.dtypes == "object"])
	sites_table_str_columns = list(sites_table.dtypes.index[sites_table.dtypes == "object"])

	origin_table = origin_table.apply(lambda x: x.astype(str).str.upper())
	sites_table = sites_table.apply(lambda x: x.astype(str).str.upper())

	#Establishing matched columns
	matching = []
	for sites_column in sites_table_str_columns:		
		sites_column_content = set(sites_table[sites_column])

		for origin_column in origin_table_str_columns:
			origin_column_content = set(origin_table[origin_column])

			#Overlap
			overlap = len(origin_column_content & sites_column_content)
			matching.append((sites_column, origin_column, overlap))

	sorted_matching = sorted(matching, key = lambda tup: -tup[-1])

	#Columns for matching
	source_id = 3
	source_id_origin = [match[1] for match in sorted_matching if match[0] == 3][0] 	#source id in the origin table
	target_id = [match[0] for match in sorted_matching if match[0] != 3][0]			#target id in from sites
	target_id_origin = [match[1] for match in sorted_matching if match[0] != 3][0]	#target id in the origin table

	print(source_id_origin)
	print(target_id)
	print(target_id_origin)
	
	#Intersect of sources and targets
	all_source_ids = set(sites_table[source_id])
	all_target_ids = set(origin_table[source_id_origin])
	common_ids = all_source_ids & all_target_ids
	print(common_ids)

	#Add target motif id to sites
	origin_table_sub = origin_table[origin_table[source_id_origin].isin(common_ids)]
	sites_table_convert = sites_table.merge(origin_table_sub[[source_id_origin, target_id_origin]], left_on=target_id, right_on=target_id_origin)


	#--------------------- --------------------#

	##### Write out edges #####
	edges_f = os.path.join(args.output, "edges.txt")
	sites_table_convert.to_csv(edges_f, sep="\t", index=False)

	###### Create adjacency list ####
	print("Creating adjacency mat")
	adjacency = {source: {"targets":[]} for source in common_ids}
	for index, row in sites_table_convert.iterrows():
		source, target = row[source_id], row[source_id_origin]
		if target not in adjacency[source]["targets"]:
			adjacency[source]["targets"].append(target)

	print("Writing out adjacency list")
	adjacency_f = os.path.join(args.output, "adjacency.txt")
	with open(adjacency_f, "w") as f:
		for source in sorted(adjacency):
			f.write("{0}\t{1}\n".format(source, ", ".join(adjacency[source]["targets"])))

	#Create possible paths through graph
	print("Create possible paths through graph")
	paths_f = os.path.join(args.output, "paths.txt")
	paths_out = open(paths_f, "w")
	paths_out.write("Regulatory_path\tn_nodes\tn_timepoints\n")

	#Start node
	if args.start != None:
		start_nodes = [one_id for one_id in common_ids if args.start.upper() in one_id]
		print("Starting nodes are: {0}".format(start_nodes))
	else:
		start_nodes = common_ids
	
	#Find paths per source_id
	for source_id in start_nodes:
		print("Paths starting from {0}".format(source_id))

		paths = dfs(adjacency = adjacency, path = [source_id], options={"max_length": args.max_len})

		print("Writing paths")
		path_edges = []
		for path in paths:
			
			#String formatting of path
			str_paths = ""
			for i, node in enumerate(path[:-1]):
				str_paths += "{0} --> ".format(node)
			str_paths += path[-1]

			n_nodes = len(path)
			paths_out.write("{0}\t{1}\n".format(str_paths, n_nodes))

			#Make pairwise edges across path
			path_edges += ["\t".join([path[i], path[j], str(j)]) for i,j in zip(range(0, len(path)-1), range(1,len(path)))]
		
		source_paths_f = os.path.join(args.output, "{0}_path_edges.txt".format(source_id))
		with open(source_paths_f, "w") as f:
			f.write("Source\tTarget\tLevel\n")
			for path in set(path_edges):
				f.write(path + "\n")
			f.close()
		
	paths_out.close()

	print("done")
