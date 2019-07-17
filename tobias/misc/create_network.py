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

def add_network_arguments(parser):

	parser.formatter_class = lambda prog: argparse.RawDescriptionHelpFormatter(prog, max_help_position=40, width=90)
	description = "" 
	parser.description = format_help_description("CreateNetwork", description)

	parser._action_groups.pop()	#pop -h
	
	#Required arguments
	required = parser.add_argument_group('Required arguments')
	required.add_argument('--TFBS', metavar="", help="File(s) containing TFBS to create network from")
	required.add_argument('--origin', metavar="", help="File containing gene origins of TF <-> gene")

	#Optional arguments
	optional = parser.add_argument_group('Optional arguments')
	optional.add_argument('--subset', metavar="", help="File containing subset of names to filter on")
	optional.add_argument('--TFBS_columns', metavar="", nargs="*", help="Source TF -> target gene columns", default=["TFBS_name", "gene_id"])
	optional.add_argument('--origin_columns', metavar="", nargs="*", help="Source TF, source gene", default=["TOBIAS_motif_id", "ensembl_gene_id"])

	#optional.add_argument('--additional', metavar="", help="Additional information on genes to add; for example RNA-seq")
	optional.add_argument('--output', metavar="", help="Path to output directory (default: tobias_network)", default="tobias_network") 

	return(parser)

#--------------------------------------------------------------------------------#
def dfs(adjacency, path, timeline, all_paths = {"paths":[], "timelines":[]}):
	last_node = path[-1] 
	if last_node in adjacency:
		target_nodes = adjacency[last_node].get("targets", [])
		if len(target_nodes) > 0:
			for target_node in target_nodes:
				#print("node {0} has target node {1}".format(last_node, target_node))
				if target_node in path:    #cyclic
					pass  #next target without saving path
				else:
					#if any time overlaps with this or next stage
					allowed = set([timeline[-1], timeline[-1]+1])
					observed = adjacency[last_node]["targets"][target_node]["bound_in"]
					observed_in_allowed = allowed.intersection(observed)
					if len(observed_in_allowed) > 0:
						current_timepoint = max(observed_in_allowed)
						if timeline.count(current_timepoint) < 2: #only two binding events in the same timepoint is possible
							new_path = path + [target_node]
							new_timeline = timeline + [current_timepoint]
							all_paths = dfs(adjacency, new_path, new_timeline, all_paths)
		else:
			if len(path) > 2:
				#print("no targets for {0} - ending path".format(last_node))
				all_paths["paths"] += [path] 
				all_paths["timelines"] += [timeline]
	else:
		if len(path) > 2:
			#print("{0} doesnt have targets".format(last_node))
			all_paths["paths"] += [path]
			all_paths["timelines"] += [timeline]

	return all_paths


#--------------------------------------------------------------------------------#
def run_network(args):

	make_directory(args.output)
	check_required(args, ["TFBS", "origin"])

	#-------------------------- Origin file translating motif name -> gene origin -----------------------------------#
	#translation file, where one motif can constitute more than one gene (jun::fos) 
	#and more genes can encode transcription factors with same motifs (close family members with same target sequence)
	origin_table = pd.read_csv(args.origin, sep="\t")
	#origin_table = origin_table[args.origin_columns]
	
	print(origin_table)

	TFBS_donor, TFBS_recipient = args.TFBS_columns
	origin_name, origin_gene = args.origin_columns

	#id2name = {gene_id: names_list for gene_id, names_list in origin_table.groupby(origin_gene)[origin_name].apply(list).iteritems()}
	
	#---------------------------------------------- BINDetect results --------------------------------------------#

	#Get all overview files from TFBS dir
	print("Getting files from {0}".format(args.TFBS))
	TF_folders = list(os.walk(args.TFBS))[0][1]													#list of folder names = TF names
	overview_files = [f for f in glob.glob(args.TFBS + "/*/*_overview.txt")]
	print("- Found results from {0} motifs".format(len(overview_files)))

	#Subset on names
	if args.subset != None:
		print("Subsetting on names from {0}".format(args.subset))
		subset_names = open(args.subset).read().rstrip().split()
		print("Subset: {0}".format(subset_names))

		overview_files = list(set(sum(match_lists([subset_names, overview_files])[0], [])))
		print("Subset to {0} files".format(len(overview_files)))

	#Read all edges to table
	#find donor col in origin_table
	#origin_name_col = 
	#todo: read in parallel
	#print(overview_files)

	print("Reading all overview tables")
	dataframes = []
	for fil in overview_files[:30]:
		print("- {0}".format(fil))

		df = pd.read_csv(fil, sep="\t")

		#Translating donors to origin genes
		print("Translating to origin genes")
		df = pd.merge(df, origin_table, how="left", left_on=TFBS_donor, right_on=origin_name)
		df.rename(columns={TFBS_donor: "donor_name", origin_gene: "donor_id", origin_name: "origin_table_" + origin_name}, inplace=True)

		#Check if donor could be translated -> otherwise it cannot be included in network
		not_found = df[df["donor_id"].isnull()]
		if not_found.shape[0] > 0:
			print("- Could not find origin gene id for motif {0}".format(not_found["donor_name"][0]))
			continue 

		#Translating recipient ids to names
		df = pd.merge(df, origin_table, how="left", left_on=TFBS_recipient, right_on=origin_gene)
		df.rename(columns={origin_name: "recipient_name", origin_gene: "recipient_id"}, inplace=True)

		#Only include recipients which are also TFS
		print("- Total of {0} sites".format(df.shape[0]))
		df = df[df["recipient_name"].notnull()]
		print("- {0} sites have TF targets".format(df.shape[0]))

		#Find conditions
		bound_columns =  [col for col in df.columns if "_bound" in col]
		condition_names = [col.replace("_bound", "").lower() for col in bound_columns]

		#Remove sites not bound in any condition
		print("- Removing unbound sites")
		df = df[df[bound_columns].any(axis=1)]

		#Rename bound_columns and create "donor_bound_in"
		print("- List of bound conditions")
		df.rename(columns=dict(zip(bound_columns, condition_names)), inplace=True)
		df["donor_bound_in"] = df[condition_names].apply(lambda x: ",".join(x.index[x == 1]), axis=1)

		#Select columns
		#selected = ["TFBS_chr", "TFBS_start", "TFBS_end", "donor_name", "donor_id", "recipient_name", "recipient_id", "donor_bound_in"]
		#df = df[selected]

		#Add to list
		dataframes.append(df)

	print("Joining dataframes")
	sites = pd.concat(dataframes)
	print("Total of {0} sites found".format(sites.shape))


	#------------------------------------- Expression info to subset edges -----------------------------------#
	"""
	print("Reading expression data")
	#Read expression values
	expression_threshold = 50
	if args.expression != None:
		expression_table = pd.read_csv(args.expression, sep="\t", index_col=0)
		expression_table.columns = [col.lower() for col in expression_table.columns]

		#Add expression of target genes
		sites = pd.merge(sites, expression_table, how="left", left_on="recipient_id", right_index=True)

		#set null columns to 0
		sites.fillna(0, inplace=True)
		sites["recipient_expressed_in"] = sites[expression_table.columns].apply(lambda x: ",".join(x.index[x > expression_threshold]), axis=1)
		sites.drop(columns=expression_table.columns, inplace=True)
		expression_dict = expression_table.to_dict()
	else:
		#All are expressed in all conditions
		all_conditions = ",".join(condition_names)
		sites["recipient_expressed_in"] = all_conditions

	#expression_table = "" 		#rows, cols, values > expression_threshold
	expression_dict = expression_table.to_dict()
	"""

	#--------------------------------------------#



	all_source_names = set(sites["donor_name"])
	print(all_source_names)

	print(sites.shape[0])

	sites = sites[(sites["recipient_name"].isin(all_source_names))]
	print(sites.shape[0])

	#--------------------------------------------#


	print(condition_names)
	conditions = condition_names

	#Subset edges on those where donor_bound in is part of recipient_expressed_in
	sites["donor_bound_in"] = sites["donor_bound_in"].apply(lambda x: x.split(","))
	#sites["recipient_expressed_in"] = sites["recipient_expressed_in"].apply(lambda x: x.split(","))

	#Expand the bound_in columns
	exploded = sites["donor_bound_in"].apply(pd.Series).stack().reset_index().rename(columns={0:"donor_bound_in_exploded"})
	sites = pd.merge(sites, exploded, left_index=True, right_on="level_0", how="left").drop(columns=["level_0", "level_1", "donor_bound_in"]).rename(columns={"donor_bound_in_exploded":"donor_bound_in"})

	#print(sites.shape[0])
	#sites = sites[sites.apply(lambda x: x["donor_bound_in"] in x["recipient_expressed_in"], axis=1)]
	#print(sites.shape[0])
	#print(sites)
	
	##### Write out edges
	edges_f = os.path.join(args.output, "edges.txt")
	sites.to_csv(edges_f, sep="\t", index=False)

	###### Create adjacency list ####
	print("Creating adjacency mat")
	all_donor_ids = set(list(sites["donor_id"]))

	adjacency = {donor: {"targets":{}} for donor in all_donor_ids}
	for index, row in sites.iterrows():
		donor, recipient, bound_in = row["donor_id"], row["recipient_id"], row["donor_bound_in"]
		if recipient not in adjacency[donor]["targets"]:
			adjacency[donor]["targets"][recipient] = {"bound_in": []}

		if bound_in not in adjacency[donor]["targets"][recipient]["bound_in"]:
			adjacency[donor]["targets"][recipient]["bound_in"].append(bound_in)

	#Convert donor_bound_in to integer timeline
	#print(adjacency)
	for donor in adjacency:
		for recipient in adjacency[donor]["targets"]:
			adjacency[donor]["targets"][recipient]["bound_in"] = [condition_names.index(cond) for cond in adjacency[donor]["targets"][recipient]["bound_in"]]
	#print(adjacency)

	#Create possible paths through graph
	paths_f = os.path.join(args.output, "paths.txt")
	paths_out = open(paths_f, "w")
	paths_out.write("Regulatory_path\tn_nodes\tn_timepoints\n")
	for gene_id in all_donor_ids:
		print("paths starting from {0}".format(gene_id))

		paths = dfs(adjacency = adjacency, path = [gene_id], timeline = [-1])	#first node must be bound in first condition

		print("Writing paths")
		for i, path in enumerate(paths["paths"]):
			
			timeline = paths["timelines"][i]

			#String formatting of path
			str_paths = ""
			for i, node in enumerate(path[:-1]):
				in_condition = conditions[timeline[i+1]]
				str_paths += "{0} --(Bound: {1}, target_expr: {2})--> ".format(id2name[node][0], expression_dict[node], conditions[in_condition])
			str_paths += id2name[path[-1]][0]

			#Number of nodes
			n_nodes = len(path)
			n_timepoints = len(set(timeline)) - 1 #-1 for the -1-entry in timeline
			paths_out.write("{0}\t{1}\t{2}\n".format(str_paths, n_nodes, n_timepoints))

	paths_out.close()

	#Sort paths
	paths_sorted_f = os.path.join(args.output, "paths_sorted.txt")
	call = "sort -k2,2nr -k3,3nr -t$'\\t' {0} | uniq -u > {1}".format(paths_f, paths_sorted_f)
	print(call)
	os.system(call)

	print("done")

