#!/usr/bin/env python

"""
Script to compare two sets of motifs with each other. Generating similarity matrix and a clustered heatmap.
If only one motif file is given, it will be compared with itself.

@author: René Wiegandt, Anastasiia Petrova
@contact: rene.wiegandt (at) mpi-bn.mpg.de, anastasiia.petrova (at) mpi-bn.mpg.de
@license: MIT

"""

import argparse
import sys
import numpy as np
import seaborn as sns
import pandas as pd
import datetime
import os
import yaml
import re
import math
import itertools
import warnings
from Bio import motifs
from matplotlib import pyplot as plt
from packaging import version
import platform

from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
import scipy.spatial.distance as ssd
from collections import defaultdict

from tobias.parsers import add_motifclust_arguments
from tobias.utils.utilities import check_required, check_files, make_directory
from tobias.utils.logger import TobiasLogger
from tobias.utils.motifs import *

#--------------------------------------------------------------------------------------------------------#
def subset_matrix(matrix, name_list_1, name_list_2):
	"""Subsetting matrix into two symetric matrices given two label lists

	Parameter:
	----------
	matrix :  Pandas DataFrame
		Full similarity matrix
	name_list_1 : list
		list of names
	name_list_2 : list
		list of names

	Returns:
	--------
	list:
		list with subsets of matrix
	"""

	sim_matrix_1 = matrix[name_list_1].loc[name_list_1]
	sim_matrix_2 = matrix[name_list_2].loc[name_list_2]
	sim_matrix = matrix[name_list_1].loc[name_list_2]
	return sim_matrix_1, sim_matrix_2, sim_matrix

#--------------------------------------------------------------------------------------------------------#
def get_motif_stats(motifs):
	"""Get the GC content and the length for a list of motifs

	Parameter:
	----------
	motifs : list
		A list of gimmemotifs motif objects.
	
	Return:
	-------
	motif_stats : dict
		A dictionary containing the motif name as key and
		a list with GC content and motif length as value.
	"""

	motif_stats = dict()
	for m in motifs:
		# Remove first row
		m_woh = re.sub('^>.*\n','',m.to_pwm())
		m_pwm = m_woh.split()
		motif_length = len(m)
		gc_sum = 0
		# Get sum of G and C ratios
		for x in range(0, len(m_pwm)-1, 4):
			gc_sum += float(m_pwm[x+1]) + float(m_pwm[x+2])

		# Calc GC ratio of the motif
		gc_content = gc_sum/motif_length
		motif_name = m.id.replace('\t', ' ')
		motif_stats[motif_name] = [round(gc_content,4), motif_length]

	return motif_stats

#--------------------------------------------------------------------------------------------------------#
def write_motif_stats(stats, out_file):
	"""Write stats to file

	Parameter:
	----------
	stats: dict
		dictionary with motif name as key and list with stats as value
	out_file: string
		path and name of out-file
	"""

	with open(out_file, "w") as f: 
		f.write("Motif\tGC content\tMotif length\n")
		for motif_name,stats in stats.items():
			stats_string = "\t".join(map(str, stats))
			f.write("{}\t{}\n".format(motif_name, stats_string))


#--------------------------------------------------------------------------------------------------------#
def scaling(axis_len):
	"""Scaling of plot figure size

	Parameter:
	----------
	axis_len : int
		Number of rows or columns

	Return:
	----------
	scale : float
		scale for axis
	"""

	return 0.5/math.log10(axis_len)


#--------------------------------------------------------------------------------------------------------#
def write_yaml(dictionary, yml_out):
	"""Writes dict to yaml
	
	Parameter:
	----------
	dictionary : dict
		Dictionary containing information to be written
	yml_out : string
		Path to outfile
	"""

	with open(yml_out, 'w') as outfile:
		yaml.dump(dictionary, outfile, default_flow_style=False)

#--------------------------------------------------------------------------------------------------------#
def plot_dendrogram(label, linkage, font_size, out, title, threshold, dpi):
	"""Plot dendrogram with highlighted threshold 
	Parameter:
	----------
	label : list
		List of labels
	linkage : ndarray
		The hierarchical clustering of rows or cols encoded as a linkage matrix.
	font_size : int
		font size
	out : String
		Output path
	title : String
		Plot title
	threshold : float
		dendrogram cluster threshold
	dpi : int
		dpi of plot
	"""

	x = 10.0
	y = x * len(label)/(x*3)    #ensure good aspect ratio
								#set cap on y axis (prevent errors from too large figure)

	plt.figure(figsize=(x, y))
	plt.title(title, fontsize=20)
	plt.axvline(x=threshold, color="red")
	dendrogram(linkage, color_threshold=threshold, labels=label, leaf_font_size=font_size, orientation="right")
	try:
		plt.tight_layout()
		plt.savefig(out, dpi=dpi)

	except ValueError as e:
		print("Skipped plotting of dendrogram.")
		print("Error: " + str(e))

#--------------------------------------------------------------------------------------------------------#
def plot_heatmap(similarity_matrix, out, col_linkage, row_linkage, dpi, x_name, y_name, color, col_clust, row_clust, zscore):
	"""Plots clustered similarity matrix as a heatmap
	Parameter:
	----------
	similarity_matrix :  DataTable
		a DataTable contaning the similarity score
	out : string
		Prefix to output file
	col_linkage : ndarray
		The hierarchical clustering of cols encoded as a linkage matrix.
	row_linkage : ndarray
		The hierarchical clustering of rows encoded as a linkage matrix.
	dpi : int
		dpi of plot
	x_name : string
		label for x-axis
	y_name : string
		label for y-axis
	color : string
		colorpalette
	col_clust : boolean
		if True, cluster the columns
	row_clust : boolean
		if True, cluster the rows
	"""

	vmin, vmax = None, None
	if zscore == "row":
		zs = 0
	elif zscore == "col":
		zs = 1
	else:
		zs = None
		vmin, vmax = 0, 1

	#Establish figsize
	x_len = len(similarity_matrix.columns)
	y_len = len(similarity_matrix.index.values)
	x = 30  
	y = x * y_len/x_len    #Ensure correct aspect ratio
	
	try:
		plot = sns.clustermap(similarity_matrix,
			row_linkage=row_linkage,
			col_linkage=col_linkage,
			col_cluster= True if row_linkage is not None else False, #not col_clust,
			row_cluster= True if col_linkage is not None else False, #not row_clust,
			z_score=zs,
			cmap=color,
			vmin=vmin, 
			vmax=vmax,
			xticklabels=similarity_matrix.columns,
			yticklabels=similarity_matrix.index.values,
			figsize=(x,y))
		plot.ax_heatmap.set_xlabel(x_name)
		plot.ax_heatmap.set_ylabel(y_name)
		plot.savefig(out, bbox_inches='tight', dpi=dpi)

	except ValueError as e:
		print("Skipped plotting of Heatmap.")
		print("Error: " + str(e))

#--------------------------------------------------------------------------------------------------------#
def run_motifclust(args):

	###### Check input arguments ######
	check_required(args, ["motifs"])			#Check input arguments
	check_files([args.motifs])					#Check if files exist
	out_cons_img = os.path.join(args.outdir, "consensus_motifs_img")
	make_directory(out_cons_img)
	out_prefix = os.path.join(args.outdir, args.prefix)

	###### Create logger and write argument overview ######
	logger = TobiasLogger("ClusterMotifs", args.verbosity)
	logger.begin()
	
	parser = add_motifclust_arguments(argparse.ArgumentParser())
	logger.arguments_overview(parser, args)
	#logger.output_files([])

	out_prefix = os.path.join(args.outdir, args.prefix)

	#----------------------------------------- Check for gimmemotifs ----------------------------------------#

	try:
		from gimmemotifs.motif import Motif
		from gimmemotifs.comparison import MotifComparer
		sns.set_style("ticks")	#set style back to ticks, as this is set globally during gimmemotifs import

	except ModuleNotFoundError:
		logger.error("MotifClust requires the python package 'gimmemotifs'. You can install it using 'pip install gimmemotifs' or 'conda install gimmemotifs'.")
		sys.exit(1)

	except ImportError as e: #gimmemotifs was installed, but there was an error during import
		
		pandas_version = pd.__version__
		python_version = platform.python_version()

		if e.name == "collections" and (version.parse(python_version) >= version.parse("3.10.0")): #collections error from norns=0.1.5 and from other packages on python=3.10
			logger.error("Due to package dependency errors, 'TOBIAS ClusterMotifs' is not available for python>=3.10. Current python version is '{0}'. Please downgrade python in order to use this tool.".format(python_version))
			sys.exit(1)

		elif e.name == "pandas.core.indexing" and (version.parse(pandas_version) >= version.parse("1.3.0")):
			logger.error("Package 'gimmemotifs' version < 0.17.0 requires 'pandas' version < 1.3.0. Current pandas version is {0}.".format(pandas_version))
			sys.exit(1) 
		
		else: #other import error
			logger.error("Tried to import package 'gimmemotifs' but failed with error: '{0}'".format(repr(e)))
			logger.error("Traceback:")
			raise e

	except Exception as e:
		logger.error("Tried to import package 'gimmemotifs' but failed with error: '{0}'".format(repr(e)))
		logger.error("Please check that 'gimmemotifs' was successfully installed.")
		sys.exit(1)
	
	#Check gimmemotifs version vs. metric given
	import gimmemotifs
	gimme_version = gimmemotifs.__version__
	if gimme_version == "0.17.0" and args.dist_method in ["pcc", "akl"]:
		logger.warning("The dist_method given ('{0}') is invalid for gimmemotifs version 0.17.0. Please choose another --dist_method. See also: https://github.com/vanheeringen-lab/gimmemotifs/issues/243".format(args.dist_method))
		sys.exit(1)

	#---------------------------------------- Reading motifs from file(s) -----------------------------------#
	logger.info("Reading input file(s)")

	motif_list = MotifList() #list containing OneMotif objects
	motif_dict = {} 		 #dictionary containing separate motif lists per file

	if sys.version_info < (3,7,0): # workaround for deepcopy with python version < 3.5
		copy._deepcopy_dispatch[type(re.compile(''))] = lambda r, _: r

	for f in args.motifs:
		logger.debug("Reading {0}".format(f))

		motif_format = get_motif_format(open(f).read())
		sub_motif_list = MotifList().from_file(f)	#MotifList object

		logger.stats("- Read {0} motifs from {1} (format: {2})".format(len(sub_motif_list), f, motif_format))

		motif_list.extend(sub_motif_list)
		motif_dict[f] = sub_motif_list

	#Check whether ids are unique
	#TODO

	#---------------------------------------- Motif stats ---------------------------------------------------#
	logger.info("Creating matrix statistics")

	gimmemotifs_list = [motif.get_gimmemotif().gimme_obj for motif in motif_list]
	
	#Stats for all motifs
	full_motifs_out = out_prefix + "_stats_motifs.txt"
	motifs_stats = get_motif_stats(gimmemotifs_list)
	write_motif_stats(motifs_stats, full_motifs_out)


	#---------------------------------------- Motif clustering ----------------------------------------------#
	logger.info("Clustering motifs")

	clusters = motif_list.cluster(threshold=args.threshold, metric=args.dist_method, clust_method=args.clust_method)
	logger.stats("- Identified {0} clusters".format(len(clusters)))
	
	#Write out overview of clusters
	cluster_dict = {cluster_id: [motif.get_gimmemotif().gimme_obj.id for motif in clusters[cluster_id]] for cluster_id in clusters}
	cluster_f = out_prefix + "_" + "clusters.yml"
	logger.info("- Writing clustering to {0}".format(cluster_f))
	write_yaml(cluster_dict, cluster_f)

	# Save similarity matrix to file
	matrix_out = out_prefix + "_matrix.txt"
	logger.info("- Saving similarity matrix to the file: " + str(matrix_out))
	motif_list.similarity_matrix.to_csv(matrix_out, sep = '\t')

	#Plot dendrogram
	logger.info("Plotting clustering dendrogram")
	dendrogram_f = out_prefix + "_" + "dendrogram." + args.type     #plot format pdf/png
	plot_dendrogram(motif_list.similarity_matrix.columns, motif_list.linkage_mat, 12, dendrogram_f, "Clustering", args.threshold, args.dpi)


	#---------------------------------------- Consensus motif -----------------------------------------------#
	logger.info("Building consensus motif for each cluster")

	consensus_motifs = MotifList()
	for cluster_id in clusters:
		consensus = clusters[cluster_id].create_consensus(metric=args.dist_method) 	#MotifList object with create_consensus method
		consensus.id = cluster_id if len(clusters[cluster_id]) > 1 else clusters[cluster_id][0].id 	#set original motif id if cluster length = 1

		consensus_motifs.append(consensus)

	#Write out consensus motif file
	out_f = out_prefix + "_consensus_motifs." + args.cons_format
	logger.info("- Writing consensus motifs to: {0}".format(out_f))
	consensus_motifs.to_file(out_f, args.cons_format)

	#Create logo plots
	out_cons_img = os.path.join(args.outdir, "consensus_motifs_img")
	logger.info("- Making logo plots for consensus motifs (output folder: {0})".format(out_cons_img))
	for motif in consensus_motifs:
		filename = os.path.join(out_cons_img, motif.id + "_consensus." + args.type)
		motif.logo_to_file(filename)


	#---------------------------------------- Plot heatmap --------------------------------------------------#

	logger.info("Plotting similarity heatmap")
	logger.info("Note: Can take a while for --type=pdf. Try \"--type png\" for speed up.")
	args.nrc = False
	args.ncc = False
	args.zscore = "None"
	clust_linkage = motif_list.linkage_mat
	similarity_matrix = motif_list.similarity_matrix

	pdf_out = out_prefix + "_heatmap_all." + args.type
	x_label = "All motifs"
	y_label = "All motifs"
	plot_heatmap(similarity_matrix, pdf_out, clust_linkage, clust_linkage, args.dpi, x_label, y_label, args.color, args.ncc, args.nrc, args.zscore)
	
	# Plot heatmaps for each combination of motif files
	comparisons = itertools.combinations(args.motifs, 2)
	for i, (motif_file_1, motif_file_2) in enumerate(comparisons):
		
		pdf_out = out_prefix + "_heatmap" + str(i) +"." + args.type
		logger.info("Plotting comparison of {0} and {1} motifs to the file ".format(motif_file_1, motif_file_2) + str(pdf_out))

		x_label, y_label = motif_file_1, motif_file_2

		#Create subset of matrices for row/col clustering
		motif_names_1 = [motif.get_gimmemotif().gimme_obj.id for motif in motif_dict[motif_file_1]]
		motif_names_2 = [motif.get_gimmemotif().gimme_obj.id for motif in motif_dict[motif_file_2]]

		m1_matrix, m2_matrix, similarity_matrix_sub = subset_matrix(similarity_matrix, motif_names_1, motif_names_2)

		col_linkage = linkage(ssd.squareform(m1_matrix)) if (len(motif_names_1) > 1 and len(motif_names_2) > 1) else None
		row_linkage = linkage(ssd.squareform(m2_matrix)) if (len(motif_names_1) > 1 and len(motif_names_2) > 1) else None

		#Plot similarity heatmap between file1 and file2
		plot_heatmap(similarity_matrix_sub, pdf_out, col_linkage, row_linkage, args.dpi, x_label, y_label, args.color, args.ncc, args.nrc, args.zscore)
	
	# ClusterMotifs finished
	logger.end()


#--------------------------------------------------------------------------------------------------------#
if __name__ == "__main__":
	
	parser = argparse.ArgumentParser()
	parser = add_motifclust_arguments(parser)
	args = parser.parse_args()

	if len(sys.argv[1:]) == 0:
		parser.print_help()
		sys.exit()

	run_motifclust(args)
