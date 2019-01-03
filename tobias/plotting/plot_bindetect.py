#!/usr/bin/env python

import os
import sys
import argparse
import numpy as np
import re

#Clustering
from sklearn.cluster import DBSCAN
import sklearn.preprocessing as preprocessing
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from scipy.spatial.distance import squareform

#Plotting
import matplotlib
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages
from cycler import cycler
import pandas as pd
from matplotlib.lines import Line2D
from adjustText import adjust_text

#Internal functions and classes
from tobias.utils.utilities import *

import warnings

#---------------------------------------------------------------------------------------------------#
def add_diffplot_arguments(parser):

	parser.formatter_class = lambda prog: argparse.RawDescriptionHelpFormatter(prog, max_help_position=40, width=90)
	description = ""
	parser.description = format_help_description("PlotBINDetect", description)
	
	parser._action_groups.pop()	#pop -h

	args = parser.add_argument_group('Arguments')
	args.add_argument('-m', '--matrix', metavar="", help="Distance matrix from TF clustering")
	args.add_argument('-b', '--bindetect', metavar="", help="Differential binding scores from bindetect_results.txt file")
	args.add_argument('-o', '--output', metavar="", help="(default: diff_bind.pdf)", default="diff_bind.pdf")

	args.add_argument('--cluster_threshold', metavar="", help="Clustering threshold (default: 0.5)", type=float, default=0.5)
	args.add_argument('--change_threshold', metavar="", help="Set change threshold (default: None)", type=float, default=None)
	args.add_argument('--pvalue_threshold', metavar="", help="Set pvalue threshold (default: set at 5%% percentile)", type=float, default=None)

	return(parser)


class Tcluster(): #

	def __init__(self):

		self.similarity_mat = None
		self.IDS = None

		#Created by clustering

	def cluster(threshold=0.5):

		self.dist_mat = squareform(similarity_mat)
		self.linkage_mat = linkage(self.dist_mat, "average")
		self.labels = fcluster(linkage_mat, threshold, criterion="distance")



def cluster_matrix(similarity_mat):

	dist_mat = squareform(similarity_mat)
	linkage_mat = linkage(dist_mat, "average")
	labels = fcluster(linkage_mat, 0.5, criterion="distance")

	return(linkage_mat, labels, clusters)

#---------------------------------------------------------------------------------------------------#
def plot_bindetect(IDS, similarity_mat, changes, pvalues, conditions, change_threshold=None, pvalue_threshold=None):
	""" Conditions refer to the order of the fold_change divison, meaning condition1/condition2 """
	warnings.filterwarnings("ignore")

	comparison = conditions
	IDS = np.array(IDS)
	no_IDS = len(IDS)

	names = np.array([name.split("_")[0] for name in IDS]) #without id

	dist_mat = squareform(similarity_mat)
	linkage_mat = linkage(dist_mat, "average")
	labels = fcluster(linkage_mat, 0.5, criterion="distance")

	diff_scores = {IDS[i]:{"change": changes[i], "pvalue": pvalues[i]} for i in range(len(IDS))}
	pvalues = np.array(pvalues)
	xvalues = np.array(changes)

	#Find min pvalue that is not zero -> set 0 values to this level
	minval = np.min(pvalues[np.nonzero(pvalues)])
	pvalues[pvalues == 0] = minval
	yvalues = -np.log10(pvalues)

	if pvalue_threshold is None:
		all_but_null = pvalues[pvalues != minval]
		if len(all_but_null) > 0:
			pvalue_threshold = np.percentile(all_but_null, 5)
		else:
			pvalue_threshold = 0

	if change_threshold is None:
		all_but_null = np.abs(xvalues[xvalues != 0])
		if len(all_but_null) > 0:
			change_threshold = np.percentile(all_but_null, 95)
		else:
			change_threshold = 0

	#>>>>>>>>>>> Clustering 
	#Find clusters below threshold
	clusters = dict(zip(range(no_IDS), [[num] for num in range(no_IDS)]))

	threshold = 0.5
	for i, row in enumerate(linkage_mat):
		ID1 = int(row[0])
		ID2 = int(row[1])
		new = no_IDS + i
		dist = row[2]

		if dist <= threshold:
			clusters[new] = clusters[ID1] + clusters[ID2] + [new]
			del clusters[ID1]
			del clusters[ID2]

	colorlist = ["blue", "green", "red", "orange"]
	node_color = ["black"] * (2*no_IDS-1)
	i = 0
	for cluster in sorted(list(clusters.keys())):
		if len(clusters[cluster]) > 1:
			color = colorlist[i]
			for node in clusters[cluster]:
				node_color[node] = color
			i += 1 

			if i == len(colorlist):
				i = 0

	#<<<<<<<<<<< Clustering
	
	#--------------------------------------- Figure -------------------------------- #
	#Make figure
	no_rows, no_cols = 2,2	
	h_ratios = [1,max(1,no_IDS/25)]
	figsize = (8,10+7*(no_IDS/25))
	
	fig = plt.figure(figsize = figsize)
	gs = gridspec.GridSpec(no_rows, no_cols, height_ratios=h_ratios)
	gs.update(hspace=0.0001, bottom=0.00001, top=0.999999)

	ax1 = fig.add_subplot(gs[0,:])	#volcano
	ax2 = fig.add_subplot(gs[1,0])	#long scatter overview
	ax3 = fig.add_subplot(gs[1,1])  #dendrogram
	
	######### Volcano plot on top of differential values ########
	
	ax1.set_title("BINDetect volcano plot", fontsize=16, pad=20)

	row, col = 0,0
	labels = IDS
	ax1.scatter(xvalues, yvalues, color="black", s=5)

	#Add +/- 10% to make room for labels
	ylim = ax1.get_ylim()
	y_extra = (ylim[1] - ylim[0]) * 0.1
	ax1.set_ylim(ylim[0], ylim[1] + y_extra)

	xlim = ax1.get_xlim()
	x_extra = (xlim[1] - xlim[0]) * 0.1
	lim = np.max([np.abs(xlim[0]-x_extra), np.abs(xlim[1]+x_extra)])
	ax1.set_xlim(-lim, lim)

	x0,x1 = ax1.get_xlim()
	y0,y1 = ax1.get_ylim()
	ax1.set_aspect((x1-x0)/(y1-y0))		#square volcano plot

	#Plot in pvalue threshold as dotted line
	if pvalue_threshold is not 0:
		threshold = -np.log10(pvalue_threshold)
		ax1.axhline(y=threshold, color="grey", linestyle="--", linewidth=0.5)
		ax1.annotate(" {0:.2E}".format(pvalue_threshold), (lim, threshold), horizontalalignment='left',verticalalignment='center', color="grey")

	if change_threshold is not 0:
		ax1.axvline(-abs(change_threshold), color="grey", linestyle="--", linewidth=0.5)
		ax1.axvline(abs(change_threshold), color="grey", linestyle="--", linewidth=0.5)			#Plot thresholds 

	#Decorate plot
	ax1.set_xlabel("Differential binding score")
	ax1.set_ylabel("-log10(pvalue)")


	########### Dendrogram over similarities of TFs #######
	#row, col = 1,1
	dendro_dat = dendrogram(linkage_mat, labels=IDS, no_labels=True, orientation="right", ax=ax3, above_threshold_color="black", link_color_func=lambda k: node_color[k])
	labels = dendro_dat["ivl"]
	ax3.set_xlabel("Transcription factor similarities\n(Clusters below threshold are colored)")

	ax3.set_ylabel("Transcription factor clustering based on TFBS overlap", rotation=270, labelpad=20)
	ax3.yaxis.set_label_position("right")

	#set aspect
	x0,x1 = ax3.get_xlim()
	y0,y1 = ax3.get_ylim()

	ax3.set_aspect(((x1-x0)/(y1-y0)) * no_IDS/10)		#square volcano plot

	########## Differential binding scores per TF ##########
	#row, col = 1,0
	ax2.set_xlabel("Differential binding score\n" + "(" + comparison[1] + r' $\leftarrow$' + r'$\rightarrow$ ' + comparison[0] + ")") #First position in comparison equals numerator in log2fc division
	ax2.xaxis.set_label_position('bottom') 
	ax2.xaxis.set_ticks_position('bottom') 

	no_labels = len(labels)
	ax2.set_ylim(0.5, no_labels+0.5)
	ax2.set_ylabel("Transcription factors")

	ax2.set_yticks(range(1,no_labels+1))
	ax2.set_yticklabels(labels)
	ax2.axvline(0) #Plot line at mean

	#Plot scores per TF
	colors = []
	for y, TF in enumerate(labels):
		
		idx = np.where(IDS == TF)[0][0]
		score = diff_scores[TF]["change"]
		pvalue = diff_scores[TF]["pvalue"]

		#Set size based on pvalue
		if pvalue <= pvalue_threshold or score < -change_threshold or score > change_threshold:
			fill = "full"
		else:
			fill = "none"

		ax2.axhline(y+1, color="grey", linewidth=1)
		ax2.plot(score, y+1, marker='o', color=node_color[idx], fillstyle=fill)
		ax2.yaxis.get_ticklabels()[y].set_color(node_color[idx])

	#Set x-axis ranges
	lim = np.max(np.abs(ax2.get_xlim()))
	ax2.set_xlim((-lim, lim))	#center on 0

	#set aspect
	x0,x1 = ax2.get_xlim()
	y0,y1 = ax2.get_ylim()
	ax2.set_aspect(((x1-x0)/(y1-y0)) * no_IDS/10)		#square volcano plot

	plt.tight_layout()    #tight layout before setting ids in volcano plot

	# Color points and set labels in volcano
	txts = []

	#points with -fc -> blue
	included = np.logical_and(xvalues < 0, np.logical_or(xvalues <= -change_threshold, pvalues <= pvalue_threshold))
	chosen_x, chosen_y, chosen_l = xvalues[included], yvalues[included], names[included]
	ax1.scatter(chosen_x, chosen_y, color="blue", s=4.5)
	for i, txt in enumerate(chosen_l):
		txts.append(ax1.text(chosen_x[i], chosen_y[i], txt, fontsize=7))

	#points with +fc -> red
	included = np.logical_and(xvalues > 0, np.logical_or(xvalues >= change_threshold, pvalues <= pvalue_threshold))
	chosen_x, chosen_y, chosen_l = xvalues[included], yvalues[included], names[included]
	ax1.scatter(chosen_x, chosen_y, color="red", s=4.5)
	for i, txt in enumerate(chosen_l):
		txts.append(ax1.text(chosen_x[i], chosen_y[i], txt, fontsize=7))

	adjust_text(txts, ax=ax1, text_from_points=True, arrowprops=dict(arrowstyle='-', color='black', lw=0.5))  #, expand_text=(0.1,1.2), expand_objects=(0.1,0.1)) #, arrowprops=dict(arrowstyle='-', color='black', lw=0.5))
	
	#Plot custom legend for colors
	legend_elements = [Line2D([0],[0], marker='o', color='w', markerfacecolor="red", label="More bound in {0}".format(conditions[0])),
						Line2D([0],[0], marker='o', color='w', markerfacecolor="blue", label="More bound in {0}".format(conditions[1]))]
	ax1.legend(handles=legend_elements, bbox_to_anchor=(1.05, 1), loc='upper left')

	return(fig)


#################################################################################################
#################################################################################################
#################################################################################################

def run_diffplot(args):

	#Test input
	check_required(args, ["matrix", "bindetect"]) 
	check_files([args.matrix, args.bindetect], "r")
	check_files([args.output], "w")

	#------------------------ Read bindetect scores -------------------------#

	#Read table, get number of comparison plots to make
	table = pd.read_csv(args.bindetect, sep="\t")
	header = list(table)

	comparisons = [col.replace("_change", "") for col in header if col.endswith("change")] 	#list of "cond1_cond2" pairs
	
	diff_scores = {}
	for comparison in comparisons:
		diff_scores[comparison] = {}
		for index, row in table.iterrows():
			TF = row["TF_name"]
			diff_scores[comparison][TF] = {"change": row[comparison + "_change"], "pvalue": row[comparison + "_pvalue"]}

	#-------------------- Read similarity mat from file ---------------------#

	#read matrix w/header
	IDS = open(args.matrix, "r").readline().rstrip().split("\t")
	distance_mat = np.loadtxt(args.matrix, skiprows=1)

	#Check dimensions of mat
	no_rows, no_columns = distance_mat.shape
	if no_columns != no_rows:
		print("ERROR: Wrong matrix format (({0},{1}))".format(no_rows, no_columns))
		sys.exit()

	#Check column ids
	#column_ids = IDS
	no_IDS = len(IDS)

	#Convert mat to numpy
	#distance_mat = np.array([[float(element) for element in columns[1:]] for columns in mat[1:]])

	#---------------------- Subset similarity mat if needed ------------------#
	# needed when bindetect file has been subset to only some rows

	#Columns/rows to keep
	comparison = list(diff_scores.keys())[0]
	diff_ids = list(diff_scores[comparison].keys())

	to_keep = [i for i in range(no_IDS) if IDS[i] in diff_ids]

	IDS = [IDS[i] for i in to_keep]
	no_IDS = len(IDS)
	distance_mat = distance_mat[np.array(to_keep),:][:,np.array(to_keep)]	 #select rows, then columns
	

	#----------------------------- Make plots -------------------------------#

	#Set up multipage pdf
	figure_pdf = PdfPages(args.output, keep_empty=True)

	#IDS / similarity same across all
	for comparison in comparisons:
	
		changes = [diff_scores[comparison][TF]["change"] for TF in IDS]
		pvalues = [diff_scores[comparison][TF]["pvalue"] for TF in IDS]
		conditions = comparison.split("_")
		print(conditions)

		fig = plot_bindetect(IDS, distance_mat, changes, pvalues, conditions, change_threshold=args.change_threshold, pvalue_threshold=args.pvalue_threshold)
		figure_pdf.savefig(fig,  bbox_inches='tight')

	figure_pdf.close()

	print("Plot written to: {0}".format(args.output))