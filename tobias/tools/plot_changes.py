#!/usr/bin/env python
"""
PlotChanges: Make a plot of changes in TF binding across different conditions (for example timeseries)

@author: Mette Bentsen
@contact: mette.bentsen (at) mpi-bn.mpg.de
@license: MIT

"""

import sys
import os
import argparse
import pandas as pd
import matplotlib
matplotlib.use("Agg")	#non-interactive backend
import matplotlib.pyplot as plt
import numpy as np
import itertools
from matplotlib.backends.backend_pdf import PdfPages

from tobias.parsers import add_plotchanges_arguments
from tobias.utils.utilities import *
from tobias.utils.logger import TobiasLogger

#-------------------------------------------------------------------------------------------------#
def run_plotchanges(args):

	#------------------------------------ Get ready ------------------------------------#
	logger = TobiasLogger("PlotChanges", args.verbosity)
	logger.begin()

	check_required(args, ["bindetect"])
	check_files([args.bindetect, args.TFS], "r")
	check_files([args.output], "w")

	#------------------------------------ Read data ------------------------------------#

	logger.info("Reading data from bindetect file")

	# Read in bindetect file
	bindetect = pd.read_csv(args.bindetect, sep="\t")
	bindetect.set_index("output_prefix", inplace=True, drop=False)
	
	all_TFS = list(bindetect["output_prefix"])
	logger.info("{0} TFS found in bindetect file".format(len(all_TFS)))

	#Read in TF names from --TFS:
	if args.TFS != None:
	
		given_TFS = open(args.TFS, "r").read().split()
		logger.info("TFS given in --TFS: {0}".format(given_TFS))

		#Find matches between all and given
		logger.info("Matching given TFs with bindetect file...")
		lofl = [given_TFS, all_TFS]
		matches = match_lists(lofl)

		for i, TF in enumerate(given_TFS):
			logger.info("- {0} matched with: {1}".format(TF, matches[0][i]))
		
		#Get tfs
		chosen_TFS = list(set(flatten_list(matches)))
		logger.info("Chosen TFS to view in plot: {0}".format(chosen_TFS))
	else:
		logger.info("Showing all TFS in plot. Please use --TFS to subset output.")
		chosen_TFS = all_TFS

	# Get order of conditions
	header = list(bindetect.columns.values)
	conditions_file = [element.replace("_bound", "") for element in header if "bound" in element]
	if args.conditions == None:	
		args.conditions = conditions_file
	else:
		if not all([z in conditions_file for z in args.conditions]):
			logger.info("ERROR: --conditions {0} is not a subset of bindetect conditions ({1})".format(args.conditions, conditions_file))
			sys.exit()
		
	logger.info("Conditions in order: {0}".format(args.conditions))

	#------------------------------------ Make plots --------------------------------#

	logger.info("Plotting figure")

	fig_out = os.path.abspath(args.output)
	figure_pdf = PdfPages(fig_out, keep_empty=True)

	#Changes over time for different measures
	for cluster_flag in [False, True]:
		#logger.info("- Use clusters: {0}".format(cluster_flag))

		#Choose whether to show individual TFs or clusters
		if cluster_flag == True:
			table = bindetect.loc[chosen_TFS,].groupby("cluster").mean() #mean of each column
		else:
			table = bindetect.loc[chosen_TFS]
		
		#Get colors ready
		cmap = matplotlib.cm.get_cmap('rainbow')
		colors = cmap(np.linspace(0,1,len(table)))

		xvals = np.arange(0,len(args.conditions))
		for measure in ["n_bound", "percent_bound", "mean_score"]:
			#logger.info("-- {0}".format(measure))
			fig, ax = plt.subplots(figsize=(10,5))
			for i, TF in enumerate(table.index):

				if measure == "n_bound":
					yvals = np.array([table.at[TF, "{0}_bound".format(cond)] for cond in args.conditions])
				elif measure == "percent_bound":
					n_bound = np.array([table.at[TF, "{0}_bound".format(cond)] for cond in args.conditions])
					yvals = n_bound / table.at[TF, "total_tfbs"] * 100.0	#percent bound
				elif measure == "mean_score":
					yvals = np.array([table.at[TF, "{0}_mean_score".format(cond)] for cond in args.conditions])

				ax.plot(xvals, yvals, color=colors[i], marker="o", label=TF)
				ax.annotate(TF, (xvals[0]-0.1, yvals[0]), color=colors[i], horizontalalignment="right", verticalalignment="center", fontsize=6)

			#General
			title = "Changes in TF binding across conditions"
			if cluster_flag == True:
				title += " (for TF clusters)"
			plt.title(title)
			plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=3, markerscale=0.5)
			plt.xticks(xvals, args.conditions)
			plt.xlabel("Conditions")
			
			if measure == "n_bound":
				plt.ylabel("Number of sites predicted bound", color="black")
			elif measure == "percent_bound":
				plt.ylabel("Percent of sites predicted bound", color="black")
			elif measure == "mean_score":
				plt.ylabel("Mean binding score", color="black")
			ax.tick_params('y', colors='black')

			plt.xlim(xvals[0]-0.5, xvals[-1]+0.25)	#make room for names
			figure_pdf.savefig(fig, bbox_inches='tight')
			plt.close()


	#Change between conditions
	#condition_comparison = list(itertools.combinations(args.conditions, 2))
	"""
		diffs = []
		for (cond1, cond2) in condition_comparison:
			try:
				diffs.append(-table.at[TF, "{0}_{1}_change".format(cond1, cond2)])		#positive means cond1 > cond2, meaning cond1->cond2 change should be negated
			except KeyError:
				diffs.append(table.at[TF, "{1}_{0}_change".format(cond1, cond2)])

		diffs = np.cumsum(diffs)
	
	ax2 = ax1.twinx()

	xvals_shift = np.arange(0.5,len(condition_comparison),1)
	ax2.plot(xvals_shift, diffs, color="r", marker="o")
	ax2.set_ylabel("Difference between conditions", color="r", rotation=270)
	ax2.tick_params('y', colors='r')
	"""

	figure_pdf.close()
	logger.end()
	logger.info("Saved figure to {0}".format(args.output))

#--------------------------------------------------------------------------------------------------------#
if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser = add_plotchanges_arguments(parser)
	args = parser.parse_args()

	if len(sys.argv[1:]) == 0:
		parser.print_help()
		sys.exit()

	run_plotchanges(args)
