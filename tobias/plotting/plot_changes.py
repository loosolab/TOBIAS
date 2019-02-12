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
import matplotlib.pyplot as plt
import numpy as np
import itertools

from tobias.utils.utilities import *
from tobias.utils.logger import *

#-------------------------------------------------------------------------------------------------#
def add_plotchanges_arguments(parser):

	parser.formatter_class = lambda prog: argparse.RawDescriptionHelpFormatter(prog, max_help_position=35, width=90)
	description = "PlotChanges is a utility to plot the changes in TF binding across multiple conditions as predicted by TOBIAS BINdetect.\n\n"
	description += "Example usage:\n$ echo CTCF GATA > TFS.txt\n$ TOBIAS PlotChanges --bindetect <bindetect_results.txt> --TFS TFS.txt\n\n"

	parser.description = format_help_description("PlotChanges", description)

	parser._action_groups.pop()	#pop -h

	required_arguments = parser.add_argument_group('Required arguments')
	required_arguments.add_argument('--bindetect', metavar="", help='Bindetect_results.txt file from BINDetect run')
	required_arguments.add_argument('--TFS', metavar="", help='Text file containing names of TFs to show in plot (one per line)') 

	#All other arguments are optional
	optional_arguments = parser.add_argument_group('Optional arguments')
	optional_arguments.add_argument('--output', metavar="", help='Output file for plot (default: bindetect_changes.pdf)', default="bindetect_changes.pdf")
	optional_arguments.add_argument('--conditions', metavar="", help="Ordered list of conditions to show (default: conditions are ordered as within the bindetect file)", nargs="*")
	optional_arguments = add_logger_args(optional_arguments)
	
	return(parser)


#-------------------------------------------------------------------------------------------------#
def run_plotchanges(args):

	#------------------------------------ Get ready ------------------------------------#
	logger = TobiasLogger("PlotChanges", args.verbosity)
	logger.begin()

	check_required(args, ["bindetect", "TFS"])
	check_files([args.bindetect, args.TFS], "r")
	check_files([args.output], "w")

	#------------------------------------ Read data ------------------------------------#

	logger.info("Reading data from bindetect file")

	# Read in bindetect file
	table = pd.read_csv(args.bindetect, sep="\t", index_col=0)
	all_TFS = list(table.index)
	logger.info("{0} TFS found in bindetect file".format(len(all_TFS)))

	#Read in TF names from --TFS:
	given_TFS = open(args.TFS, "r").read().split()
	logger.info("TFS given in --TFS: {0}".format(given_TFS))

	#Find matches between all and given
	logger.info("Matching given TFs with bindetect file...")
	lofl = [given_TFS, all_TFS]
	matches = match_lists(lofl)

	for i, TF in enumerate(given_TFS):
		logger.info("- {0} matched with: {1}".format(TF, matches[0][i]))
	
	#Get tfs
	chosen_TFS = list(flatten_list(matches))
	logger.info("Chosen TFS to view in plot: {0}".format(chosen_TFS))

	# Get order of conditions
	header = list(table.columns.values)
	conditions_file = [element.replace("_bound", "") for element in header if "bound" in element]
	if args.conditions == None:	
		args.conditions = conditions_file
	else:
		if not all([z in conditions_file for z in args.conditions]):
			logger.info("ERROR: --conditions {0} is not a subset of bindetect conditions ({1})".format(args.conditions, conditions_file))
			sys.exit()
		
	#condition_comparison = list(itertools.combinations(args.conditions, 2))
	logger.info("Conditions in order: {0}".format(args.conditions))


	#------------------------------------ Make plot --------------------------------#

	logger.info("Plotting figure")
	cmap = matplotlib.cm.get_cmap('rainbow')
	colors = cmap(np.linspace(0,1,len(chosen_TFS)))

	fig, ax1 = plt.subplots(figsize=(10,5))

	#Make lineplot per TF
	for i, TF in enumerate(chosen_TFS):
		no_bound = np.array([table.at[TF, "{0}_bound".format(cond)] for cond in args.conditions])

		"""
		diffs = []
		for (cond1, cond2) in condition_comparison:
			try:
				diffs.append(-table.at[TF, "{0}_{1}_change".format(cond1, cond2)])		#positive means cond1 > cond2, meaning cond1->cond2 change should be negated
			except KeyError:
				diffs.append(table.at[TF, "{1}_{0}_change".format(cond1, cond2)])

		diffs = np.cumsum(diffs)
		"""
		#diffs = [table.at[TF, "{0}_{1}_change".format(cond1, cond2)] for (cond1, cond2) in condition_comparison]

		percent_bound = no_bound / table.at[TF, "total_tfbs"] * 100.0

		#Number of bound sites
		xvals = np.arange(0,len(args.conditions))
		ax1.plot(xvals, percent_bound, color=colors[i], marker="o", label=TF)

		#Annotate
		ax1.annotate(TF, (xvals[0]-0.1, percent_bound[0]), color=colors[i], horizontalalignment="right", verticalalignment="center")

	ax1.set_ylabel("Percent of total sites predicted bound", color="black")
	ax1.tick_params('y', colors='black')

	#Change between conditions
	"""
	ax2 = ax1.twinx()

	xvals_shift = np.arange(0.5,len(condition_comparison),1)
	ax2.plot(xvals_shift, diffs, color="r", marker="o")
	ax2.set_ylabel("Difference between conditions", color="r", rotation=270)
	ax2.tick_params('y', colors='r')
	"""

	#General
	plt.title("Changes in TF binding across conditions")
	plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
	plt.xticks(xvals, args.conditions)
	plt.xlabel("Conditions")

	plt.xlim(xvals[0]-2, xvals[-1]+0.5)
		
	#plt.tight_layout()
	plt.savefig(args.output, format="pdf", bbox_inches='tight')
	#plt.show()


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
