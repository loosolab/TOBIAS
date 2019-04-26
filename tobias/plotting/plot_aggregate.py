#!/usr/bin/env python

"""
Plot aggregate signals from TFBS across different bigwigs

@author: Mette Bentsen
@contact: mette.bentsen (at) mpi-bn.mpg.de
@license: MIT
"""

import os
import sys
import argparse
import logging
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from sklearn import preprocessing

import itertools
from datetime import datetime
import scipy
import sklearn

#Bio-stuff
import pyBigWig
import pybedtools as pb

#Internal classes
from tobias.utils.utilities import *
from tobias.utils.regions import *
from tobias.utils.logger import *


def add_aggregate_arguments(parser):

	#---------------- Parse arguments ---------------#
	parser.formatter_class = lambda prog: argparse.RawDescriptionHelpFormatter(prog, max_help_position=40, width=90)
	description = ""
	parser.description = format_help_description("PlotAggregate", description)

	parser._action_groups.pop()	#pop -h

	IO = parser.add_argument_group('Input / output arguments')
	IO.add_argument('--TFBS', metavar="", nargs="*", help="TFBS sites (*required)")
	IO.add_argument('--signals', metavar="", nargs="*", help="Signals in bigwig format (*required)")
	IO.add_argument('--regions', metavar="", nargs="*", help="Regions to overlap with TFBS", default=[])
	IO.add_argument('--whitelist', metavar="", nargs="*", help="Only plot sites within whitelist (.bed)", default=[])
	IO.add_argument('--blacklist', metavar="", nargs="*", help="Exclude sites within blacklist (.bed)", default=[])
	IO.add_argument('--output', metavar="", help="Path to output (default: TOBIAS_aggregate.pdf)", default="TOBIAS_aggregate.pdf")

	PLOT = parser.add_argument_group('Plot arguments')
	PLOT.add_argument('--title', metavar="", help="Title of plot (default: \"Aggregated signals\")", default="Aggregated signals")
	PLOT.add_argument('--flank', metavar="", help="Flanking basepairs (+/-) to show in plot (counted from middle of motif) (default: 60)", default="60", type=int)
	PLOT.add_argument('--TFBS_labels', metavar="", help="Labels used for each TFBS file (default: prefix of each --TFBS)", nargs="*")
	PLOT.add_argument('--signal_labels', metavar="", help="Labels used for each signal file (default: prefix of each --signals)", nargs="*")
	PLOT.add_argument('--region_labels', metavar="", help="Labels used for each regions file (default: prefix of each --regions)", nargs="*")
	PLOT.add_argument('--share_y', metavar="", help="Share y-axis range across plots (none/signals/sites/both). Use \"--share_y signals\" if bigwig signals have similar ranges. Use \"--share_y sites\" if sites per bigwig are comparable, but bigwigs themselves aren't comparable. (default: none)", choices=["none", "signals", "sites", "both"], default="none")

	#signals / regions
	PLOT.add_argument('--norm_comparisons', action='store_true', help="Normalize the aggregate signal in comparison column/row to the same range (default: the true range is shown)")
	PLOT.add_argument('--negate', action='store_true', help="Negate overlap with regions")
	PLOT.add_argument('--log_transform', help="", action="store_true")
	PLOT.add_argument('--plot_boundaries', help="Plot TFBS boundaries", action='store_true')
	#PLOT.add_argument('--outliers', help="")

	RUN = parser.add_argument_group("Run arguments")
	RUN = add_logger_args(RUN)

	return(parser)


def run_aggregate(args):
	""" Function to make aggregate plot given input from args """

	#########################################################################################
	############################## Setup logger/input/output ################################
	#########################################################################################

	logger = TobiasLogger("PlotAggregate", args.verbosity)
	logger.begin()

	parser = add_aggregate_arguments(argparse.ArgumentParser())
	logger.arguments_overview(parser, args)
	logger.output_files([args.output])

	#Check input parameters
	check_required(args, ["TFBS", "signals"])
	check_files([args.TFBS, args.signals, args.regions, args.whitelist, args.blacklist], action="r")
	check_files([args.output], action="w")
	
	#### Test input ####
	if args.TFBS_labels != None and (len(args.TFBS) != len(args.TFBS_labels)):
		logger.error("ERROR --TFBS and --TFBS_labels have different lengths ({0} vs. {1})".format(len(args.TFBS), len(args.TFBS_labels)))
	if args.region_labels != None and (len(args.regions) != len(args.region_labels)):
		logger.error("ERROR: --regions and --region_labels have different lengths ({0} vs. {1})".format(len(args.regions), len(args.region_labels)))
		sys.exit()
	if args.signal_labels != None and (len(args.signals) != len(args.signal_labels)):
		logger.error("ERROR: --signals and --signal_labels have different lengths ({0} vs. {1})".format(len(args.signals), len(args.signal_labels)))
		sys.exit()

	#### Format input ####
	#args.TFBS = [os.path.abspath(f) for f in args.TFBS]
	#args.regions = [os.path.abspath(f) for f in args.regions]
	#args.signals = [os.path.abspath(f) for f in args.signals]
	args.TFBS_labels = [os.path.splitext(os.path.basename(f))[0] for f in args.TFBS] if args.TFBS_labels == None else args.TFBS_labels
	args.region_labels = [os.path.splitext(os.path.basename(f))[0] for f in args.regions] if args.region_labels == None else args.region_labels
	args.signal_labels = [os.path.splitext(os.path.basename(f))[0] for f in args.signals] if args.signal_labels == None else args.signal_labels 
	#args.output = os.path.abspath(args.output)


	#########################################################################################
	############################ Get input regions/signals ready ############################
	#########################################################################################

	logger.info("---- Processing input ----")
	logger.info("Reading information from .bed-files")

	#Make combinations of TFBS / regions
	column_names = []

	if len(args.regions) > 0:
		logger.info("Overlapping sites to --regions")
		regions_dict = {}

		combis = itertools.product(range(len(args.TFBS)), range(len(args.regions)))
		for (i,j) in combis:
			TFBS_f = args.TFBS[i]
			region_f = args.regions[j]

			#Make overlap
			pb_tfbs = pb.BedTool(TFBS_f)
			pb_region = pb.BedTool(region_f)
			#todo: write out lengths

			overlap = pb_tfbs.intersect(pb_region, u=True)
			#todo: length after overlap

			name = args.TFBS_labels[i] + " <OVERLAPPING> " + args.region_labels[j]	#name for column 
			column_names.append(name)
			regions_dict[name] = RegionList().from_bed(overlap.fn)

			if args.negate == True:
				overlap_neg = pb_tfbs.intersect(pb_region, v=True)

				name = args.TFBS_labels[i] + " <NOT OVERLAPPING> " + args.region_labels[j]
				column_names.append(name)
				regions_dict[name] = RegionList().from_bed(overlap_neg.fn)

	else:
		column_names = args.TFBS_labels
		regions_dict = {args.TFBS_labels[i]: RegionList().from_bed(args.TFBS[i]) for i in range(len(args.TFBS))}

		for name in regions_dict:
			logger.stats("COUNT {0}: {1} sites".format(name, len(regions_dict[name]))) #length of RegionList obj


	#-------- Do overlap of regions if whitelist / blacklist -------#
	if len(args.whitelist) > 0 or len(args.blacklist) > 0:
		logger.info("Subsetting regions on whitelist/blacklist")
		for regions_id in regions_dict:
			sites = pb.BedTool(regions_dict[regions_id].as_bed(), from_string=True)
			logger.stats("Found {0} sites in {1}".format(len(regions_dict[regions_id]), regions_id))
			
			if len(args.whitelist) > 0:
				for whitelist_f in args.whitelist:
					whitelist = pb.BedTool(whitelist_f)
					sites_tmp = sites.intersect(whitelist, u = True)
					sites = sites_tmp
					logger.stats("Overlapped to whitelist -> {0}".format(len(sites)))

			if len(args.blacklist) > 0:
				for blacklist_f in args.blacklist:
					blacklist = pb.BedTool(blacklist_f)
					sites_tmp = sites.intersect(blacklist, v = True)
					sites = sites_tmp
					logger.stats("Removed blacklist -> {0}".format(format(len(sites))))

			regions_dict[regions_id] = RegionList().from_bed(sites.fn)

	# Estimate motif width per region
	site_list = regions_dict[list(regions_dict.keys())[0]]
	if len(site_list) > 0:
		motif_width = site_list[0].get_width()
	else:
		motif_width = 0

	# Set width (centered on mid)
	args.width = args.flank*2
	for regions_id in regions_dict:
		regions_dict[regions_id].apply_method(OneRegion.set_width, args.width)


	#########################################################################################
	############################ Read signal for bigwig per site ############################
	#########################################################################################

	logger.info("Reading signal from bigwigs")

	signal_dict = {} 
	for i, signal_f in enumerate(args.signals):

		signal_name = args.signal_labels[i]
		signal_dict[signal_name] = {}

		#Open pybw to read signal
		pybw = pyBigWig.open(signal_f, "rb")

		logger.info("- Reading signal from {0}".format(signal_name))

		for regions_id in regions_dict:
			for one_region in regions_dict[regions_id]:
				tup = one_region.tup()	#(chr, start, end, strand)
				if tup not in signal_dict[signal_name]:	#only get signal if it was not already read previously
					signal_dict[signal_name][tup] = one_region.get_signal(pybw) 	#returns signal

		pybw.close()

	


	#########################################################################################
	################################## Calculate aggregates #################################
	#########################################################################################
	
	logger.comment("")
	logger.info("---- Analysis ----")

	#Calculate aggregate per signal/region comparison
	logger.info("Calculating aggregate signals")
	aggregate_dict = {signal_name:{region_name: [] for region_name in regions_dict} for signal_name in args.signal_labels}
	for row, signal_name in enumerate(args.signal_labels):	
		for col, region_name in enumerate(column_names):
			
			signalmat = np.array([signal_dict[signal_name][reg.tup()] for reg in regions_dict[region_name]])

			#Exclude outliers from each column
			lower_limit, upper_limit = -np.inf, np.inf 
			#lower_limit, upper_limit = np.percentile(signalmat[signalmat != 0], [0.5,99.5])
			logical = np.logical_and(np.min(signalmat, axis=1) > lower_limit, np.max(signalmat, axis=1) < upper_limit)
			#logger.debug("Lower limits: {0}".format(lower_limit))
			#logger.debug("Upper limits: {0}".format(upper_limit))
			signalmat = signalmat[logical]
			
			if args.log_transform:
				signalmat_abs = np.abs(signalmat)
				signalmat_log = np.log2(signalmat_abs + 1)
				signalmat_log[signalmat < 0] *= -1	 #original negatives back to <0
				signalmat = signalmat_log

			aggregate = np.nanmean(signalmat, axis=0)
			aggregate_dict[signal_name][region_name] = aggregate

	#Measure of footprint depth in comparison to baseline
	logger.info("Calculating footprint depth measure")
	logger.info("FPD (signal,regions): footprint_width baseline middle FPD")
	for row, signal_name in enumerate(args.signal_labels):	
		for col, region_name in enumerate(column_names):

			agg = aggregate_dict[signal_name][region_name]

			#Estimation of possible footprint width
			FPD_results = []
			for fp_flank in range(int(motif_width/2), min([25, args.flank])):	#motif width for this bed

				#Baseline level
				baseline_indices = list(range(0,args.flank-fp_flank)) + list(range(args.flank+fp_flank, len(agg)))
				baseline = np.mean(agg[baseline_indices])

				#Footprint level
				middle_indices = list(range(args.flank-fp_flank, args.flank+fp_flank))
				middle = np.mean(agg[middle_indices]) 	#within the motif

				#Footprint depth
				depth = middle - baseline
				FPD_results.append([fp_flank*2, baseline, middle, depth])

			#Estimation of possible footprint width
			all_fpds = [result[-1] for result in FPD_results]
			FPD_results_best = FPD_results #[result + ["  "] if result[-1] != min(all_fpds) else result + ["*"] for result in FPD_results]

			for result in FPD_results_best:
				logger.stats("FPD ({0},{1}): {2} {3:.3f} {4:.3f} {5:.3f}".format(signal_name, region_name, result[0], result[1], result[2], result[3]))

	#Compare pairwise to calculate chance of footprint
	logger.comment("")
	logger.info("Calculating pairwise aggregate pearson correlation")
	logger.info("CORRELATION (signal1,region1) VS (signal2,region2): PEARSONR")
	plots = itertools.product(args.signal_labels, column_names)
	combis = itertools.combinations(plots, 2)

	for ax1, ax2 in combis:

		signal1, region1 = ax1
		signal2, region2 = ax2
		agg1 = aggregate_dict[signal1][region1]
		agg2 = aggregate_dict[signal2][region2]

		pearsonr, pval = scipy.stats.pearsonr(agg1, agg2)

		logger.stats("CORRELATION ({0},{1}) VS ({2},{3}): {4:.5f}".format(signal1, region1, signal2, region2, pearsonr))


	#########################################################################################
	################################ Set up plotting grid ###################################
	#########################################################################################

	logger.comment("")
	logger.info("---- Plotting aggregates ----")
	logger.info("Setting up plotting grid")

	no_rows = len(args.signals) + 1 if len(args.signals) > 1 else len(args.signals)
	no_cols = len(column_names) + 1 if len(column_names) > 1 else len(column_names)
	row_compare = True if no_rows > 1 else False
	col_compare = True if no_cols > 1 else False

	#Set grid
	fig, axarr = plt.subplots(no_rows, no_cols, figsize = (no_cols*5, no_rows*5))
	axarr = np.array(axarr).reshape((-1, 1)) if no_cols == 1 else axarr		#Fix indexing for one column figures
	axarr = np.array(axarr).reshape((1, -1)) if no_rows == 1 else axarr		#Fix indexing for one row figures

	#Title of plot and grid
	plt.suptitle(args.title, fontsize=16)

	row_names = args.signal_labels + ["Comparison"] if row_compare else args.signal_labels
	col_names = column_names + ["Comparison"] if col_compare else column_names

	for col in range(no_cols):
		axarr[0,col].set_title(col_names[col].replace(" ","\n"))

	for row in range(no_rows):
		axarr[row,0].set_ylabel(row_names[row], fontsize=12)

	#Colors
	colors = mpl.cm.brg(np.linspace(0, 1, len(args.signals) + len(column_names)))

	#xvals
	flank = int(args.width/2.0)
	xvals = np.arange(-flank,flank+1)
	xvals = np.delete(xvals, flank)

	#Settings for each subplot
	for row in range(no_rows):
		for col in range(no_cols):
			axarr[row, col].set_xlim(-flank, flank)
			axarr[row, col].set_xlabel('bp from center')
			#axarr[row, col].set_ylabel('Mean aggregated signal')

			#Motif boundaries
			if args.plot_boundaries:
				width = motif_width
				mstart = - np.floor(width/2.0)
				mend = np.ceil(width/2.0) - 1 #as it spans the "0" nucleotide
				axarr[row, col].axvline(mstart, color="grey", linestyle="dashed", linewidth=1)
				axarr[row, col].axvline(mend, color="grey", linestyle="dashed", linewidth=1)
				
			minor_ticks = np.arange(-flank, flank, args.width/10.0)


	#Settings for comparison plots
	a = [axarr[-1, col].set_facecolor("0.9") if row_compare == True else 0 for col in range(no_cols)]
	a = [axarr[row, -1].set_facecolor("0.9") if col_compare == True else 0 for row in range(no_rows)]

	#Air between subplots
	plt.subplots_adjust(hspace=0.3, top=0.93)



	#########################################################################################
	############################## Fill in with bigwig scores ###############################
	#########################################################################################

	for row, signal_name in enumerate(args.signal_labels):
		for col, region_name in enumerate(column_names):

			logger.info("Plotting regions {0} from signal {1}".format(region_name, signal_name))

			#If there are any regions:
			if len(regions_dict[region_name]) > 0:
				
				#Signal in region
				aggregate = aggregate_dict[signal_name][region_name] 
				axarr[row, col].plot(xvals, aggregate, color=colors[col+row], linewidth=1, label=signal_name)

				aggregate_norm = preprocessing.minmax_scale(aggregate)

				#Compare across rows and cols
				if col_compare: 	#compare between different regions
					if args.norm_comparisons:
						aggregate_compare = aggregate_norm
					else:
						aggregate_compare = aggregate		
					axarr[row, -1].plot(xvals, aggregate_compare, color=colors[row+col], linewidth=1, alpha=0.8, label=region_name)

				if row_compare:	#compare between different bigwigs
					if args.norm_comparisons:
						aggregate_compare = aggregate_norm
					else:
						aggregate_compare = aggregate
					axarr[-1, col].plot(xvals, aggregate_compare, color=colors[row+col], linewidth=1, alpha=0.8, label=signal_name)
					axarr[-1, col].legend(loc="lower right")

				#Diagonal comparison
				if no_rows == no_cols and col_compare and row_compare and col == row:
					axarr[-1, -1].plot(xvals, aggregate_compare, color=colors[row+col], linewidth=1, alpha=0.8, label=signal_name)

				#Add number of sites to plot
				axarr[row, col].text(0.98,0.98,str(len(regions_dict[region_name])), transform = axarr[row,col].transAxes, fontsize=12, va="top", ha="right")


	#------------- Finishing up plots ---------------#

	logger.info("Adjusting final details")

	#remove lower-right corner if not applicable
	if no_rows != no_cols:
		axarr[-1,-1] = None
	
	if args.share_y == "none":
		pass 	#Do not set ylim for plots

	#Signals are comparable (for example normalized signal between two conditions)
	elif args.share_y == "signals":	
		for col in range(no_cols):
			lims = np.array([ax.get_ylim() for ax in axarr[:,col] if ax is not None])
			ymin, ymax = np.min(lims), np.max(lims)

			for row in range(no_rows):
				if axarr[row, col] is not None:
					axarr[row, col].set_ylim(ymin, ymax)

	#Regions are comparable (for example bound/unbound)
	elif args.share_y == "sites":
		for row in range(no_rows):
			lims = np.array([ax.get_ylim() for ax in axarr[row,:] if ax is not None])
			ymin, ymax = np.min(lims), np.max(lims)
			for col in range(no_cols):
				if axarr[row, col] is not None:
					axarr[row, col].set_ylim(ymin, ymax)

	#Comparable on both rows/columns
	elif args.share_y == "both":
		global_ymin, global_ymax = np.inf, -np.inf
		for row in range(no_rows):
			for col in range(no_cols):
				if axarr[row, col] is not None:
					local_ymin, local_ymax = axarr[row, col].get_ylim()
					global_ymin = local_ymin if local_ymin < global_ymin else global_ymin
					global_ymax = local_ymax if local_ymax > global_ymax else global_ymax


		for row in range(no_rows):
			for col in range(no_cols):
				if axarr[row, col] is not None:
					axarr[row, col].set_ylim(global_ymin, global_ymax)
	
	plt.savefig(args.output, bbox_inches='tight')
	plt.close()

	logger.end()

#--------------------------------------------------------------------------------------------------------#
if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser = add_aggregate_arguments(parser)
	args = parser.parse_args()

	if len(sys.argv[1:]) == 0:
		parser.print_help()
		sys.exit()
		
	run_aggregate(args)
