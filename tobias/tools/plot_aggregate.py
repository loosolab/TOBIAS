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
import numpy as np
import copy
import itertools

import matplotlib as mpl
mpl.use("Agg")	#non-interactive backend
import matplotlib.pyplot as plt
import scipy
import sklearn
from sklearn import preprocessing

#Bio-stuff
import pyBigWig
import pybedtools as pb

#Internal classes
from tobias.parsers import add_aggregate_arguments
from tobias.utils.utilities import check_required, check_files
from tobias.utils.regions import OneRegion, RegionList
from tobias.utils.logger import TobiasLogger, add_logger_args
from tobias.utils.signals import fast_rolling_math

def forceSquare(ax):
	""" Force axes to be square regardless of data limits """
	if ax is not None:
		x0,x1 = ax.get_xlim()
		y0,y1 = ax.get_ylim()
		ax.set_aspect((x1-x0)/(y1-y0))

def fontsize_func(l):
	""" Function to set the fontsize based on the length (l) of the label """

	#Empirically defined thresholds
	lmin = 35
	lmax = 90

	if l < lmin:
		return(12)	#fontsize 12
	elif l > lmax:
		return(5)	#fontsize 5
	else:

		#Map lengths between min/max with linear equation
		p1 = (lmin,12)
		p2 = (lmax,5)

		a = (p2[1] - p1[1]) / (p2[0] - p1[0])
		b = (p2[1] - (a * p2[0]))
		return(a * l + b)
	
#----------------------------------------------------------------------------------------------------#
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
		logger.error("ERROR --TFBS and --TFBS-labels have different lengths ({0} vs. {1})".format(len(args.TFBS), len(args.TFBS_labels)))
		sys.exit(1)
	if args.region_labels != None and (len(args.regions) != len(args.region_labels)):
		logger.error("ERROR: --regions and --region-labels have different lengths ({0} vs. {1})".format(len(args.regions), len(args.region_labels)))
		sys.exit(1)
	if args.signal_labels != None and (len(args.signals) != len(args.signal_labels)):
		logger.error("ERROR: --signals and --signal-labels have different lengths ({0} vs. {1})".format(len(args.signals), len(args.signal_labels)))
		sys.exit(1)

	#### Format input ####
	args.TFBS_labels = [os.path.splitext(os.path.basename(f))[0] for f in args.TFBS] if args.TFBS_labels == None else args.TFBS_labels
	args.region_labels = [os.path.splitext(os.path.basename(f))[0] for f in args.regions] if args.region_labels == None else args.region_labels
	args.signal_labels = [os.path.splitext(os.path.basename(f))[0] for f in args.signals] if args.signal_labels == None else args.signal_labels 

	#TFBS labels cannot be the same
	if len(set(args.TFBS_labels)) < len(args.TFBS_labels):	#this indicates duplicates
		logger.error("ERROR: --TFBS-labels are not allowed to contain duplicates. Note that '--TFBS-labels' are created automatically from the '--TFBS'-files if no input was given." +
					  "Please check that neither contain duplicate names")
		sys.exit(1)

	#########################################################################################
	############################ Get input regions/signals ready ############################
	#########################################################################################

	logger.info("---- Processing input ----")
	logger.info("Reading information from .bed-files")

	#Make combinations of TFBS / regions
	region_names = []

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
			region_names.append(name)
			regions_dict[name] = RegionList().from_bed(overlap.fn)

			if args.negate == True:
				overlap_neg = pb_tfbs.intersect(pb_region, v=True)

				name = args.TFBS_labels[i] + " <NOT OVERLAPPING> " + args.region_labels[j]
				region_names.append(name)
				regions_dict[name] = RegionList().from_bed(overlap_neg.fn)

	else:
		region_names = args.TFBS_labels
		regions_dict = {args.TFBS_labels[i]: RegionList().from_bed(args.TFBS[i]) for i in range(len(args.TFBS))}

		for name in region_names:
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

	# Estimate motif width per --TFBS
	motif_widths = {}
	for regions_id in regions_dict:
		site_list = regions_dict[regions_id]
		if len(site_list) > 0:
			motif_widths[regions_id] = site_list[0].get_width()
		else:
			motif_widths[regions_id] = 0


	#########################################################################################
	############################ Read signal for bigwig per site ############################
	#########################################################################################

	logger.info("Reading signal from bigwigs")

	args.width = args.flank*2	#output regions will be of args.width

	signal_dict = {} 
	for i, signal_f in enumerate(args.signals):

		signal_name = args.signal_labels[i]
		signal_dict[signal_name] = {}

		#Open pybw to read signal
		pybw = pyBigWig.open(signal_f)
		boundaries = pybw.chroms()	#dictionary of {chrom: length}

		logger.info("- Reading signal from {0}".format(signal_name))
		for regions_id in regions_dict:

			original = copy.deepcopy(regions_dict[regions_id])

			# Set width (centered on mid)
			regions_dict[regions_id].apply_method(OneRegion.set_width, args.width)

			#Check that regions are within boundaries and remove if not
			invalid = [i for i, region in enumerate(regions_dict[regions_id]) if region.check_boundary(boundaries, action="remove") == None] 
			for invalid_idx in invalid[::-1]:	#idx from higher to lower
				logger.warning("Region '{reg}' ('{orig}' before flank extension) from bed regions '{id}' is out of chromosome boundaries. This region will be excluded from output.".format(
																									reg=regions_dict[regions_id][invalid_idx].pretty(),
																									orig=original[invalid_idx].pretty(),
																									id=regions_id))
				del regions_dict[regions_id][invalid_idx]

			#Get signal from remaining regions
			for one_region in regions_dict[regions_id]:
				tup = one_region.tup()	#(chr, start, end, strand)
				if tup not in signal_dict[signal_name]:	#only get signal if it was not already read previously
					signal_dict[signal_name][tup] = one_region.get_signal(pybw, logger=logger, key=signal_name) 	#returns signal

		pybw.close()


	#########################################################################################
	################################## Calculate aggregates #################################
	#########################################################################################
	
	signal_names = args.signal_labels

	#Calculate aggregate per signal/region comparison
	logger.info("Calculating aggregate signals")
	aggregate_dict = {signal_name:{region_name: [] for region_name in regions_dict} for signal_name in signal_names}
	for row, signal_name in enumerate(signal_names):	
		for col, region_name in enumerate(region_names):
			
			signalmat = np.array([signal_dict[signal_name][reg.tup()] for reg in regions_dict[region_name]])

			#Check shape of signalmat
			if signalmat.shape[0] == 0: #no regions
				logger.warning("No regions left for '{0}'. The aggregate for this signal will be set to 0.".format(signal_name))
				aggregate = np.zeros(args.width)
			else:	

				#Exclude outlier rows 
				max_values = np.max(signalmat, axis=1)
				upper_limit = np.percentile(max_values, [100*args.remove_outliers])[0]	#remove-outliers is a fraction
				logical = max_values <= upper_limit 
				logger.debug("{0}:{1}\tUpper limit: {2} (regions removed: {3})".format(signal_name, region_name, upper_limit, len(signalmat) - sum(logical)))
				signalmat = signalmat[logical]
							
				#Log-transform values before aggregating
				if args.log_transform:
					signalmat_abs = np.abs(signalmat)
					signalmat_log = np.log2(signalmat_abs + 1)
					signalmat_log[signalmat < 0] *= -1	 #original negatives back to <0
					signalmat = signalmat_log
				
				aggregate = np.nanmean(signalmat, axis=0)

				#normalize between 0-1
				if args.normalize:
					aggregate = preprocessing.minmax_scale(aggregate)

				if args.smooth > 1:
					aggregate_extend = np.pad(aggregate, args.smooth, "edge")
					aggregate_smooth = fast_rolling_math(aggregate_extend.astype('float64'), args.smooth, "mean")
					aggregate = aggregate_smooth[args.smooth:-args.smooth]

			aggregate_dict[signal_name][region_name] = aggregate
			signalmat = None	#free up space

	signal_dict = None #free up space


	#########################################################################################
	############################## Write aggregates to file #################################
	#########################################################################################

	if args.output_txt is not None:

		#Open file for writing
		f_out = open(args.output_txt, "w")
		f_out.write("### AGGREGATE\n")
		f_out.write("# Signal\tRegions\tAggregate\n")
		for row, signal_name in enumerate(signal_names):
			for col, region_name in enumerate(region_names):
				
				agg = aggregate_dict[signal_name][region_name]
				agg_txt = ",".join(["{:.4f}".format(val) for val in agg])

				f_out.write("{0}\t{1}\t{2}\n".format(signal_name, region_name, agg_txt))

		f_out.close()

	#########################################################################################
	################################## Footprint measures ###################################
	#########################################################################################
	
	logger.comment("")
	logger.info("---- Analysis ----")

	#Measure of footprint depth in comparison to baseline
	logger.info("Calculating footprint depth measure")
	logger.info("FPD (signal,regions): footprint_width baseline middle FPD")
	for row, signal_name in enumerate(signal_names):	
		for col, region_name in enumerate(region_names):

			agg = aggregate_dict[signal_name][region_name]
			motif_width = motif_widths[region_name]

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
			FPD_results_best = FPD_results 	#[result + ["  "] if result[-1] != min(all_fpds) else result + ["*"] for result in FPD_results]

			for result in FPD_results_best:
				logger.stats("FPD ({0},{1}): {2} {3:.5f} {4:.5f} {5:.5f}".format(signal_name, region_name, result[0], result[1], result[2], result[3]))

	#Compare pairwise to calculate correlation of signals
	logger.comment("")
	logger.info("Calculating measures for comparing pairwise aggregates")
	logger.info("CORRELATION (signal1,region1) VS (signal2,region2): PEARSONR\tSUM_DIFF")
	plots = itertools.product(signal_names, region_names)
	combis = itertools.combinations(plots, 2)

	for ax1, ax2 in combis:
		signal1, region1 = ax1
		signal2, region2 = ax2
		agg1 = aggregate_dict[signal1][region1]
		agg2 = aggregate_dict[signal2][region2]

		pearsonr, pval = scipy.stats.pearsonr(agg1, agg2)

		diff = np.sum(np.abs(agg1 - agg2))	#Sum of difference between agg1 and agg2

		logger.stats("CORRELATION ({0},{1}) VS ({2},{3}): {4:.5f}\t{5:.5f}".format(signal1, region1, signal2, region2, pearsonr, diff))


	#########################################################################################
	################################ Set up plotting grid ###################################
	#########################################################################################

	logger.comment("")
	logger.info("---- Plotting aggregates ----")
	logger.info("Setting up plotting grid")

	n_signals = len(signal_names) 
	n_regions = len(region_names) #regions are set of sites

	signal_compare = True if n_signals > 1 else False
	region_compare = True if n_regions > 1 else False
	
	#Define whether signal is on x/y
	if args.signal_on_x:
		
		#x-axis
		n_cols = n_signals
		col_compare = signal_compare
		col_names = signal_names

		#y-axis
		n_rows = n_regions 
		row_compare = region_compare
		row_names = region_names
	else:
		#x-axis
		n_cols = n_regions
		col_compare = region_compare
		col_names = region_names 

		#y-axis
		n_rows = n_signals
		row_compare = signal_compare
		row_names = signal_names 

	#Compare across rows/cols?
	if row_compare:
		n_rows += 1
		row_names += ["Comparison"]
	if col_compare:
		n_cols += 1
		col_names += ["Comparison"]

	#Set grid
	fig, axarr = plt.subplots(n_rows, n_cols, figsize = (n_cols*5, n_rows*5), constrained_layout=True)
	axarr = np.array(axarr).reshape((-1, 1)) if n_cols == 1 else axarr		#Fix indexing for one column figures
	axarr = np.array(axarr).reshape((1, -1)) if n_rows == 1 else axarr		#Fix indexing for one row figures

	#X axis / Y axis labels
	#mainax = fig.add_subplot(111, frameon=False)
	#mainax.set_xlabel("X label", labelpad=30, fontsize=16)
	#mainax.set_ylabel("Y label", labelpad=30, fontsize=16)
	#mainax.xaxis.set_label_position('top') 

	#Title of plot and grid
	plt.suptitle(" "*7 + args.title, fontsize=16)	#Add a little whitespace to center the title on the plot; not the frame

	#Titles per column
	for col in range(n_cols):
		title = col_names[col].replace(" ","\n")
		l = max([len(line) for line in title.split("\n")]) 	#length of longest line in title
		s = fontsize_func(l) 								#decide fontsize based on length
		axarr[0, col].set_title(title, fontsize=s)

	#Titles (ylabels) per row
	for row in range(n_rows):
		label = row_names[row]
		l = max([len(line) for line in label.split("\n")])
		axarr[row, 0].set_ylabel(label, fontsize=fontsize_func(l))

	#Colors
	colors = mpl.cm.brg(np.linspace(0, 1, len(signal_names) + len(region_names)))

	#xvals
	flank = int(args.width/2.0)
	xvals = np.arange(-flank,flank+1)
	xvals = np.delete(xvals, flank)

	#Settings for each subplot
	for row in range(n_rows):
		for col in range(n_cols):
			axarr[row, col].set_xlim(-flank, flank)
			axarr[row, col].set_xlabel('bp from center')
			#axarr[row, col].set_ylabel('Mean aggregated signal')
			minor_ticks = np.arange(-flank, flank, args.width/10.0)

	#Settings for comparison plots
	a = [axarr[-1, col].set_facecolor("0.9") if row_compare == True else 0 for col in range(n_cols)]
	a = [axarr[row, -1].set_facecolor("0.9") if col_compare == True else 0 for row in range(n_rows)]


	#########################################################################################
	####################### Fill in grid with aggregate bigwig scores #######################
	#########################################################################################

	for si in range(n_signals):
		signal_name = signal_names[si]
		for ri in range(n_regions):
			region_name = region_names[ri]

			logger.info("Plotting regions {0} from signal {1}".format(region_name, signal_name))

			row, col = (ri, si) if args.signal_on_x else (si, ri)

			#If there are any regions:
			if len(regions_dict[region_name]) > 0:
				
				#Signal in region
				aggregate = aggregate_dict[signal_name][region_name] 
				axarr[row, col].plot(xvals, aggregate, color=colors[col+row], linewidth=1, label=signal_name)

				#Compare across rows and cols
				if col_compare: 	#compare between different columns by adding one more column	
					axarr[row, -1].plot(xvals, aggregate, color=colors[row+col], linewidth=1, alpha=0.8, label=col_names[col])

					s = min([ax.title.get_fontproperties()._size for ax in axarr[0,:]])	#smallest fontsize of all columns
					axarr[row, -1].legend(loc="lower right", fontsize=s)

				if row_compare:	#compare between different rows by adding one more row

					axarr[-1, col].plot(xvals, aggregate, color=colors[row+col], linewidth=1, alpha=0.8, label=row_names[row])

					s = min([ax.yaxis.label.get_fontproperties()._size for ax in axarr[:,0]])	#smallest fontsize of all rows
					axarr[-1, col].legend(loc="lower right", fontsize=s)

				#Diagonal comparison
				if n_rows == n_cols and col_compare and row_compare and col == row:
					axarr[-1, -1].plot(xvals, aggregate, color=colors[row+col], linewidth=1, alpha=0.8)

				#Add number of sites to plot
				axarr[row, col].text(0.98, 0.98, str(len(regions_dict[region_name])), transform = axarr[row, col].transAxes, fontsize=12, va="top", ha="right")
				
				#Motif boundaries (can only be compared across sites)
				if args.plot_boundaries:
				
					#Get motif width for this list of TFBS
					width = motif_widths[region_names[min(row,n_rows-2)]] if args.signal_on_x else motif_widths[region_names[min(col,n_cols-2)]]

					mstart = - np.floor(width/2.0)
					mend = np.ceil(width/2.0) - 1 	#as it spans the "0" nucleotide
					axarr[row, col].axvline(mstart, color="grey", linestyle="dashed", linewidth=1)
					axarr[row, col].axvline(mend, color="grey", linestyle="dashed", linewidth=1)

	#------------- Finishing up plots ---------------#

	logger.info("Adjusting final details")

	#remove lower-right corner if not applicable 
	if n_rows != n_cols and n_rows > 1 and n_cols > 1:
		axarr[-1,-1].axis('off')
		axarr[-1,-1] = None

	#Check whether share_y is set
	if args.share_y == "none":
		pass

	#Comparable rows (rowcompare = same across all)
	elif (args.share_y == "signals" and args.signal_on_x == False) or (args.share_y == "sites" and args.signal_on_x == True):	
		for col in range(n_cols):
			lims = np.array([ax.get_ylim() for ax in axarr[:,col] if ax is not None])
			ymin, ymax = np.min(lims), np.max(lims)

			#Set limit across rows for this col
			for row in range(n_rows):
				if axarr[row, col] is not None:
					axarr[row, col].set_ylim(ymin, ymax)

	#Comparable columns (colcompare = same across all)
	elif (args.share_y == "sites" and args.signal_on_x == False) or (args.share_y == "signals" and args.signal_on_x == True):	
		for row in range(n_rows):
			lims = np.array([ax.get_ylim() for ax in axarr[row,:] if ax is not None])
			ymin, ymax = np.min(lims), np.max(lims)

			#Set limit across cols for this row
			for col in range(n_cols):
				if axarr[row, col] is not None:
					axarr[row, col].set_ylim(ymin, ymax)

	#Comparable on both rows/columns
	elif args.share_y == "both":
		global_ymin, global_ymax = np.inf, -np.inf
		for row in range(n_rows):
			for col in range(n_cols):
				if axarr[row, col] is not None:
					local_ymin, local_ymax = axarr[row, col].get_ylim()
					global_ymin = local_ymin if local_ymin < global_ymin else global_ymin
					global_ymax = local_ymax if local_ymax > global_ymax else global_ymax

		for row in range(n_rows):
			for col in range(n_cols):
				if axarr[row, col] is not None:
					axarr[row, col].set_ylim(global_ymin, global_ymax)
	
	#Force plots to be square
	for row in range(n_rows):
		for col in range(n_cols):
			forceSquare(axarr[row, col])

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
