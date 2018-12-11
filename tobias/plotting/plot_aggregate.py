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
import pyBigWig
import logging
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from sklearn import preprocessing
import pybedtools as pb
import itertools
from datetime import datetime

from tobias.utils.utilities import *
from tobias.utils.regions import *
#from footprinting.signals import *


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
	IO.add_argument('--whitelist', metavar="", help="Only plot sites within whitelist (.bed)")
	IO.add_argument('--blacklist', metavar="", help="Exclude sites within blacklist (.bed)")
	IO.add_argument('--output', metavar="", help="Path to output (default: TOBIAS_aggregate.pdf)", default="TOBIAS_aggregate.pdf")

	OPT = parser.add_argument_group('Plot options')
	OPT.add_argument('--negate', action='store_true', help="Negate overlap with regions")
	OPT.add_argument('--title', metavar="", help="Title of plot", default="Aggregated signals")
	OPT.add_argument('--width', metavar="", help="", type=int, default=150)
	OPT.add_argument('--region_names', metavar="", nargs="*")
	OPT.add_argument('--signal_names', metavar="", nargs="*")
	OPT.add_argument('--share_y', metavar="", help="Share y-axis range across plots (none/rows/cols/both) (default: none)", choices=["none", "rows", "cols", "both"], default="none")
	#OPT.add_argument('--normalize')

	#OPT.add_argument('--norm_signals')
	OPT.add_argument('--log_transform', help="", action="store_true")
	OPT.add_argument('--norm_regions', help="Normalize aggregate across regions", action='store_true')
	OPT.add_argument('--norm_signals', help="Normalize aggregate across signals", action='store_true')
	OPT.add_argument('--plot_boundaries', help="Plot TFBS boundaries", action='store_true')
	#todo:negate regions

	#RUN = parser.add_argument_group("Run options")
	#RUN.add_argument('-l', '--log', metavar="", help="Full path of logfile (default: logfile is written to stdout)")
	#RUN.add_argument('--silent', help="Prevents info printing to stdout", action='store_true')

	return(parser)


def run_aggregate(args):

	begin_time = datetime.now()

	check_required(args, ["TFBS", "signals"])
	
	#### Test input ####
	if args.region_names != None and (len(args.regions) != len(args.region_names)):
		print("ERROR: --regions and --region_names have different lengths")
		#print(len(args.regions))
		sys.exit()
	if args.signal_names != None and (len(args.signals) != len(args.signal_names)):
		print("ERROR: --signals and --signal_names have different lengths.")
		sys.exit()


	#### Format input ####
	args.TFBS = [os.path.abspath(f) for f in args.TFBS]
	args.regions = [os.path.abspath(f) for f in args.regions]
	args.signals = [os.path.abspath(f) for f in args.signals]
	args.TFBS_names = [os.path.splitext(os.path.basename(f))[0] for f in args.TFBS]
	args.region_names = [os.path.splitext(os.path.basename(f))[0] for f in args.regions] if args.region_names == None else args.region_names 
	args.signal_names = [os.path.splitext(os.path.basename(f))[0] for f in args.signals] if args.signal_names == None else args.signal_names 
	#args.output = os.path.abspath("signal_aggregate.pdf") if args.output == None else os.path.abspath(args.output)


	#########################################################################################
	##################################### Logger info #######################################
	#########################################################################################

	# Create logger
	logger = create_logger()

	#Print info on run
	logger.comment("#TOBIAS PlotAggregate (run started {0})\n".format(begin_time))
	logger.comment("#Command line call: {0}\n".format(" ".join(sys.argv)))

	parser = add_aggregate_arguments(argparse.ArgumentParser())
	logger.comment(arguments_overview(parser, args))
	

	#########################################################################################
	############################ Get input regions/signals ready ############################
	#########################################################################################

	#Make combinations of TFBS / regions
	column_names = []

	if len(args.regions) > 0:
		logger.info("Overlapping regions")
		regions_dict = {}

		combis = itertools.product(range(len(args.TFBS)),range(len(args.regions)))
		for (i,j) in combis:
			TFBS_f = args.TFBS[i]
			region_f = args.regions[j]

			#Make overlap
			pb_tfbs = pb.BedTool(TFBS_f)
			pb_region = pb.BedTool(region_f)
			overlap = pb_tfbs.intersect(pb_region, u=True)

			name = args.TFBS_names[i] + " <OVERLAPPING> " + args.region_names[j]
			column_names.append(name)
			regions_dict[name] = RegionList().from_bed(overlap.fn)

			if args.negate == True:
				overlap_neg = pb_tfbs.intersect(pb_region, v=True)

				name = args.TFBS_names[i] + " <NOT OVERLAPPING> " + args.region_names[j]
				column_names.append(name)
				regions_dict[name] = RegionList().from_bed(overlap_neg.fn)

	else:
		column_names = args.TFBS_names
		regions_dict = {args.TFBS_names[i]: RegionList().from_bed(args.TFBS[i]) for i in range(len(args.TFBS))}



	#-------- Do overlap of regions if whitelist / blacklist -------#
	if args.whitelist != None or args.blacklist != None:
		logger.info("Subsetting regions...")
		for regions_id in regions_dict:
			sites = pb.BedTool(regions_dict[regions_id].as_bed(), from_string=True)
			logger.info("Found {0} sites in {1}".format(len(regions_dict[regions_id]), regions_id))
			
			if args.whitelist != None:
				whitelist = pb.BedTool(args.whitelist)
				sites_tmp = sites.intersect(whitelist, u = True)
				sites = sites_tmp
				logger.info("Overlapped to whitelist -> {0}".format(len(sites)))

			if args.blacklist != None:
				blacklist = pb.BedTool(args.blacklist)
				sites_tmp = sites.intersect(blacklist, v = True)
				sites = sites_tmp
				logger.info("Removed blacklist -> {0}".format(format(len(sites))))

			regions_dict[regions_id] = RegionList().from_bed(sites.fn)

	# Estimate motif width
	site_list = regions_dict[list(regions_dict.keys())[0]]
	if len(site_list) > 0:
		motif_width = site_list[0].get_width()
	else:
		motif_width = 0

	# Set width (centered on mid)
	for regions_id in regions_dict:
		regions_dict[regions_id].apply_method(OneRegion.set_width, args.width)

	# Print number of sites
	for regions_id in regions_dict:
		logger.info("{0}: {1}".format(regions_id, len(regions_dict[regions_id])))



	#########################################################################################
	############################ Read signal for bigwig per site ############################
	#########################################################################################

	signal_dict = {} 
	for i, signal_f in enumerate(args.signals):

		signal_name = args.signal_names[i]
		signal_dict[signal_name] = {}

		#Open pybw to read signal
		pybw = pyBigWig.open(signal_f, "rb")

		logger.info("Reading signal from {0}".format(signal_name))

		for regions_id in regions_dict:
			for one_region in regions_dict[regions_id]:
				tup = one_region.tup()	#(chr, start, end, strand)
				if tup not in signal_dict[signal_name]:
					signal_dict[signal_name][tup] = one_region.get_signal(pybw) 	#returns signal

		pybw.close()


	#########################################################################################
	################################ Set up plotting grid ###################################
	#########################################################################################

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

	row_names = args.signal_names + ["Comparison"] if row_compare else args.signal_names
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

	for row, signal_name in enumerate(args.signal_names):
		for col, region_name in enumerate(column_names):

			logger.info("Plotting regions {0} from signal {1}".format(region_name, signal_name))

			#If there are any regions:
			if len(regions_dict[region_name]) > 0:
				
				#Signal in region
				signalmat = np.array([signal_dict[signal_name][reg.tup()] for reg in regions_dict[region_name]])

				if args.log_transform:
					signalmat_abs = np.abs(signalmat)
					signalmat_log = np.log2(signalmat_abs + 1)
					signalmat_log[signalmat < 0] *= -1	 #original negatives back to <0
					signalmat = signalmat_log

				aggregate = np.nanmean(signalmat, axis=0)
				axarr[row, col].plot(xvals, aggregate, color=colors[col+row], linewidth=1, label=signal_name)

				aggregate_norm = preprocessing.minmax_scale(aggregate)

				#Compare across rows and cols
				if col_compare: 	#compare between different regions
					if args.norm_regions:
						aggregate_compare = aggregate_norm
					else:
						aggregate_compare = aggregate		
					axarr[row, -1].plot(xvals, aggregate_compare, color=colors[row+col], linewidth=1, alpha=0.8, label=region_name)

				if row_compare:	#compare between different bigwigs
					if args.norm_signals:
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

	#remove lower-right corner if not applicable
	if no_rows != no_cols:
		axarr[-1,-1] = None
	
	if args.share_y == "none":
		pass 	#Do not set ylim for plots

	elif args.share_y == "cols":
		for col in range(no_cols):
			lims = np.array([ax.get_ylim() for ax in axarr[:,col] if ax is not None])
			ymin, ymax = np.min(lims), np.max(lims)

			for row in range(no_rows):
				if axarr[row, col] is not None:
					axarr[row, col].set_ylim(ymin, ymax)

	elif args.share_y == "rows":
		for row in range(no_rows):
			lims = np.array([ax.get_ylim() for ax in axarr[row,:] if ax is not None])
			ymin, ymax = np.min(lims), np.max(lims)
			for col in range(no_cols):
				if axarr[row, col] is not None:
					axarr[row, col].set_ylim(ymin, ymax)

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


#--------------------------------------------------------------------------------------------------------#
if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser = add_aggregate_arguments(parser)
	args = parser.parse_args()

	if len(sys.argv[1:]) == 0:
		parser.print_help()
		sys.exit()
		
	run_aggregate(args)
