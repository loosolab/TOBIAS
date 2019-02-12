#!/usr/bin/env python

"""
Plot signal heatmaps from TFBS across different bigwigs

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
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.gridspec as gridspec
from datetime import datetime

from sklearn import preprocessing

import pyBigWig
import pysam
import pybedtools as pb

from tobias.utils.regions import *
from tobias.utils.utilities import *


def add_heatmap_arguments(parser):

	#---------------- Parse arguments ---------------#
	parser.formatter_class = lambda prog: argparse.RawDescriptionHelpFormatter(prog, max_help_position=40, width=90)
	description = "PlotHeatmap plots a heatmap of signals from bigwig(s) (each row is one site) as well as the aggregate signal across all sites."
	parser.description = format_help_description("PlotHeatmap", description)
	
	parser._action_groups.pop()	#pop -h
	
	IO = parser.add_argument_group('Input / output arguments')
	IO.add_argument('--TFBS', metavar="", nargs="*", action='append', help="TFBS sites per column (*required)")	#if more than one, set to next column
	IO.add_argument('--signals', metavar="", nargs="*", help="Signals in bigwig format (*required)")
	IO.add_argument('--output',  metavar="", help="Output filename (default: TOBIAS_heatmap.pdf)", default="TOBIAS_heatmap.pdf")

	PLOT = parser.add_argument_group('Plot arguments')
	PLOT.add_argument('--plot_boundaries', help="Plot TFBS boundaries", action='store_true')
	PLOT.add_argument('--share_colorbar', help="Share colorbar across all bigwigs (default: estimate colorbar per bigwig)", action='store_true')
	PLOT.add_argument('--flank', metavar="", help="", type=int, default=75)
	
	PLOT.add_argument('--title', metavar="", default="TOBIAS heatmap")
	PLOT.add_argument('--TFBS_labels', metavar="", nargs="*", action='append', help="Labels of TFBS (default: basename of --TFBS)")
	PLOT.add_argument('--signal_labels', metavar="", nargs="*", help="Labels of signals (default: basename of --signals)")

	PLOT.add_argument('--show_columns', nargs="*", metavar="", type=int, help="Show scores from TFBS column besides heatmap. Set to 0-based python coordinates (for example -1 for last column) (default: None)", default=[])
	PLOT.add_argument('--sort_by', metavar="", help="Columns in .bed to sort heatmap by (default: input .beds are not sorted)", type=int)

	RUN = parser.add_argument_group('Run arguments')
	RUN = add_logger_args(RUN)

	return(parser)


#----------------------------------------------------------------------------------------#
def run_heatmap(args):

	#Start logger
	logger = TobiasLogger("PlotHeatmap", args.verbosity)
	logger.begin()

	parser = add_heatmap_arguments(argparse.ArgumentParser())
	logger.arguments_overview(parser, args)
	logger.output_files([args.output])

	check_required(args, ["TFBS", "signals"])
	
	#Setup TFBS names if not yet
	if args.TFBS_labels == None:
		args.TFBS_labels = [[os.path.basename(fil) for fil in args.TFBS[i]] for i in range(len(args.TFBS))] 

	if args.signal_labels == None:
		args.signal_labels = [os.path.basename(fil) for fil in args.signals]


	########################################################

	#Check valid input parameters (number of input TFBS vs. bigwig etc.)
	no_signals = len(args.signals)
	no_columns = len(args.show_columns)
	no_TFBS_col = len(args.TFBS)

	if no_TFBS_col > 1 and len(args.show_columns) > 0:
		sys.exit("Error: option --show_columns is not available for multiple --TFBS inputs.")

	if no_TFBS_col > 1 and no_signals != no_TFBS_col:
		sys.exit("Error: Number of --TFBS does not match number of signals")

	elif no_TFBS_col == 1 and no_signals > 1:
		#copy bed_f to other columns

		logger.info("Using bedfiles: {0} across all bigwigs".format(args.TFBS))	
		for i in range(no_signals-1):
			args.TFBS.append(args.TFBS[0])
			args.TFBS_labels.append(args.TFBS_labels[0])

	else:
		for i, signal in enumerate(args.signals):
			logger.info("Using {0} with signal from {1}".format(args.TFBS[i], signal))

	#todo: logger overview of bedfiles per column?
	
	######################################################################################
	##################################### INPUT DATA #####################################
	######################################################################################

	#Setup info dict
	heatmap_info = {col:{row:{"bigwig_f": args.signals[col], "bed_f":args.TFBS[col][row]} for row in range(len(args.TFBS[col]))} for col in range(len(args.signals))}

	#Add extra columns
	for i, bed_column in enumerate(args.show_columns):
		heatmap_info[no_signals+i] = {row:{"column": bed_column, "bed_f":args.TFBS[0][row]} for row in range(len(args.TFBS[0]))}


	#------------------------------------------------------------------------------------#
	#------------------------ Read input files to RegionLists ---------------------------#
	#------------------------------------------------------------------------------------#

	seen_bed = []

	#Read regions per heatmap in grid
	logger.comment("")
	logger.info("Reading bedfiles")
	for col in range(len(heatmap_info)):
		for row in range(len(heatmap_info[col])):
		
			heatmap_info[col][row]["regions"] = RegionList().from_bed(heatmap_info[col][row]["bed_f"])

			#Estimate region width
			distri = heatmap_info[col][row]["regions"].get_width_distri()
			if len(distri) > 1:
				sys.exit(distri)
			
			heatmap_info[col][row]["width"] = list(distri.keys())[0]

			#Extend to flank
			heatmap_info[col][row]["regions"] = heatmap_info[col][row]["regions"].apply_method(OneRegion.set_width, 2*args.flank)
						
			#Sort if chosen 
			if args.sort_by != None:
				try:
					heatmap_info[col][row]["regions"].sort(key=lambda region: float(region[args.sort_by]), reverse=True)
				except:
					heatmap_info[col][row]["regions"].sort(key=lambda region: region[args.sort_by], reverse=True)

			#Get scores from file
			invalid = []
			for i, bed_column in enumerate(args.show_columns):
				heatmap_info[no_signals+i][row]["column_{0}".format(bed_column)] = [region[bed_column] for region in heatmap_info[col][row]["regions"]]
				
				try:
					heatmap_info[no_signals+i][row]["column_{0}".format(bed_column)] = [float(element) for element in heatmap_info[no_signals+i][row]["column_{0}".format(bed_column)]]
				except:
					logger.info("Column {0} cannot be converted to float - excluding".format(bed_column))
					del heatmap_info[no_signals+i][row]["column_{0}".format(bed_column)]
					invalid.append(bed_column)
			
			for bed_column in invalid:
				args.show_columns.remove(bed_column)
				


			#Logger info about bedfile
			if heatmap_info[col][row]["bed_f"] not in seen_bed:
				logger.info("- Read {1} sites from {0} of width {2}".format(heatmap_info[col][row]["bed_f"], len(heatmap_info[col][row]["regions"]), heatmap_info[col][row]["width"]))
			seen_bed.append(heatmap_info[col][row]["bed_f"])


	#------------------------------------------------------------------------------------#
	#------------------------------ Signals from all sites ------------------------------#
	#------------------------------------------------------------------------------------#

	logger.comment("")
	logger.info("Reading signals from bigwigs")
	for col in range(len(args.TFBS)):

		bigwig_f = heatmap_info[col][0]["bigwig_f"]		#bigwig is the same for all rows, therefore row == 0
		pybw = pyBigWig.open(bigwig_f, "rb") 			
		
		for row in heatmap_info[col]:
		
			logger.info("- Reading {0} from {1}".format(heatmap_info[col][row]["bed_f"], bigwig_f))

			if len(heatmap_info[col][row]["regions"]) > 0:
				heatmap_info[col][row]["signal_mat"] = np.array([region.get_signal(pybw) for region in heatmap_info[col][row]["regions"]]) 
				heatmap_info[col][row]["aggregate"] = np.mean(heatmap_info[col][row]["signal_mat"], axis=0) 
			else:
				heatmap_info[col][row]["signal_mat"] = None
				heatmap_info[col][row]["aggregate"] = None
		
		pybw.close()

	logger.comment("")

	#------------------------------------------------------------------------------------#
	#---------------------------------- Colorbar min/max --------------------------------#
	#------------------------------------------------------------------------------------#

	#Estimate min/max from all matrices
	if args.share_colorbar == True:

		mats = []
		for col, bigwig in enumerate(args.signals):
			for row in heatmap_info[col]:
				if heatmap_info[col][row]["signal_mat"] is not None:
					mats.append(heatmap_info[col][row]["signal_mat"])
		
		vmin, vmax = (0,0)
		if len(mats) > 0:
			joined = np.vstack(mats)
			vmin, vmax = np.percentile(joined, [1, 99])

		#Set vmin/vmax for all plots
		for col, bigwig in enumerate(args.signals):
			for row in heatmap_info[col]:
				heatmap_info[col][row].update({"vmin":vmin, "vmax":vmax})


	# Estimate min/max for each bigwig
	else:

		for col, bigwig in enumerate(args.signals):
			mats = [heatmap_info[col][row]["signal_mat"] for row in heatmap_info[col] if heatmap_info[col][row]["signal_mat"] is not None]
			vmin, vmax = (0,0)
			if len(mats) > 0:
				joined = np.vstack(mats)

				vmin, vmax = np.percentile(joined, [1, 99])
			
			for row in heatmap_info[col]:
				heatmap_info[col][row].update({"vmin":vmin, "vmax":vmax})

	del mats
	del joined

	# Estimate min/max for extra columns			
	for i, name in enumerate(args.show_columns):
		col = no_signals + i
		glob_values = []
		for row in range(len(args.TFBS[0])):
			glob_values.extend(heatmap_info[col][row]["column_{0}".format(name)])

		vmin, vmax = np.percentile(glob_values, [1, 99])
		for row in range(len(args.TFBS[0])):
			heatmap_info[col][row]["vmin"] = vmin
			heatmap_info[col][row]["vmax"] = vmax

		del glob_values

	######################################################################################
	##################################### PLOTTING #######################################
	######################################################################################


	#------------------------------------------------------------------------------------#
	#------------------------------------ Set up plots ----------------------------------#
	#------------------------------------------------------------------------------------#

	logger.info("Setting up plotting grid")

	total_columns = no_signals + no_columns
	xvals = np.arange(-args.flank, args.flank)

	fig = plt.figure(figsize = (no_signals*5, 5*5))
	h_ratios = [2,10,0.1]
	w_ratios = [1]*no_signals + [0.1]*no_columns
	gs = gridspec.GridSpec(3, total_columns, height_ratios=h_ratios, width_ratios=w_ratios, hspace=0.1, wspace=0.3)	#aggregate + heatmaps (with sub heatmaps) + colorbar

	#Setup axarr fitting to grid
	axdict = {col:{row:"ax" for row in ["aggregate"] + list(heatmap_info[col]) + ["colorbar"]} for col in range(no_signals)}
	axdict.update({col:{row:"ax" for row in ["aggregate"] + list(heatmap_info[col]) + ["colorbar"]} for col in range(no_signals, no_signals+no_columns)})

	#Per signal column
	xvals = np.arange(-args.flank, args.flank)
	for col in range(no_signals):

		#Aggregates
		axdict[col]["aggregate"] = fig.add_subplot(gs[0,col])
		axdict[col]["aggregate"].set_xlim(left=-args.flank, right=args.flank)
		axdict[col]["aggregate"].set_xlabel('bp from center')
		axdict[col]["aggregate"].set_ylabel('Mean aggregate signal')
		axdict[col]["aggregate"].set_title("{0}".format(args.signal_labels[col]))

		#Heatmaps
		no_beds = len(args.TFBS[col]) 
		h_ratios = [len(heatmap_info[col][row]["regions"]) for row in heatmap_info[col]]
		h_ratios = [max(num,1) for num in h_ratios] 	#deal with empty beds
		gs_sub = gridspec.GridSpecFromSubplotSpec(no_beds, 1, subplot_spec=gs[1,col], height_ratios=h_ratios, hspace=0.05)

		for row in range(no_beds):
			axdict[col][row] = plt.Subplot(fig, gs_sub[row,0])
			fig.add_subplot(axdict[col][row])

			#Appearance
			plt.setp(axdict[col][row].get_yticklabels(), visible=False)  #Hide y-axis ticks
			plt.setp(axdict[col][row].get_xticklabels(), visible=False)  #Hide x-axis ticks
			axdict[col][row].tick_params(direction="in")
			axdict[col][row].set_ylabel("{0} ({1})".format(args.TFBS_labels[col][row], len(heatmap_info[col][row]["regions"])))

			#Last row
			if row == no_beds-1:
				axdict[col][row].set_xlabel('bp from center')
		#Colorbar
		axdict[col]["colorbar"] = fig.add_subplot(gs[2,col])	#row number 3

	for col in range(no_signals, no_signals + no_columns):
		gs_sub = gridspec.GridSpecFromSubplotSpec(no_beds, 1, subplot_spec=gs[1,col], height_ratios=h_ratios, hspace=0.05)
		for row in range(no_beds):
			axdict[col][row] = plt.Subplot(fig, gs_sub[row,0])

			plt.setp(axdict[col][row].get_yticklabels(), visible=False)  #Hide y-axis ticks
			plt.setp(axdict[col][row].get_xticklabels(), visible=False)  #Hide x-axis ticks
			axdict[col][row].tick_params(direction="in")
			fig.add_subplot(axdict[col][row])


	#------------------------------------------------------------------------------------#
	#--------------------------------- Fill in plots ------------------------------------#
	#------------------------------------------------------------------------------------#

	logger.info("Filling in grid")

	#Colormaps
	for col, bigwig in enumerate(args.signals):
		colors = mpl.cm.jet(np.linspace(0, 1, len(heatmap_info[col]))) 	#colors for aggregate plots

		for row in heatmap_info[col]:
			if heatmap_info[col][row]["signal_mat"] is not None:
				
				#Aggregate
				axdict[col]["aggregate"].plot(xvals, heatmap_info[col][row]["aggregate"], color=colors[row], linewidth=2, label=args.TFBS_labels[col][row])

				#Heatmap
				lim = np.max([np.abs(heatmap_info[col][row]["vmin"]),np.abs(heatmap_info[col][row]["vmax"])])
				heatmap_info[col][row]["vmin"] = -lim
				heatmap_info[col][row]["vmax"] = lim
				heatmap = axdict[col][row].imshow(heatmap_info[col][row]["signal_mat"], aspect="auto", cmap="seismic", norm=mpl.colors.Normalize(vmin=heatmap_info[col][row]["vmin"], vmax=heatmap_info[col][row]["vmax"]))

				#Insert colorbar (inserted multiple times for each bigwig, but since it is shared for the same bigwig, it doesn't matter)
				fig.colorbar(heatmap, cax=axdict[col]["colorbar"], orientation="horizontal")


	#Extra columns w/ scores from bed
	for i, col in enumerate(range(no_signals, no_signals + no_columns)):
		bed_column = args.show_columns[i]
		for row in heatmap_info[col]:
			values = np.array(heatmap_info[col][row]["column_{0}".format(bed_column)])
			values = values.reshape(-1,1)

			vmin, vmax = np.percentile(values, [1, 99])
			lim = np.max([abs(vmin), abs(vmax)])

			axdict[col][row].imshow(values, aspect="auto", cmap="seismic", norm=mpl.colors.Normalize(vmin=-lim, vmax=lim))


	#------------------------------------------------------------------------------------#
	#-------------------------------- Plot decorations ----------------------------------#
	#------------------------------------------------------------------------------------#

	if args.plot_boundaries:
		for col in heatmap_info:

			motif_len = heatmap_info[col][0]["width"] 
			mstart = int(-np.floor(motif_len/2.0))
			mend = int(np.ceil(motif_len/2.0))

			axdict[col]["aggregate"].axvline(mstart, color="black", linestyle="dashed", linewidth=1)
			axdict[col]["aggregate"].axvline(mend, color="black", linestyle="dashed", linewidth=1)

			for row in heatmap_info[col]:

				motif_len = heatmap_info[col][row]["width"] 
				mstart = int(-np.floor(motif_len/2.0))
				mend = int(np.ceil(motif_len/2.0))

				axdict[col][row].axvline(mstart+args.flank, color="black", linestyle="dashed", linewidth=1)
				axdict[col][row].axvline(mend+args.flank, color="black", linestyle="dashed", linewidth=1)


	#Add legend to aggregate plots
	for col in range(len(args.signals)):
		axdict[col]["aggregate"].legend(loc=1, prop={"size":6})
	
	
	if args.share_colorbar == True:
		ymin = min([axdict[col]["aggregate"].get_ylim()[0] for col in range(no_signals)])
		ymax = max([axdict[col]["aggregate"].get_ylim()[1] for col in range(no_signals)])
		for col in range(no_signals):
			axdict[col]["aggregate"].set_ylim([ymin, ymax])


	#------------------------------------------------------------------------------------#
	#----------------------------- Finish off and output --------------------------------#
	#------------------------------------------------------------------------------------#

	"""
	#For each heatmap
	for row in [1,2]:
		plt.setp(axarr[row].get_yticklabels(), visible=False)  #Hide y-axis ticks
		plt.setp(axarr[row].get_xticklabels(), visible=False)  #Hide x-axis ticks
		axarr[row].tick_params(direction="in")
	"""

	plt.subplots_adjust(top=0.95)
	plt.suptitle(args.title, fontsize=25)

	logger.info("Writing output file")
	plt.savefig(args.output, bbox_inches='tight')
	plt.close()

	logger.end()

#--------------------------------------------------------------------------------------------------------#
if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser = add_heatmap_arguments(parser)
	args = parser.parse_args()

	if len(sys.argv[1:]) == 0:
		parser.print_help()
		sys.exit()
		
	run_heatmap(args)
