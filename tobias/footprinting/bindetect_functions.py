#!/usr/bin/env python

"""
BINDetect_functions: Functions to be called from main BINDetect script

@author: Mette Bentsen
@contact: mette.bentsen (at) mpi-bn.mpg.de
@license: MIT

"""

import numpy as np
import pandas as pd
import scipy
from datetime import datetime
import itertools
import xlsxwriter
import random

#Plotting
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import NullFormatter
from cycler import cycler
from matplotlib.lines import Line2D
from adjustText import adjust_text
from scipy.optimize import curve_fit

#Bio-specific packages
import pyBigWig
import pysam
import MOODS.scan
import MOODS.tools
import MOODS.parsers

#Internal functions and classes
from tobias.utils.regions import *
from tobias.utils.sequences import *
from tobias.utils.utilities import *
from tobias.utils.motifs import *
from tobias.utils.signals import *

import warnings
#np.seterr(all='raise')

#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#

def sigmoid(x, a, b, L, shift):
	""" a is the xvalue at the sigmoid midpoint """

	y = L / (1 + np.exp(-b*(x-a))) + shift
	return y

class ArrayNorm:

	def __init__(self, popt):
		self.popt = popt

	def normalize(self, arr):
		return(arr * sigmoid(arr, *self.popt))

def dict_to_tab(dict_list, fname, chosen_columns, header=False):

	#Establish header
	if header == True:
		out_str = "\t".join(chosen_columns) + "\n"
	else:
		out_str = ""
	
	#Add lines
	out_str += "\n".join(["\t".join([str(line_dict[column]) for column in chosen_columns]) for line_dict in dict_list]) + "\n"
	
	#Write file
	f = open(fname, "w")
	f.write(out_str)
	f.close()

#Quantile normalization 
def quantile_normalization(list_of_arrays): #lists paired values to normalize

	n = len(list_of_arrays) #number of arrays to normalize

	#Calculate 
	quantiles = np.linspace(0.2,0.999,500)
	array_quantiles = [np.quantile(arr[arr > 0], quantiles) for arr in list_of_arrays]
	mean_array_quantiles = [np.mean([array_quantiles[i][j] for i in range(n)]) for j in range(len(quantiles))]

	norm_objects = []
	for i in range(n):
		
		#Plot q-q
		#f, ax = plt.subplots()
		#plt.scatter(array_quantiles[i], mean_array_quantiles)
		#ax.set_xlim(0,1)
		#ax.set_ylim(0,1)
		#ax.plot([0, 1], [0, 1], transform=ax.transAxes)
		#plt.close()
		#plt.show()
		
		#Plot normalizyation factors
		#f, ax = plt.subplots()
		xdata = array_quantiles[i]
		ydata = mean_array_quantiles/xdata
		#plt.scatter(xdata, ydata)
		#plt.close()

		popt, pcov = curve_fit(sigmoid, xdata, ydata, bounds=((0,-np.inf, 0, 0), (np.inf, np.inf, np.inf, np.inf)))
		norm_objects.append(ArrayNorm(popt))

		y_est = sigmoid(xdata, *popt)
		plt.plot(xdata, y_est)

		plt.close()

	#Normalize each array using means per rank
	normed = []
	for i in range(n):	#for each array
		normed.append(norm_objects[i].normalize(list_of_arrays[i]))

	return((normed, norm_objects))

def plot_score_distribution(list_of_arr, labels=[], title="Score distribution"):

	fig, ax = plt.subplots(1, 1)
	xlim = []
	for i in range(len(list_of_arr)):

		values = np.array(list_of_arr[i])
		x_max = np.percentile(values, [99]) 
		values = values[values < x_max]
		xlim.append(x_max)

		plt.hist(values, bins=100, alpha=.4, density=True, label=labels[i])

	ax.set_xlabel("Scores")
	ax.set_ylabel("Density")
	ax.set_xlim(0, min(xlim))
	plt.legend()
	plt.title(title)

	return(fig)


def get_gc_content(regions, fasta):
	""" Get GC content from regions in fasta """
	nuc_count = {"T":0, "t":0, "A":0, "a":0, "G":1, "g":1, "C":1, "c":1}

	gc = 0
	total = 0
	fasta_obj = pysam.FastaFile(fasta)
	for region in regions:
		seq = fasta_obj.fetch(region.chrom, region.start, region.end)
		gc += sum([nuc_count.get(nuc, 0.5) for nuc in seq])
		total += region.end - region.start
	fasta_obj.close()
	gc_content = gc / float(total)

	return(gc_content)


#---------------------------------------------------------------------------------------------------------#
#------------------------------------------- Main functions ----------------------------------------------#
#---------------------------------------------------------------------------------------------------------#
def scan_and_score(regions, motifs_obj, args, log_q, qs):
	""" Scanning and scoring runs in parallel for subsets of regions """
	
	logger = TobiasLogger("", args.verbosity, log_q)	#sending all logger calls to log_q

	logger.debug("Setting up scanner/bigwigs/fasta")
	motifs_obj.setup_moods_scanner()	#MotifList object

	pybw = {condition:pyBigWig.open(args.signals[i], "rb") for i, condition in enumerate(args.cond_names)}
	fasta_obj = pysam.FastaFile(args.genome)
	chrom_boundaries = dict(zip(fasta_obj.references, fasta_obj.lengths))

	rand_window = 200

	background_signal = {"gc":[], "signal":{condition:[] for condition in args.cond_names}}

	######## Scan for motifs in each region ######
	logger.debug("Scanning for motif occurrences")
	all_TFBS = {TF: RegionList() for TF in motifs_obj.names} 	# Dict for saving sites before writing
	for i, region in enumerate(regions):
		logger.spam("Processing region: {0}".format(region.tup()))
	
		extra_columns = region
		
		#Check whether region is within boundaries
		if region.end >= chrom_boundaries[region.chrom]:
			print("ERROR IN REGION: {0}".format(region))

		#Random positions for sampling
		reglen = region.get_length()
		random.seed(reglen)		#Each region is processed identifically regardless of order in file
		rand_positions = random.sample(range(reglen), max(1,int(reglen/rand_window)))		#theoretically one in every 500 bp
		logger.spam("Random indices: {0} for region length {1}".format(rand_positions, reglen))

		#Read footprints in region
		footprints = {}
		for condition in args.cond_names:
			try:
				footprints[condition] = pybw[condition].values(region.chrom, region.start, region.end, numpy=True)
				footprints[condition] = np.nan_to_num(footprints[condition])	#nan to 0
			except:
				logger.error("Error reading footprints from region: {0}".format(region))
				continue

			if len(footprints[condition]) == 0:
				logger.error("ERROR IN REGION: {0}".format(region))

			#Read random positions for background
			for pos in rand_positions:
				background_signal["signal"][condition].append(footprints[condition][pos])

		#Scan for motifs across sequence from fasta
		seq = fasta_obj.fetch(region.chrom, region.start, region.end)
		region_TFBS = motifs_obj.scan_sequence(seq, region)		#RegionList of TFBS

		#Extend all TFBS with extra columns from peaks and bigwigs 
		extra_columns = region
		for TFBS in region_TFBS:
			motif_length = TFBS.end - TFBS.start 
			pos = TFBS.start - region.start + int(motif_length/2.0) #middle of TFBS
			
			TFBS.extend(extra_columns)

			#Assign scores from bigwig
			for bigwig in args.cond_names:
				bigwig_score = footprints[bigwig][pos]
				TFBS.append("{0:.5f}".format(bigwig_score))

		#Split regions to single TFs
		for TFBS in region_TFBS:
			all_TFBS[TFBS.name].append(TFBS)

	####### All input regions have been scanned #######
	global_TFBS = RegionList()	#across all TFs

	#Sent sites to writer
	for name in all_TFBS:	
		all_TFBS[name] = all_TFBS[name].resolve_overlaps()
		no_sites = len(all_TFBS[name])

		logger.spam("Sending {0} sites from {1} to bed-writer queue".format(no_sites, name))
		bed_content = all_TFBS[name].as_bed()	#string 
		qs[name].put((name, bed_content))

		global_TFBS.extend(all_TFBS[name])
		all_TFBS[name] = []

	overlap = global_TFBS.count_overlaps()

	#Close down open file handles
	fasta_obj.close()
	for bigwig_f in pybw:
		pybw[bigwig_f].close()
			
	logger.stop()
	logger.total_time

	return(background_signal, overlap)


#-----------------------------------------------------------------------------------------------#
def process_tfbs(TF_name, args, log2fc_params): 	#per tf
	""" Processes single TFBS to split into bound/unbound and write out overview file """

	#begin_time = datetime.now()
	logger = TobiasLogger("", args.verbosity, args.log_q) 	#sending all logger calls to log_q

	#Pre-scanned sites to read
	bed_outdir = os.path.join(args.outdir, TF_name, "beds")
	filename = os.path.join(bed_outdir, TF_name + ".tmp")
	no_cond = len(args.cond_names)
	comparisons = args.comparisons

	#Read file to list of dicts
	header = ["TFBS_chr", "TFBS_start", "TFBS_end", "TFBS_name", "TFBS_score", "TFBS_strand"] + args.peak_header_list + ["{0}_score".format(condition) for condition in args.cond_names]
	with open(filename) as f:
		bedlines = [dict(zip(header, line.rstrip().split("\t"))) for line in f.readlines()]
	n_rows = len(bedlines)

	############################## Local effects ###############################
	
	#Sort, scale and calculate log2fc
	bedlines = sorted(bedlines, key=lambda line: (line["TFBS_chr"], int(line["TFBS_start"]), int(line["TFBS_end"])))
	for line in bedlines:
	
		#Condition specific
		for condition in args.cond_names:
			threshold = args.thresholds[condition]
			line[condition + "_score"] = float(line[condition + "_score"])
			line[condition + "_score"] = np.round(args.norm_objects[condition].normalize(line[condition + "_score"] ), 5)
			line[condition + "_bound"] = 1 if line[condition + "_score"] > threshold else 0

		#Comparison specific
		for i, (cond1, cond2) in enumerate(comparisons):
			base = "{0}_{1}".format(cond1, cond2)
			line[base + "_log2fc"] = round(np.log2((line[cond1 + "_score"] + args.pseudo) / (line[cond2 + "_score"] + args.pseudo)), 5)

	#### Write _all file ####
	outfile = os.path.join(bed_outdir, TF_name + "_all.bed")
	dict_to_tab(bedlines, outfile, header)

	#### Write _bound/_unbound files ####
	for condition in args.cond_names:
		chosen_columns = header[:-no_cond-1] + [condition + "_score"]

		#Subset bedlines per state
		for state in ["bound", "unbound"]:
			outfile = os.path.join(bed_outdir, "{0}_{1}_{2}.bed".format(TF_name, condition, state))
			chosen_bool = 1 if state == "bound" else 0
			bedlines_subset = [bedline for bedline in bedlines if bedline[condition + "_bound"] == chosen_bool]
			bedlines_subset = sorted(bedlines_subset, key= lambda line: line[condition + "_score"], reverse=True)
			dict_to_tab(bedlines_subset, outfile, chosen_columns)

	##### Write overview with scores, bound and log2fcs ####
	overview_columns = header + [condition + "_bound" for condition in args.cond_names] + ["{0}_{1}_log2fc".format(cond1, cond2) for (cond1, cond2) in comparisons]
	overview_txt = os.path.join(args.outdir, TF_name, TF_name + "_overview.txt")
	dict_to_tab(bedlines, overview_txt, overview_columns, header=True)
	
	#Write xlsx overview
	bed_table = pd.DataFrame(bedlines)
	try:
		overview_excel = os.path.join(args.outdir, TF_name, TF_name + "_overview.xlsx")
		writer = pd.ExcelWriter(overview_excel, engine='xlsxwriter')
		bed_table.to_excel(writer, index=False, columns=overview_columns)
		
		worksheet = writer.sheets['Sheet1']
		no_rows, no_cols = bed_table.shape
		worksheet.autofilter(0,0,no_rows, no_cols-1)
		writer.save()

	except:
		print("Error writing excelfile for TF {0}".format(TF_name))
		sys.exit()


	############################## Global effects ##############################

	#Get info table ready
	info_columns = ["total_tfbs"]
	info_columns.extend(["{0}_{1}".format(cond, metric) for (cond, metric) in itertools.product(args.cond_names, ["mean_score", "bound"])])
	info_columns.extend(["{0}_{1}_{2}".format(comparison[0], comparison[1], metric) for (comparison, metric) in itertools.product(comparisons, ["change", "pvalue"])])
	rows, cols = 1, len(info_columns)
	info_table = pd.DataFrame(np.zeros((rows, cols)), columns=info_columns, index=[TF_name])

	#Fill in info table
	info_table.at[TF_name, "total_tfbs"] = n_rows

	for condition in args.cond_names:
		info_table.at[TF_name, condition + "_mean_score"] = round(np.mean(bed_table[condition + "_score"]), 5)
		info_table.at[TF_name, condition + "_bound"] = np.sum(bed_table[condition + "_bound"].values) #_bound contains bool 0/1
		
	#### Calculate statistical test for binding in comparison to background ####
	fig_out = os.path.abspath(os.path.join(args.outdir, TF_name, "plots", TF_name + "_log2fcs.pdf"))
	log2fc_pdf = PdfPages(fig_out, keep_empty=True)

	for i, (cond1, cond2) in enumerate(comparisons):
		base = "{0}_{1}".format(cond1, cond2)

		# Compare log2fcs to background log2fcs
		included = np.logical_or(bed_table[cond1 + "_score"].values > 0, bed_table[cond2 + "_score"].values > 0)
		subset = bed_table[included].copy() 		#included subset 
		subset.loc[:,"peak_id"] = ["_".join([chrom, str(start), str(end)]) for (chrom, start, end) in zip(subset["peak_chr"].values, subset["peak_start"].values, subset["peak_end"].values)]	
		
		observed_log2fcs = subset.groupby('peak_id')[base + '_log2fc'].mean().reset_index()[base + "_log2fc"].values		#if more than one TFBS per peak -> take mean value

		#Estimate mean/std
		bg_params = log2fc_params[(cond1, cond2)]
		obs_params = scipy.stats.norm.fit(observed_log2fcs)

		obs_mean, obs_std = obs_params
		bg_mean, bg_std = bg_params
		obs_no = np.min([len(observed_log2fcs), 50000])		#Set cap on obs_no to prevent super small p-values

		#If there was any change found at all (0 can happen if two bigwigs are the same)
		if obs_mean != bg_mean: 
			info_table.at[TF_name, base + "_change"] = (obs_mean - bg_mean) / np.mean([obs_std, bg_std])  #effect size
			info_table.at[TF_name, base + "_change"] = np.round(info_table.at[TF_name, base + "_change"], 5)
		
			#pval = scipy.stats.mannwhitneyu(observed_log2fcs, bg_log2fcs, alternative="two-sided")[1]
			pval = scipy.stats.ttest_ind_from_stats(obs_mean, obs_std, obs_no, bg_mean, bg_std, obs_no, equal_var=False)[1] 	#pvalue is second in tup
			info_table.at[TF_name, base + "_pvalue"] = pval
		
		#Else not possible to compare groups
		else:
			info_table.at[TF_name, base + "_change"] = 0
			info_table.at[TF_name, base + "_pvalue"] = 1

		#### Plot comparison ###
		fig, ax = plt.subplots(1,1)
		ax.hist(observed_log2fcs, bins='auto', label="Observed log2fcs", density=True)
		xvals = np.linspace(plt.xlim()[0], plt.xlim()[1], 100)
		
		#Observed distribution
		pdf = scipy.stats.norm.pdf(xvals, *obs_params)
		ax.plot(xvals, pdf, label="Observed distribution (fit)", color="red", linestyle="--")
		ax.axvline(obs_mean, color="red", label="Observed mean")
		
		#Background distribution
		pdf = scipy.stats.norm.pdf(xvals, *bg_params)
		ax.plot(xvals, pdf, label="Background distribution (fit)", color="Black", linestyle="--")
		ax.axvline(bg_mean, color="black", label="Background mean")

		#Set size
		x0,x1 = ax.get_xlim()
		y0,y1 = ax.get_ylim()
		ax.set_aspect(((x1-x0)/(y1-y0)) / 1.5)

		#Decorate
		ax.legend()
		plt.xlabel("Log2 fold change", fontsize=8)
		plt.ylabel("Density", fontsize=8)
		plt.title("Differential binding for TF \"{0}\"\nbetween ({1} / {2})".format(TF_name, cond1, cond2), fontsize=10)
		ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
		
		plt.tight_layout()
		log2fc_pdf.savefig(fig, bbox_inches='tight')
		plt.close(fig)

	log2fc_pdf.close()	
	
	#################### Remove temporary file ######################
	try:
		os.remove(filename)
	except:
		logger.error("Could not remove temporary file {0} - this does not effect the results of BINDetect.".format(filename) )

	return(info_table)


#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#

def plot_bindetect(motifs, cluster_obj, conditions, args):
	"""
		Conditions refer to the order of the fold_change divison, meaning condition1/condition2 
		- Clusters is a RegionCluster object 
		- conditions is a tup of condition names (cond1, cond2)
	"""
	warnings.filterwarnings("ignore")

	cond1, cond2 = conditions
	no_IDS = cluster_obj.n

	#Link information from motifs / clusters
	diff_scores = {}
	for motif in motifs:
		diff_scores[motif.prefix] = {"change": motif.change,
									"pvalue": motif.pvalue,
									"log10pvalue": -np.log10(motif.pvalue) if  motif.pvalue > 0 else -np.log10(1e-308),	#smallest possible number before python underflows
									"volcano_label": motif.name,	#shorter name
									"overview_label": "{0} ({1})".format(motif.name, motif.id) 		#the name which was output used in bindetect output
									}
	
	xvalues = np.array([diff_scores[TF]["change"] for TF in diff_scores])
	yvalues = np.array([diff_scores[TF]["log10pvalue"] for TF in diff_scores])

	#### Define the TFs to plot IDs for ####
	y_min = np.percentile(yvalues[yvalues < -np.log10(1e-300)], 95)	
	x_min, x_max = np.percentile(xvalues, [5,95])

	for TF in diff_scores:
		if diff_scores[TF]["change"] < x_min or diff_scores[TF]["change"] > x_max or diff_scores[TF]["log10pvalue"] > y_min:
			diff_scores[TF]["show"] = True
			if diff_scores[TF]["change"] < 0:
				diff_scores[TF]["color"] = "blue"
			elif diff_scores[TF]["change"] > 0:
				diff_scores[TF]["color"] = "red"
		else:
			diff_scores[TF]["show"] = False 
			diff_scores[TF]["color"] = "black"

	node_color = cluster_obj.node_color
	IDS = np.array(cluster_obj.names)
	
	"""
	#Set cluster names
	for motif_name in diff_scores:
		for cluster in cluster_obj.clusters:

			if motif_name in cluster_obj.clusters[cluster]["member_names"]:
				diff_scores[motif_name]["cluster_name"] = cluster_obj.clusters[cluster]["cluster_name"]

			if motif_name == cluster_obj.clusters[cluster]["representative"]:
				diff_scores[TF]["show"] = True
				diff_scores[motif_name]["representative"] = True
	"""

	#--------------------------------------- Figure --------------------------------#

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

	#Decorate plot
	ax1.set_xlabel("Differential binding score")
	ax1.set_ylabel("-log10(pvalue)")

	########### Dendrogram over similarities of TFs #######
	
	dendro_dat = dendrogram(cluster_obj.linkage_mat, labels=IDS, no_labels=True, orientation="right", ax=ax3, above_threshold_color="black", link_color_func=lambda k: cluster_obj.node_color[k])
	labels = dendro_dat["ivl"]	#Now sorted for the order in dendrogram
	ax3.set_xlabel("Transcription factor similarities\n(Clusters below threshold are colored)")

	ax3.set_ylabel("Transcription factor clustering based on TFBS overlap", rotation=270, labelpad=20)
	ax3.yaxis.set_label_position("right")

	#Set aspect of dendrogram/changes
	x0,x1 = ax3.get_xlim()
	y0,y1 = ax3.get_ylim()
	ax3.set_aspect(((x1-x0)/(y1-y0)) * no_IDS/10)		

	########## Differential binding scores per TF ##########
	ax2.set_xlabel("Differential binding score\n" + "(" + cond2 + r' $\leftarrow$' + r'$\rightarrow$ ' + cond1 + ")") #First position in comparison equals numerator in log2fc division
	ax2.xaxis.set_label_position('bottom') 
	ax2.xaxis.set_ticks_position('bottom') 

	no_labels = len(labels)
	ax2.set_ylim(0.5, no_labels+0.5)
	ax2.set_ylabel("Transcription factors")

	ax2.set_yticks(range(1,no_labels+1))
	ax2.set_yticklabels([diff_scores[TF]["overview_label"] for TF in labels])
	ax2.axvline(0, color="grey", linestyle="--") 	#Plot line at middle

	#Plot scores per TF
	for y, TF in enumerate(labels):	#labels are the output motif names from output
		

		idx = np.where(IDS == TF)[0][0]
		score = diff_scores[TF]["change"]

		#Set coloring based on change/pvalue
		if diff_scores[TF]["show"] == True:
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

	######### Color points and set labels in volcano ########
	txts = []
	for TF in diff_scores:
		coord = [diff_scores[TF]["change"], diff_scores[TF]["log10pvalue"]]
		ax1.scatter(coord[0], coord[1], color=diff_scores[TF]["color"], s=4.5)

		if diff_scores[TF]["show"] == True:
			txts.append(ax1.text(coord[0], coord[1], diff_scores[TF]["volcano_label"], fontsize=7))

	adjust_text(txts, ax=ax1, text_from_points=True, arrowprops=dict(arrowstyle='-', color='black', lw=0.5))  #, expand_text=(0.1,1.2), expand_objects=(0.1,0.1))
	
	"""
	#Add arrows to other cluster members
	print(txts[0].__dict__)
	label_positions = {text._text:text for text in txts}
	print(label_positions)
	for TF in diff_scores:
		if diff_scores[TF]["show"]:
			cluster_name = diff_scores[TF]["cluster_name"]
			
			if cluster_name in label_positions: 
				print(cluster_name)

				point_x, point_y = diff_scores[TF]["change"], diff_scores[TF]["log10pvalue"]
				text_x, text_y = label_positions[cluster_name]._x, label_positions[cluster_name]._y
				len_x, len_y = text_x - point_x, text_y - point_y

				ax1.arrow(point_x, point_y, len_x, len_y, linestyle="-", color="black", lw=0.5)
	"""
	#print(txts)

	#Plot custom legend for colors
	legend_elements = [Line2D([0],[0], marker='o', color='w', markerfacecolor="red", label="Higher scores in {0}".format(conditions[0])),
						Line2D([0],[0], marker='o', color='w', markerfacecolor="blue", label="Higher scores in {0}".format(conditions[1]))]
	ax1.legend(handles=legend_elements, bbox_to_anchor=(1.05, 1), loc='upper left')

	return(fig)

"""
def multi_annotate(ax, s, xy_arr=[], *args, **kwargs):
	ans = []
	an = ax.annotate(s, xy_arr[0], *args, **kwargs)
	ans.append(an)
	d = {}
	try:
		d['xycoords'] = kwargs['xycoords']
	except KeyError:
    pass
	  try:
    d['arrowprops'] = kwargs['arrowprops']
  except KeyError:
    pass
  for xy in xy_arr[1:]:
    an = ax.annotate(s, xy, alpha=0.0, xytext=(0,0), textcoords=an, **d)
    ans.append(an)
  return ans
"""