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
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

#Bio-specific packages
import pyBigWig
import pysam
import MOODS.scan
import MOODS.tools
import MOODS.parsers

from tobias.utils.regions import *
from tobias.utils.sequences import *
from tobias.utils.utilities import *
from tobias.utils.motifs import *
from tobias.utils.signals import *
from tobias.plotting.plot_bindetect import *

#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#

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

def find_nearest_idx(array, value):
    idx = np.abs(array - value).argmin()
    return(idx)

def is_nan(x):
    return (x is np.nan or x != x)

### Logger runs in main process
def main_logger_process(queue, logger):

	logger.debug("Started main logger process")
	while True:
		try:
			record = queue.get()
			if record is None:
				break
			logger.handle(record) 

		except Exception:
			import sys, traceback
			print('Problem in main logger process:', file=sys.stderr)
			traceback.print_exc(file=sys.stderr)
			break

	return(1)

"""
def norm_fit(x, loc, scale, size):
	size scales the height of the pdf-curve 
	return(scipy.stats.norm.pdf(x, loc, scale) * size)
	

def lognorm_fit(x, s, loc, scale, size):
	size scales the height of the pdf-curve 

	y = scipy.stats.lognorm.pdf(x, s, loc=loc, scale=scale) * size
	return(y) 
"""

#---------------------------------------------------------------------------------------------------------#
#------------------------------------------- Main functions ----------------------------------------------#
#---------------------------------------------------------------------------------------------------------#

def scan_and_score(regions, motifs_obj, args, log_q, qs):
	""" Scanning and scoring runs in parallel for subsets of regions """
	
	logger = create_mp_logger(args.verbosity, log_q)	#sending all logger calls to log_q

	logger.debug("Setting scanner")
	motifs_obj.setup_moods_scanner()	#MotifList object

	logger.debug("Opening bigwigs and fasta")
	pybw = {condition:pyBigWig.open(args.signals[i], "rb") for i, condition in enumerate(args.cond_names)}
	fasta_obj = pysam.FastaFile(args.genome)
	chrom_boundaries = dict(zip(fasta_obj.references, fasta_obj.lengths))

	gc_window = 500
	rand_window = 500
	extend = int(np.ceil(gc_window / 2.0))

	background_signal = {"gc":[], "signal":{condition:[] for condition in args.cond_names}}

	#TODO: Estimate number of background positions sampled to pre-allocate space

	######## Scan for motifs in each region ######
	all_TFBS = {TF: RegionList() for TF in motifs_obj.names} 	# Dict for saving sites before writing
	for i, region in enumerate(regions):
		logger.debug("Processing region: {0}".format(region.tup()))
	
		extra_columns = region

		#Random positions for sampling
		reglen = region.get_length()
		random.seed(reglen)		#Each region is processed identifically regardless of order in file
		rand_positions = random.sample(range(reglen), max(1,int(reglen/rand_window)))		#theoretically one in every 500 bp
		logger.debug("Random indices: {0} for region length {1}".format(rand_positions, reglen))

		#Read footprints in region
		footprints = {}
		for condition in args.cond_names:
			try:
				footprints[condition] = pybw[condition].values(region.chrom, region.start, region.end, numpy=True)
				footprints[condition] = np.nan_to_num(footprints[condition])	#nan to 0
			except:
				logger.critical("ERROR reading footprints from region: {0}".format(region))
				continue

			#Read random positions for background
			for pos in rand_positions:
				background_signal["signal"][condition].append(footprints[condition][pos])

		#Scan for motifs across sequence from fasta
		extended_region = copy.copy(region).extend_reg(extend)	 #extend to calculate gc
		extended_region.check_boundary(chrom_boundaries, action="cut")

		seq = fasta_obj.fetch(region.chrom, extended_region.start, extended_region.end)

		#Calculate GC content for regions
		num_sequence = nuc_to_num(seq) 
		Ns = num_sequence == 4
		boolean = 1 * (num_sequence > 1)		# Convert to 0/1 gc
		boolean[Ns] = 0.5						# replace Ns 0.5 - with neither GC nor AT
		boolean = boolean.astype(np.float64)	# due to input of fast_rolling_math
		gc = fast_rolling_math(boolean, gc_window, "mean")
		gc = gc[extend:-extend]
		background_signal["gc"].extend([gc[pos] for pos in rand_positions])

		region_TFBS = motifs_obj.scan_sequence(seq[extend:-extend], region)		#RegionList of TFBS

		#Extend all TFBS with extra columns from peaks and bigwigs 
		extra_columns = region
		for TFBS in region_TFBS:
			motif_length = TFBS.end - TFBS.start 
			pos = TFBS.start - region.start + int(motif_length/2.0) #middle of TFBS
			
			TFBS.extend(extra_columns)
			TFBS.append(gc[pos])

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

		if name in args.new_motifs: 	# Only write if name was in the requested motifs
			logger.debug("Sending {0} sites from {1} to bed-writer queue".format(no_sites, name))
			bed_content = all_TFBS[name].as_bed()	#string 
			qs[name].put((name, bed_content))

		global_TFBS.extend(all_TFBS[name])
		all_TFBS[name] = []

	overlap = global_TFBS.count_overlaps()

	#Close down open file handles
	fasta_obj.close()
	for bigwig_f in pybw:
		pybw[bigwig_f].close()
			
	end_time = datetime.now()

	return(background_signal, overlap)


#-----------------------------------------------------------------------------------------------#
def process_tfbs(TF_name, args, log2fc_params): 	#per tf
	""" Processes single TFBS to split into bound/unbound and write out overview file """

	begin_time = datetime.now()

	bed_outdir = os.path.join(args.outdir, TF_name, "beds")
	filename = os.path.join(bed_outdir, TF_name + ".tmp")
	no_cond = len(args.cond_names)
	comparisons = list(itertools.combinations(args.cond_names, 2))

	#Get info table ready
	info_columns = ["total_tfbs"]
	info_columns.extend(["{0}_{1}".format(cond, metric) for (cond, metric) in itertools.product(args.cond_names, ["bound"])])
	info_columns.extend(["{0}_{1}_{2}".format(comparison[0], comparison[1], metric) for (comparison, metric) in itertools.product(comparisons, ["change", "pvalue"])])
	rows, cols = 1, len(info_columns)
	info_table = pd.DataFrame(np.zeros((rows, cols)), columns=info_columns, index=[TF_name])

	#Read file to pandas
	arr = np.genfromtxt(filename, dtype=None, delimiter="\t", names=None, encoding="utf8", autostrip=True)	#Read using genfromtxt to get automatic type
	bed_table = pd.DataFrame(arr, index=None, columns=None)
	no_rows, no_cols = bed_table.shape

	#no_rows, no_cols = overview_table.shape
	info_table.at[TF_name, "total_tfbs"] = no_rows

	#Set header in overview
	header = [""]*no_cols
	header[:6] = ["TFBS_chr", "TFBS_start", "TFBS_end", "TFBS_name", "TFBS_score", "TFBS_strand"]
	
	if args.peak_header_list != None:
		header[6:6+len(args.peak_header_list)] = args.peak_header_list
	else:
		no_peak_col = len(header[6:])
		header[6:6+no_peak_col] = ["peak_chr", "peak_start", "peak_end"] + ["additional_" + str(num + 1) for num in range(no_peak_col-3)]

	header[-no_cond:] = ["{0}_score".format(condition) for condition in args.cond_names] 	#signal scores
	header[-no_cond-1] = "GC"
	bed_table.columns = header

	#Sort and format
	bed_table = bed_table.sort_values(["TFBS_chr", "TFBS_start", "TFBS_end"])
	for condition in args.cond_names:
		bed_table[condition + "_score"] = bed_table[condition + "_score"].round(5)

	#### Write all file ####
	chosen_columns = [col for col in header if col != "GC"]
	outfile = os.path.join(bed_outdir, TF_name + "_all.bed")
	bed_table.to_csv(outfile, sep="\t", index=False, header=False, columns=chosen_columns)

	#### Estimate bound/unbound split ####
	for condition in args.cond_names:

		threshold = args.thresholds[condition]
		bed_table[condition + "_bound"] = np.where(bed_table[condition + "_score"] > threshold, 1, 0).astype(int)
		info_table.at[TF_name, condition + "_bound"] = np.sum(bed_table[condition + "_bound"].values)	#_bound contains bool 0/1

	#Write bound/unbound
	for (condition, state) in itertools.product(args.cond_names, ["bound", "unbound"]):

		outfile = os.path.join(bed_outdir, "{0}_{1}_{2}.bed".format(TF_name, condition, state))

		#Subset bed table
		chosen_bool = 1 if state == "bound" else 0
		bed_table_subset = bed_table.loc[bed_table[condition + "_bound"] == chosen_bool]
		bed_table_subset = bed_table_subset.sort_values([condition + "_score"], ascending=False)

		#Write out subset with subset of columns
		chosen_columns = header[:-no_cond-1] + [condition + "_score"]
		bed_table_subset.to_csv(outfile, sep="\t", index=False, header=False, columns=chosen_columns)

	#### Calculate statistical test in comparison to background ####
	fig_out = os.path.abspath(os.path.join(args.outdir, TF_name, "plots", TF_name + "_log2fcs.pdf"))
	log2fc_pdf = PdfPages(fig_out, keep_empty=True)

	for i, (cond1, cond2) in enumerate(comparisons):
		base = "{0}_{1}".format(cond1, cond2)

		#Calculate log2fcs of TFBS for this TF
		cond1_values = bed_table[cond1 + "_score"].values
		cond2_values = bed_table[cond2 + "_score"].values
		bed_table[base + "_log2fc"] = np.log2(np.true_divide(cond1_values + args.pseudo, cond2_values + args.pseudo))
		bed_table[base + "_log2fc"] = bed_table[base + "_log2fc"].round(5)
		
		# Compare log2fcs to background log2fcs
		excluded = np.logical_and(np.isclose(bed_table[cond1 + "_score"].values, 0), np.isclose(bed_table[cond2 + "_score"].values, 0))		#exclude 0's from both conditions
		subset = bed_table[np.logical_not(excluded)].copy() 		#included subset 
		subset.loc[:,"peak_id"] = ["_".join([chrom, str(start), str(end)]) for (chrom, start, end) in zip(subset["peak_chr"].values, subset["peak_start"].values, subset["peak_end"].values)]	
		observed_log2fcs = subset.groupby('peak_id')[base + '_log2fc'].mean().reset_index()[base + "_log2fc"].values		#if more than one TFBS per peak -> take mean value
		observed_gcs = subset.groupby('peak_id')["GC"].mean().reset_index()["GC"].values									#if more than one TFBS per peak -> take mean value

		#Resample from background
		log2fc_params[(cond1, cond2)].set_params(random_state=len(observed_log2fcs))
		sampling, labels = log2fc_params[(cond1, cond2)].sample(int(info_table.at[TF_name, "total_tfbs"]*2))
		sampled_log2fcs, sampled_gcs = sampling[:,0], sampling[:,1]

		indices = []
		for val in observed_gcs:
			idx = find_nearest_idx(sampled_gcs, val)
			indices.append(idx)
			sampled_gcs[idx] = np.inf

		bg_log2fcs = sampled_log2fcs[indices]

		#Estimate mean/std
		bg_params = scipy.stats.norm.fit(bg_log2fcs)
		obs_params = scipy.stats.norm.fit(observed_log2fcs)

		obs_mean, obs_std = obs_params
		bg_mean, bg_std = bg_params
		obs_no = np.min([len(observed_log2fcs), 50000])		#Set cap on obs_no to prevent super small p-values

		#If there was any change found at all (0 can happen if two bigwigs are the same)
		if obs_mean != bg_mean: 
			info_table.at[TF_name, base + "_change"] = (obs_mean - bg_mean) / bg_std  #effect size
			info_table.at[TF_name, base + "_change"] = np.round(info_table.at[TF_name, base + "_change"], 5)
		
			#pval = scipy.stats.mannwhitneyu(observed_log2fcs, bg_log2fcs, alternative="two-sided")[1]
			pval = scipy.stats.ttest_ind_from_stats(obs_mean, obs_std, obs_no, bg_mean, bg_std, obs_no, equal_var=False)[1] 	#pvalue is second in tup
			info_table.at[TF_name, base + "_pvalue"] = pval
		
		#Else not possible to compare groups
		else:
			info_table.at[TF_name, base + "_change"] = 0
			info_table.at[TF_name, base + "_pvalue"] = 1

		# Plot differences of distributions
		fig, ax = plt.subplots(1, 1)
		ax.hist(observed_log2fcs, density=True, bins=50, color="red", label="Observed log2fcs", alpha=0.5)
		ax.hist(bg_log2fcs, density=True, bins=50, color="black", label="Background log2fcs", alpha=0.6)
		ax.axvline(obs_mean, color="red", label="Observed mean")
		ax.axvline(bg_mean, color="black", label="Background mean")
		
		x0,x1 = ax.get_xlim()
		y0,y1 = ax.get_ylim()
		ax.set_aspect(((x1-x0)/(y1-y0)) / 1.5)		#square volcano plot

		#Add text showing change
		ax.text(0.05, 0.95, "Diff. score: {0:.3f}".format(info_table.at[TF_name, base + "_change"]), va="top", transform = ax.transAxes)

		plt.xlabel("Log2(fold change)")
		plt.ylabel("Density")
		plt.title("Differential binding for TF \"{0}\"\nbetween ({1} / {2})".format(TF_name, cond1, cond2))
		ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

		plt.tight_layout()
		log2fc_pdf.savefig(fig, bbox_inches='tight')
		plt.close(fig)

	log2fc_pdf.close()	


	########## Write overview with scores, bound and log2fcs ##############
	chosen_columns = [col for col in bed_table.columns if col != "GC"]
	overview_txt = os.path.join(args.outdir, TF_name, TF_name + "_overview.txt")
	bed_table.to_csv(overview_txt, sep="\t", index=False, header=True, columns=chosen_columns)

	#Write xlsx overview
	try:
		overview_excel = os.path.join(args.outdir, TF_name, TF_name + "_overview.xlsx")
		writer = pd.ExcelWriter(overview_excel, engine='xlsxwriter')
		bed_table.to_excel(writer, index=False, columns=chosen_columns)
		
		worksheet = writer.sheets['Sheet1']
		no_rows, no_cols = bed_table.shape
		worksheet.autofilter(0,0,no_rows, no_cols-1)
		writer.save()

	except:
		print("Error writing excelfile for TF {0}".format(TF_name))
		sys.exit() #logger.critical("Error writing excelfile for TF {0}".format(TF_name)
	
	#### Remove temporary file ####
	try:
		os.remove(filename)
	except:
		print("Error removing temporary file {0} - this does not effect the results of BINDetect.".format(filename) )

	return(info_table)
