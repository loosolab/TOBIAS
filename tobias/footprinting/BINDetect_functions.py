import numpy as np
from sklearn import mixture
#import matplotlib.pyplot as plt
#from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
import scipy
from datetime import datetime
import itertools
import matplotlib.pyplot as plt
import xlsxwriter
import random

import logging
import logging.handlers
from decimal import Decimal

#Bio-specific packages
import pyBigWig
import pysam
import MOODS.scan
import MOODS.tools
import MOODS.parsers

from tobias.utils.regions import *
from tobias.utils.sequences import *
from tobias.utils.utilities import *


#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#

def get_gc_content(regions, fasta):
	""" Get GC content from regions in fasta """
	gc = 0
	total = 0
	fasta_obj = pysam.FastaFile(fasta)
	for region in regions:
		seq = fasta_obj.fetch(region.chrom, region.start, region.end)
		gc += sum([1 for nuc in seq if (nuc == "G" or nuc == "g" or nuc == "C" or nuc == "c")])
		total += region.end - region.start
	fasta_obj.close()
	gc_content = gc / float(total)

	return(gc_content)


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

	background_signal = {condition:[] for condition in args.cond_names}
	#TODO: Estimate number of background positions sampled to pre-allocate space

	######## Scan for motifs in each region ######
	all_TFBS = {TF: RegionList() for TF in motifs_obj.names} 	# Dict for saving sites before writing
	for i, region in enumerate(regions):
		logger.debug("Processing region: {0}".format(region.tup()))
	
		extra_columns = region

		#Random positions for sampling
		reglen = region.get_length()
		rand_positions = random.sample(range(reglen), max(1,int(reglen/200.0)))		#theoretically one in every 200 bp
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
				background_signal[condition].append(footprints[condition][pos])

		#Scan for motifs across sequence from fasta
		seq = fasta_obj.fetch(region.chrom, region.start, region.end)
		region_TFBS = motifs_obj.scan_sequence(seq, region)		#RegionList of TFBS

		#Extend all TFBS with extra columns from peaks and bigwigs 
		extra_columns = region
		for TFBS in region_TFBS:
			TFBS.extend(extra_columns)

			#Assign scores from bigwig
			for bigwig in args.cond_names:
				motif_length = TFBS.end - TFBS.start 
				bigwig_score = footprints[bigwig][TFBS.start - region.start + int(motif_length/2.0)]	#middle of TFBS
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
	filename = os.path.join(bed_outdir, TF_name + "_all.bed")
	no_cond = len(args.cond_names)
	comparisons = list(itertools.combinations(args.cond_names, 2))

	#Get overview file ready
	#comparisons = list(itertools.combinations(args.cond_names, 2))

	#Read file to pandas
	arr = np.genfromtxt(filename, dtype=None, delimiter="\t", names=None, encoding="utf8", autostrip=True)	#Read using genfromtxt to get automatic type
	bed_table = pd.DataFrame(arr, index=None, columns=None)
	no_rows, no_cols = bed_table.shape

	#Set header in overview
	header = [""]*no_cols
	header[:6] = ["TFBS_chr", "TFBS_start", "TFBS_end", "TFBS_name", "TFBS_score", "TFBS_strand"]
	
	if args.peak_header_list != None:
		header[6:6+len(args.peak_header_list)] = args.peak_header_list
	else:
		no_peak_col = len(header[6:])
		header[6:6+no_peak_col] = ["peak_chr", "peak_start", "peak_end"] + ["additional_" + str(num + 1) for num in range(no_peak_col-3)]

	header[-no_cond:] = ["{0}_score".format(condition) for condition in args.cond_names] 	#signal scores
	bed_table.columns = header
	
	#Sort and format
	bed_table = bed_table.sort_values(["TFBS_chr", "TFBS_start", "TFBS_end"])
	for condition in args.cond_names:
		bed_table[condition + "_score"] = bed_table[condition + "_score"].round(5)

	#### Estimate bound/unbound split ####
	for condition in args.cond_names:

		threshold = args.thresholds[condition]
		bed_table[condition + "_bound"] = np.where(bed_table[condition + "_score"] > threshold, 1, 0).astype(int)

	#Write bound/unbound
	for (condition, state) in itertools.product(args.cond_names, ["bound", "unbound"]):

		outfile = os.path.join(bed_outdir, "{0}_{1}_{2}.bed".format(TF_name, condition, state))

		#Subset bed table
		chosen_bool = 1 if state == "bound" else 0
		bed_table_subset = bed_table.loc[bed_table[condition + "_bound"] == chosen_bool]
		
		#Write out subset with subset of columns
		output_columns = header[:-no_cond] + [condition + "_score"]		#names of columns
		bed_table_subset = bed_table_subset.loc[:,output_columns]
		col = list(bed_table_subset)[-1]
		bed_table_subset = bed_table_subset.sort_values(col, ascending=False)
		bed_table_subset.to_csv(outfile, sep="\t", index=False, header=False)

	#### Calculate statistical test in comparison to condition comparison ####
	for i, (cond1, cond2) in enumerate(comparisons):
		base = "{0}_{1}".format(cond1, cond2)

		cond1_values = bed_table[cond1 + "_score"].values
		cond2_values = bed_table[cond2 + "_score"].values
		bed_table[base + "_log2fc"] = np.log2(np.true_divide(cond1_values + args.pseudo, cond2_values + args.pseudo))
		bed_table[base + "_log2fc"] = bed_table[base + "_log2fc"].round(5)

		"""
		#If more than one TFBS per peak -> take mean value
		df = bed_table.copy(deep=True)
		df['peak_id'] = list(zip(df['peak_chr'], df['peak_start'], df['peak_end']))
		df = df.groupby('peak_id')[base + '_log2fc'].mean().reset_index()

		##### Compare log2fcs to background #####
		observed_log2fcs = df[base + "_log2fc"]
		observed_log2fcs = observed_log2fcs[np.logical_not(np.isclose(observed_log2fcs,0))] 	#Remove 0 log2fcs from observed
		obs_mean, obs_std, obs_no = (np.mean(observed_log2fcs), np.std(observed_log2fcs), len(observed_log2fcs))

		bg_mean, bg_std, bg_no = log2fc_params[(cond1, cond2)]

		info_table.at[TF_name, base + "_change"] = (obs_mean - bg_mean) / ((bg_std + obs_std)*0.5)		#effect size
		info_table.at[TF_name, base + "_change"] = np.round(info_table.at[TF_name, base + "_change"], 5)
	
		pval = scipy.stats.ttest_ind_from_stats(obs_mean, obs_std, obs_no, bg_mean, bg_std, bg_no, equal_var=False)[1] 	#pvalue is second in tup
		"""
		#bed_table[base + "_log2fc"] = bed_table[base + "_log2fc"].round(5)
		
	#Write overview with scores, bound and log2fcs
	overview_txt = os.path.join(args.outdir, TF_name, TF_name + "_overview.txt")
	bed_table.to_csv(overview_txt, sep="\t", index=False, header=True)

	#Write xlsx overview
	try:
		overview_excel = os.path.join(args.outdir, TF_name, TF_name + "_overview.xlsx")
		writer = pd.ExcelWriter(overview_excel, engine='xlsxwriter')
		bed_table.to_excel(writer, index=False)
		
		worksheet = writer.sheets['Sheet1']
		no_rows, no_cols = bed_table.shape
		worksheet.autofilter(0,0,no_rows, no_cols-1)
		writer.save()

	except:
		print("Error writing excelfile for TF {0}".format(TF_name))
		#logger.critical("Error writing excelfile for TF {0}".format(TF_name)
		
	return(1)



def calculate_global_stats(TF_name, args, log2fc_params):
	""" Reads overview file for each TF and calculates differential score and overall stats """

	filename = os.path.join(args.outdir, TF_name, TF_name + "_overview.txt")
	no_cond = len(args.cond_names)
	conditions = args.cond_names
	comparisons = list(itertools.combinations(conditions, 2))

	#Get info table ready
	info_columns = ["total_tfbs"]
	info_columns.extend(["{0}_{1}".format(cond, metric) for (cond, metric) in itertools.product(args.cond_names, ["bound"])])
	info_columns.extend(["{0}_{1}_{2}".format(comparison[0], comparison[1], metric) for (comparison, metric) in itertools.product(comparisons, ["change", "pvalue"])])

	rows, cols = 1, len(info_columns)
	info_table = pd.DataFrame(np.zeros((rows, cols)), columns=info_columns, index=[TF_name])

	#Read overview file to pandas
	#arr = np.genfromtxt(filename, dtype=None, delimiter="\t", encoding="utf8", autostrip=True)	#Read using genfromtxt to get automatic type
	#print(arr)
	overview_table = pd.read_csv(filename, delimiter="\t") #, DataFrame(arr, index=None, columns=None)

	#Global stats
	no_rows, no_cols = overview_table.shape
	info_table.at[TF_name, "total_tfbs"] = no_rows

	#Number bound per TF 
	for condition in conditions:
		info_table.at[TF_name, condition + "_bound"] = np.sum(overview_table[condition + "_bound"].values)	#_bound contains bool 0/1

	##### Calculate global differential binding #####
	overview_table['peak_id'] = list(zip(overview_table['peak_chr'], overview_table['peak_start'], overview_table['peak_end']))		#used for resolcing multiple sites per peak
	for i, (cond1, cond2) in enumerate(comparisons):
		base = "{0}_{1}".format(cond1, cond2)

		# Compare log2fcs to background 
		observed_log2fcs = overview_table.groupby('peak_id')[base + '_log2fc'].mean().reset_index()[base + "_log2fc"].values	#if more than one TFBS per peak -> take mean value
		observed_log2fcs = observed_log2fcs[np.logical_not(np.isclose(observed_log2fcs,0))] 	#Remove 0 log2fcs from observed
		
		obs_mean, obs_std, obs_no = (np.mean(observed_log2fcs), np.std(observed_log2fcs), len(observed_log2fcs))
		bg_mean, bg_std, bg_no = log2fc_params[(cond1, cond2)]
		info_table.at[TF_name, base + "_change"] = (obs_mean - bg_mean) / ((bg_std + obs_std)*0.5)		#effect size
		info_table.at[TF_name, base + "_change"] = np.round(info_table.at[TF_name, base + "_change"], 5)

		pval = scipy.stats.ttest_ind_from_stats(obs_mean, obs_std, obs_no, bg_mean, bg_std, bg_no, equal_var=False)[1] 	#pvalue is second in tup
		info_table.at[TF_name, base + "_pvalue"] = pval

	return(info_table)