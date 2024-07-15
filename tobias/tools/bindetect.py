#!/usr/bin/env python

"""
BINDetect: Detects differential binding between conditions as well as bound transcription factors from footprints and motifs

@author: Mette Bentsen
@contact: mette.bentsen (at) mpi-bn.mpg.de
@license: MIT
"""

import os
import sys
import argparse
import numpy as np
import multiprocessing as mp
import time
from copy import deepcopy
import logging
import itertools
import pandas as pd
import seaborn as sns
from collections import Counter

#Machine learning and statistics
import sklearn
from sklearn import mixture
import scipy
from kneed import KneeLocator	

#Plotting
import matplotlib
matplotlib.use("Agg")	#non-interactive backend
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import NullFormatter

#Bio-specific packages
import pysam
import pyBigWig as pybw

#Internal functions and classes
from tobias.parsers import add_bindetect_arguments
from tobias.tools.bindetect_functions import *
from tobias.utils.utilities import *
from tobias.utils.regions import *
from tobias.utils.motifs import *
from tobias.utils.logger import TobiasLogger

#For warnings from curve_fit
import warnings
from scipy.optimize import OptimizeWarning
warnings.simplefilter("ignore", OptimizeWarning)
warnings.simplefilter("ignore", RuntimeWarning)


def norm_fit(x, mean, std, scale):
	return(scale * scipy.stats.norm.pdf(x, mean, std))

#----------------------------------------------------------------------------------------------------------------#
def run_bindetect(args):
	""" Main function to run bindetect algorithm with input files and parameters given in args """

	#Checking input and setting cond_names
	check_required(args, ["signals", "motifs", "genome", "peaks"])
	args.cond_names = [os.path.basename(os.path.splitext(bw)[0]) for bw in args.signals] if args.cond_names is None else args.cond_names
	args.outdir = os.path.abspath(args.outdir)

	#Set output files
	states = ["bound", "unbound"]
	outfiles = [os.path.abspath(os.path.join(args.outdir, "*", "beds", "*_{0}_{1}.bed".format(condition, state))) for (condition, state) in itertools.product(args.cond_names, states)]
	outfiles.append(os.path.abspath(os.path.join(args.outdir, "*", "beds", "*_all.bed")))
	outfiles.append(os.path.abspath(os.path.join(args.outdir, "*", "plots", "*_log2fcs.pdf")))
	outfiles.append(os.path.abspath(os.path.join(args.outdir, "*", "*_overview.txt")))
	outfiles.append(os.path.abspath(os.path.join(args.outdir, "*", "*_overview.xlsx")))

	outfiles.append(os.path.abspath(os.path.join(args.outdir, args.prefix + "_distances.txt")))
	outfiles.append(os.path.abspath(os.path.join(args.outdir, args.prefix + "_results.txt")))
	outfiles.append(os.path.abspath(os.path.join(args.outdir, args.prefix + "_results.xlsx")))
	outfiles.append(os.path.abspath(os.path.join(args.outdir, args.prefix + "_figures.pdf")))


	#-------------------------------------------------------------------------------------------------------------#
	#-------------------------------------------- Setup logger and pool ------------------------------------------#
	#-------------------------------------------------------------------------------------------------------------#

	logger = TobiasLogger("BINDetect", args.verbosity)
	logger.begin()

	parser = add_bindetect_arguments(argparse.ArgumentParser())
	logger.arguments_overview(parser, args)
	logger.output_files(outfiles)

	# Setup pool
	args.cores = check_cores(args.cores, logger)
	writer_cores = max(1, int(args.cores*0.1))
	worker_cores = max(1, args.cores - writer_cores)
	logger.debug("Worker cores: {0}".format(worker_cores))
	logger.debug("Writer cores: {0}".format(writer_cores))

	pool = mp.Pool(processes=worker_cores)
	writer_pool = mp.Pool(processes=writer_cores)

	#-------------------------------------------------------------------------------------------------------------#
	#-------------------------- Pre-processing data: Reading motifs, sequences, peaks ----------------------------#
	#-------------------------------------------------------------------------------------------------------------#

	logger.info("----- Processing input data -----")

	#Check that cond_names are the right length and are unique:
	if len(args.cond_names) != len(args.signals):
		logger.error("The given number of given '--cond-names' ({0}) differ from the given input '--signals' ({1}). Please enter one condition name per signal.".format(len(args.cond_names), len(args.signals)))
		sys.exit(1)

	if len(args.cond_names) != len(set(args.cond_names)):
		logger.error("The condition names are not unique ({0}). Please use --cond-names to set a unique set of condition names.".format(args.cond_names))
		sys.exit(1)

	#Check opening/writing of files
	logger.info("Checking reading/writing of files")
	check_files([args.signals, args.motifs, args.genome, args.peaks], action="r")
	check_files(outfiles[-3:], action="w")
	make_directory(args.outdir)

	#Comparisons between conditions
	no_conditions = len(args.signals)
	if args.time_series:
		comparisons = list(zip(args.cond_names[:-1], args.cond_names[1:]))
		args.comparisons = comparisons
	else:
		comparisons = list(itertools.combinations(args.cond_names, 2))	#all-against-all
		args.comparisons = comparisons

	#Pdf for debug output
	if args.debug: 
		debug_out = os.path.abspath(os.path.join(args.outdir, args.prefix + "_debug.pdf"))
		debug_pdf = PdfPages(debug_out, keep_empty=True)

	#Open figure pdf and write overview
	fig_out = os.path.abspath(os.path.join(args.outdir, args.prefix + "_figures.pdf"))
	figure_pdf = PdfPages(fig_out, keep_empty=True)

	plt.figure()
	plt.axis('off')
	plt.text(0.5,0.8, "BINDETECT FIGURES", ha="center", va="center", fontsize=20)

	#output and order
	titles = []
	titles.append("Raw score distributions")
	if no_conditions > 1 and args.norm_off == False:
		titles.append("Normalized score distributions")
	if args.debug:
		for (cond1, cond2) in comparisons:
			titles.append("Background log2FCs ({0} / {1})".format(cond1, cond2))	

	for (cond1, cond2) in comparisons:
		titles.append("BINDetect plot ({0} / {1})".format(cond1, cond2))

	plt.text(0.1, 0.6, "\n".join(["Page {0}) {1}".format(i+2, titles[i]) for i in range(len(titles))]) + "\n\n", va="top")
	figure_pdf.savefig(bbox_inches='tight')
	plt.close()

	################# Read peaks ################
	#Read peak and peak_header
	logger.info("Reading peaks")
	peaks = RegionList().from_bed(args.peaks)
	logger.info("- Found {0} regions in input peaks".format(len(peaks)))

	#Check number of columns in peaks
	n_cols = len(peaks[0])
	for i, peak in enumerate(peaks):
		if len(peak) != n_cols:
			logger.error("The lines in --peaks have a varying number of columns. Line 1 has {0} columns, but line {1} has {2} columns! Please adjust the format of this file to run TOBIAS BINDetect.".format(n_cols, i+1, len(peak)))
			sys.exit(1)

	#Merge overlapping peaks
	peaks = peaks.merge()	
	logger.info("- Merged to {0} regions".format(len(peaks)))

	if len(peaks) == 0:
		logger.error("Input --peaks file is empty!")
		sys.exit(1)
		
	#Read header and check match with number of peak columns
	peak_columns = len(peaks[0]) #number of columns
	logger.debug("--peaks have {0} columns".format(peak_columns))
	if args.peak_header != None:
		content = open(args.peak_header, "r").read()
		args.peak_header_list = content.split()
		logger.debug("Peak header: {0}".format(args.peak_header_list))

		#Check whether peak header fits with number of peak columns
		if len(args.peak_header_list) != peak_columns:
			logger.error("Length of --peak_header ({0}) does not fit number of columns in --peaks ({1}).".format(len(args.peak_header_list), peak_columns))
			sys.exit(1)
	else:
		args.peak_header_list = ["peak_chr", "peak_start", "peak_end"] + ["additional_" + str(num + 1) for num in range(peak_columns-3)]
	logger.debug("Peak header list: {0}".format(args.peak_header_list))

	################# Check for match between peaks and fasta/bigwig #################
	logger.info("Checking for match between --peaks and --fasta/--signals boundaries")
	logger.info("- Comparing peaks to {0}".format(args.genome))
	fasta_obj = pysam.FastaFile(args.genome)
	fasta_boundaries = dict(zip(fasta_obj.references, fasta_obj.lengths))
	fasta_obj.close()
	logger.debug("Fasta boundaries: {0}".format(fasta_boundaries))
	peaks = peaks.apply_method(OneRegion.check_boundary, fasta_boundaries, "exit")	#will exit if peaks are outside borders

	#Check boundaries of each bigwig signal individually
	for signal in args.signals:	
		logger.info("- Comparing peaks to {0}".format(signal))
		pybw_obj = pybw.open(signal)
		pybw_header = pybw_obj.chroms()
		pybw_obj.close()
		logger.debug("Signal boundaries: {0}".format(pybw_header))
		peaks = peaks.apply_method(OneRegion.check_boundary, pybw_header, "exit")

	##### GC content for motif scanning ######
	#Make chunks of regions for multiprocessing
	logger.info("Estimating GC content from peak sequences")
	peak_chunks = peaks.chunks(args.split)
	gc_content_pool = pool.starmap(get_gc_content, itertools.product(peak_chunks, [args.genome])) 
	gc_content = np.mean(gc_content_pool)	#fraction
	args.gc = gc_content
	bg = np.array([(1-args.gc)/2.0, args.gc/2.0, args.gc/2.0, (1-args.gc)/2.0])
	logger.info("- GC content estimated at {0:.2f}%".format(gc_content*100))

	################ Get motifs ################
	logger.info("Reading motifs from file") 
	motif_list = MotifList()
	args.motifs = expand_dirs(args.motifs)
	for f in args.motifs:
		try:
			motif_list += MotifList().from_file(f)  #List of OneMotif objects
		except Exception as e:
			logger.error("Error reading motifs from '{0}'. Error message was: {1}".format(f, e))
			sys.exit(1)

	no_pfms = len(motif_list)
	logger.info("- Read {0} motifs".format(no_pfms))

	logger.debug("Getting motifs ready")
	motif_list.bg = bg

	#Set prefixes
	for motif in motif_list:
		motif.set_prefix(args.naming)
		motif.bg = bg

		logger.spam("Getting pssm for motif {0}".format(motif.name))
		motif.get_pssm()
	
	#Check that prefixes are unique regardless of upper/lower case name
	motif_prefixes = [motif.prefix.upper() for motif in motif_list]
	name_count = Counter(motif_prefixes)
	if max(name_count.values()) > 1:

		duplicated = [key for key, value in name_count.items() if value > 1]
		logger.warning("The motif output names (as given by --naming) are not unique.")
		logger.warning("The following names occur more than once: {0}".format(duplicated))
		logger.warning("These motifs will be renamed with '_1', '_2' etc. To prevent this renaming, please make the names of the input --motifs unique")
		
		motif_count = {dup_motif: 1 for dup_motif in duplicated}
		for i, motif in enumerate(motif_list):
			if motif.prefix.upper() in duplicated:
				
				original_name = motif.prefix
				motif.prefix = motif.prefix + "_{0}".format(motif_count[motif.prefix.upper()])	#Add number to make prefix unique
				logger.debug("Renamed motif {0}: {1} -> {2}".format(i+1, original_name, motif.prefix))
				motif_count[original_name.upper()] += 1

	motif_names = [motif.prefix for motif in motif_list]

	#Get threshold for motifs
	logger.debug("Getting match threshold per motif")
	outlist = pool.starmap(OneMotif.get_threshold, itertools.product(motif_list, [args.motif_pvalue])) 

	logger.spam(motif_list)

	motif_list = MotifList(outlist)	
	for motif in motif_list:
		logger.debug("Motif {0}: threshold {1}".format(motif.name, motif.threshold))

	logger.info("Creating folder structure for each TF")
	for TF in motif_names:
		logger.spam("Creating directories for {0}".format(TF))
		make_directory(os.path.join(args.outdir, TF))
		make_directory(os.path.join(args.outdir, TF, "beds"))
		make_directory(os.path.join(args.outdir, TF, "plots"))

	#-------------------------------------------------------------------------------------------------------------#	
	#----------------------------------------- Plot logos for all motifs -----------------------------------------#
	#-------------------------------------------------------------------------------------------------------------#

	logo_filenames = {motif.prefix: os.path.join(args.outdir, motif.prefix, motif.prefix + ".png") for motif in motif_list}

	logger.info("Plotting sequence logos for each motif")
	task_list = [pool.apply_async(OneMotif.logo_to_file, (motif, logo_filenames[motif.prefix], )) for motif in motif_list]
	monitor_progress(task_list, logger)
	results = [task.get() for task in task_list]
	logger.comment("")

	logger.debug("Getting base64 strings per motif")
	for motif in motif_list:
		#motif.get_base() 
		with open(logo_filenames[motif.prefix], "rb") as png:
			motif.base = base64.b64encode(png.read()).decode("utf-8") 

	#-------------------------------------------------------------------------------------------------------------#
	#--------------------- Motif scanning: Find binding sites and match to footprint scores ----------------------#
	#-------------------------------------------------------------------------------------------------------------#

	logger.comment("")
	logger.start_logger_queue()		#start process for listening and handling through the main logger queue
	args.log_q = logger.queue 		#queue for multiprocessing logging
	manager = mp.Manager()
	logger.info("Scanning for motifs and matching to signals...")

	#Create writer queues for bed-file output
	logger.debug("Setting up writer queues")
	qs_list = []
	writer_qs = {}

	#writer_queue = create_writer_queue(key2file, writer_cores)
	#writer_queue.stop()	#wait until all are done

	manager = mp.Manager()
	TF_names_chunks = [motif_names[i::writer_cores] for i in range(writer_cores)]
	writer_tasks = []
	for TF_names_sub in TF_names_chunks:
		logger.debug("Creating writer queue for {0}".format(TF_names_sub))
		files = [os.path.join(args.outdir, TF, "beds", TF + ".tmp") for TF in TF_names_sub]

		q = manager.Queue()
		qs_list.append(q)

		writer_tasks.append(writer_pool.apply_async(file_writer, args=(q, dict(zip(TF_names_sub, files)), args)))	 #, callback = lambda x: finished.append(x) print("Writing time: {0}".format(x)))
		for TF in TF_names_sub:
			writer_qs[TF] = q
	writer_pool.close() #no more jobs applied to writer_pool

	#todo: use run_parallel
	#Start working on data
	if worker_cores == 1:
		logger.debug("Running with cores = 1")
		results = []
		for chunk in peak_chunks:
			results.append(scan_and_score(chunk, motif_list, args, args.log_q, writer_qs))
		
	else: 
		logger.debug("Sending jobs to worker pool")

		task_list = [pool.apply_async(scan_and_score, (chunk, motif_list, args, args.log_q, writer_qs, )) for chunk in peak_chunks]
		monitor_progress(task_list, logger)
		results = [task.get() for task in task_list]
	
	logger.info("Done scanning for TFBS across regions!")
	#logger.stop_logger_queue()	#stop the listening process (wait until all was written)
	
	#--------------------------------------#
	logger.info("Waiting for bedfiles to write")

	#Stop all queues for writing
	logger.debug("Stop all queues by inserting None")
	for q in qs_list:
		q.put((None, None))

	#Wait for all writer tasks to finish
	finished = 0
	while finished == 0:
		logger.debug("Writer task return status: {0}".format([task.get() if task.ready() else "NA" for task in writer_tasks]))
		if sum([task.ready() for task in writer_tasks]) == len(writer_tasks):	
			finished = 1
			return_codes = [task.get() for task in writer_tasks]
			if sum(return_codes) != 0:
				logger.error("Bedfile writer finished with an error ({0})".format())
			else:
				logger.debug("Bedfile writer(s) finished!")
		time.sleep(0.5) 

	logger.debug("Joining bed_writer queues")
	for i, q in enumerate(qs_list):
		logger.debug("- Queue {0} (size {1})".format(i, q.qsize()))
			
	#Waits until all queues are closed
	writer_pool.join() 

	#-------------------------------------------------------------------------------------------------------------#
	#---------------------------------- Process information on background scores  --------------------------------#
	#-------------------------------------------------------------------------------------------------------------#

	logger.info("Merging results from subsets")
	background = merge_dicts([result[0] for result in results])
	TF_overlaps = merge_dicts([result[1] for result in results])
	results = None

	#Add missing TF overlaps (if some TFs had zero sites)
	for TF1 in motif_list:
		if TF1.prefix not in TF_overlaps:
			TF_overlaps[TF1.prefix] = 0

		for TF2 in motif_list:
			tup = (TF1.prefix, TF2.prefix)
			if tup not in TF_overlaps:
				TF_overlaps[tup] = 0

	#Collect sampled background values
	for bigwig in args.cond_names:
		background["signal"][bigwig] = np.array(background["signal"][bigwig])

	#Check how many values were fetched from background
	n_bg_values = len(background["signal"][args.cond_names[0]])
	logger.debug("Collected {0} values from background".format(n_bg_values))
	if n_bg_values < 1000:
		err_str = "Number of background values collected from peaks is low (={0}) ".format(n_bg_values)
		err_str += "- this affects estimation of the bound/unbound threshold and the normalization between conditions. "
		err_str += "To improve this estimation, please run BINDetect with --peaks = the full peak set across all conditions."
		logger.warning(err_str) 

	#Plot score distribution
	fig = plot_score_distribution([background["signal"][bigwig] for bigwig in args.cond_names], labels=args.cond_names, title="Raw scores per condition")
	figure_pdf.savefig(fig, bbox_inches='tight')
	plt.close()

	#Normalize arrays
	args.norm_objects = {}
	if args.norm_off == True or len(args.cond_names) == 1: #if norm_off or length of cond is 1 - create constant normalization
		for bigwig in args.cond_names:
			args.norm_objects[bigwig] = ArrayNorm("constant", popt=1.0, value_min=0, value_max=1) #no normalization; min/max don't matter for constant norm

	else:
		logger.comment("")
		logger.info("Normalizing scores across conditions")

		list_of_vals = [background["signal"][bigwig] for bigwig in args.cond_names]
		if args.debug: 
			args.norm_objects = quantile_normalization(list_of_vals, args.cond_names, pdfpages=debug_pdf, logger=logger)
		else:
			args.norm_objects = quantile_normalization(list_of_vals, args.cond_names, logger=logger)

		#Normalize background and visualize score distribution
		for bigwig in args.cond_names:
			
			original = background["signal"][bigwig]

			#Check for nan
			logger.debug("Background nans ({0}): {1}".format(bigwig, sum(np.isnan(original))))
			normalized = args.norm_objects[bigwig].normalize(original)
			
			#Replace negative values with 0
			negatives = normalized < 0
			normalized[negatives] = 0

			background["signal"][bigwig] = normalized
			logger.debug("Background nans after normalization ({0}): {1}".format(bigwig, sum(np.isnan(background["signal"][bigwig]))))
		
		fig = plot_score_distribution([background["signal"][bigwig] for bigwig in args.cond_names], labels=args.cond_names, title="Normalized scores per condition")
		figure_pdf.savefig(fig, bbox_inches='tight')
		plt.close()


	#-------------------------------------------------------------------------------------------------------------#
	#-------------------------------------- Estimate bound/unbound threshold -------------------------------------#
	#-------------------------------------------------------------------------------------------------------------#

	logger.info("Estimating bound/unbound threshold")

	#Prepare scores (remove 0's etc.)
	bg_values = np.array([background["signal"][bigwig] for bigwig in args.cond_names]).flatten()	#scores from all conditions
	logger.debug("Size of background array collected: {0}".format(bg_values.size))
	bg_values = bg_values[np.logical_not(np.isclose(bg_values, 0.0))]	#only non-zero counts
	logger.debug("Size of background array after filtering > 0: {0}".format(bg_values.size))
	if len(bg_values) == 0:
		logger.error("Error processing bigwig scores from background. It could be that there are no scores in the bigwig (= all scores are 0) assigned for the peaks. Please check your input files.")
		sys.exit(1)

	x_max = np.percentile(bg_values, [99]) 
	bg_values = bg_values[bg_values < x_max]
	logger.debug("Size of background array after filtering < x_max ({0}): {1}".format(x_max, bg_values.size))

	#Fit mixture of normals
	log_vals = np.log(bg_values).reshape(-1, 1)
	lowest_bic = np.inf
	for n_components in [2]:	#2 components; one for 0's and one for true signal
		gmm = sklearn.mixture.GaussianMixture(n_components=n_components, random_state=1)
		gmm.fit(log_vals)
		
		bic = gmm.bic(log_vals)
		logger.debug("n_compontents: {0} | bic: {1}".format(n_components, bic))
		if bic < lowest_bic:
			lowest_bic = bic
			best_gmm = gmm
	gmm = best_gmm
	
	#Obtain parameters for each component
	means = gmm.means_.flatten()
	stds = np.sqrt(gmm.covariances_).flatten()	

	#Plot components for debugging
	if args.debug:

		fig, ax = plt.subplots(nrows=2, ncols=1, constrained_layout=True)

		#Plot background distribution
		ax[0].hist(log_vals, bins='auto', density=True, color="grey")  #log space
		ax[1].hist(bg_values, bins='auto', density=True, color="grey") #normal space

		#Plot components
		x_log = np.linspace(np.min(log_vals), np.max(log_vals), 1000) 
		x_norm = np.exp(x_log)
		for i in range(len(means)):
			pdf = scipy.stats.norm.pdf(x_log, loc=means[i], scale=stds[i])
			ax[0].plot(x_log, pdf, label="Component {0}".format(i+1))

			#Plot component in normal space
			log_params = scipy.stats.lognorm.fit(bg_values, f0=stds[i], fscale=np.exp(means[i]))
			pdf =  scipy.stats.lognorm.pdf(x_norm, *log_params)
			ax[1].plot(x_norm, pdf, label="Component {0}".format(i+1))

		ax[0].set_title("Background score distribution")
		ax[0].set_xlabel("log(background score)")
		ax[0].set_ylabel("Density")
		ax[0].legend()

		ax[1].set_xlabel("Background score")
		ax[1].set_ylabel("Density")
		ax[1].legend()

		debug_pdf.savefig(fig)
		plt.close()

	#Extract most-right gaussian 
	chosen_i = np.argmax(means) 	#Mixture with largest mean
	log_params = scipy.stats.lognorm.fit(bg_values, f0=stds[chosen_i], fscale=np.exp(means[chosen_i]))

	#Mode of distribution
	mode = scipy.optimize.fmin(lambda x: -scipy.stats.lognorm.pdf(x, *log_params), 0, disp=False)[0]
	logger.debug("- Mode estimated at: {0}".format(mode))
	pseudo = mode / 2.0		#pseudo is half the mode
	args.pseudo = pseudo
	logger.debug("Pseudocount estimated at: {0}".format(round(args.pseudo, 5)))
	
	# Estimate theoretical normal for threshold
	leftside_x = np.linspace(scipy.stats.lognorm(*log_params).ppf([0.01]), mode, 100)
	leftside_pdf = scipy.stats.lognorm.pdf(leftside_x, *log_params)

	#Flip over
	leftside_x_scale = leftside_x - np.min(leftside_x) #scale to min 0
	mirrored_x = np.concatenate([leftside_x, np.max(leftside_x) + leftside_x_scale]).flatten()
	mirrored_pdf = np.concatenate([leftside_pdf, leftside_pdf[::-1]]).flatten()
	popt, cov = scipy.optimize.curve_fit(lambda x, std, sc: sc * scipy.stats.norm.pdf(x, mode, std), mirrored_x, mirrored_pdf)
	norm_params = (mode, popt[0])
	logger.debug("Theoretical normal parameters: {0}".format(norm_params))

	#Set threshold for bound/unbound
	threshold = round(scipy.stats.norm.ppf(1-args.bound_pvalue, *norm_params), 5)

	args.thresholds = {bigwig: threshold for bigwig in args.cond_names}
	logger.stats("- Threshold estimated at: {0}".format(threshold))

	#Only plot if args.debug is True
	if args.debug:

		#Plot mirrored data
		fig, ax = plt.subplots(1,1)
		ax.hist(bg_values[bg_values < x_max], bins='auto', density=True, label="Observed score distribution")
		ax.plot(mirrored_x, mirrored_pdf, color="black")
		plt.xlabel("Bigwig score")
		plt.title("Theoretical normal")
		debug_pdf.savefig(fig)
		plt.close(fig)
		
		#Plot fit and threshold
		fig, ax = plt.subplots(1, 1)
		ax.hist(bg_values[bg_values < x_max], bins='auto', density=True, label="Observed score distribution")

		xvals = np.linspace(0, x_max, 1000)
		log_probas = scipy.stats.lognorm.pdf(xvals, *log_params)
		ax.plot(xvals, log_probas, label="Log-normal fit", color="orange")

		#Theoretical normal
		norm_probas = scipy.stats.norm.pdf(xvals, *norm_params)
		ax.plot(xvals, norm_probas * (np.max(log_probas) / np.max(norm_probas)), color="grey", linestyle="--", label="Theoretical normal")

		ax.axvline(threshold, color="black", label="Bound/unbound threshold")
		ymax = plt.ylim()[1]
		ax.text(threshold, ymax, "\n {0:.3f}".format(threshold), va="top")
		
		#Decorate plot
		plt.title("Score distribution")
		plt.xlabel("Bigwig score")
		plt.ylabel("Density")
		plt.legend(fontsize=8)
		plt.xlim((0,x_max))

		debug_pdf.savefig(fig)
		plt.close(fig)

	#-------------------------------------------------------------------------------------------------------------#
	#--------------------------------------- Foldchanges between conditions --------------------------------------#
	#-------------------------------------------------------------------------------------------------------------#

	logger.comment("")
	log2fc_params = {}
	if len(args.signals) > 1:
		logger.info("Calculating background log2 fold-changes between conditions")

		for (bigwig1, bigwig2) in comparisons:	#cond1, cond2
			logger.info("- {0} / {1}".format(bigwig1, bigwig2))

			#Estimate background log2fc 
			scores1 = np.copy(background["signal"][bigwig1])
			scores2 = np.copy(background["signal"][bigwig2])
			
			included = np.logical_or(scores1 > 0, scores2 > 0)
			scores1 = scores1[included]
			scores2 = scores2[included]

			#Calculate background log2fc normal disitribution
			log2fcs = np.log2(np.true_divide(scores1 + args.pseudo, scores2 + args.pseudo))
			
			lower, upper = np.percentile(log2fcs, [1,99])
			log2fcs_fit = log2fcs[np.logical_and(log2fcs >= lower, log2fcs <= upper)]
			
			#Decide on diff_dist
			diff_dist = scipy.stats.norm
			norm_params = diff_dist.fit(log2fcs_fit)

			logger.debug("({0} / {1}) Background log2fc distribution: {2}".format(bigwig1, bigwig2, norm_params))
			log2fc_params[(bigwig1, bigwig2)] = norm_params

			#If debug: plot background log2fc to figures
			if args.debug:
				fig, ax = plt.subplots(1, 1)
				plt.hist(log2fcs, density=True, bins='auto', label="Background log2fc ({0} / {1})".format(bigwig1, bigwig2))

				xvals = np.linspace(plt.xlim()[0], plt.xlim()[1], 100)
				pdf = diff_dist.pdf(xvals, *log2fc_params[(bigwig1, bigwig2)])	
				plt.plot(xvals, pdf, label="Distribution fit")
				plt.title("Background log2FCs ({0} / {1})".format(bigwig1, bigwig2))
				plt.xlabel("Log2 fold change")
				plt.ylabel("Density")

				debug_pdf.savefig(fig, bbox_inches='tight')
				plt.close()
				
				#f = open(os.path.join(args.outdir, "{0}_{1}_log2fcs.txt".format(bigwig1, bigwig2)), "w")
				#f.write("\n".join([str(val) for val in log2fcs]))
				#f.close()
			
	background = None	 #free up space 

	#-------------------------------------------------------------------------------------------------------------#
	#----------------------------- Read total sites per TF to estimate bound/unbound -----------------------------#
	#-------------------------------------------------------------------------------------------------------------#

	logger.comment("")
	logger.info("Processing scanned TFBS individually")
		
	#Getting bindetect table ready
	info_columns = ["total_tfbs"]
	info_columns.extend(["{0}_{1}".format(cond, metric) for (cond, metric) in itertools.product(args.cond_names, ["threshold", "bound"])])
	info_columns.extend(["{0}_{1}_{2}".format(comparison[0], comparison[1], metric) for (comparison, metric) in itertools.product(comparisons, ["change", "pvalue"])])

	cols = len(info_columns)
	rows = len(motif_names)
	info_table = pd.DataFrame(np.zeros((rows, cols)), columns=info_columns, index=motif_names)

	#Starting calculations
	results = []
	if args.cores == 1:
		for name in motif_names:
			logger.info("- {0}".format(name))
			results.append(process_tfbs(name, args, log2fc_params))
	else:
		logger.debug("Sending jobs to worker pool")

		task_list = [pool.apply_async(process_tfbs, (name, args, log2fc_params, )) for name in motif_names]
		monitor_progress(task_list, logger) 	#will not exit before all jobs are done
		results = [task.get() for task in task_list]

	logger.info("Concatenating results from subsets")
	info_table = pd.concat(results)	 	#pandas tables

	pool.terminate()
	pool.join()
	
	logger.stop_logger_queue()
	
	#-------------------------------------------------------------------------------------------------------------#	
	#------------------------------------------------ Cluster TFBS -----------------------------------------------#	
	#-------------------------------------------------------------------------------------------------------------#	
	
	clustering = RegionCluster(TF_overlaps)
	clustering.cluster(threshold=args.cluster_threshold)

	#Convert full ids to alt ids
	convert = {motif.prefix: motif.name for motif in motif_list}
	for cluster in clustering.clusters:
		for name in convert:
			clustering.clusters[cluster]["cluster_name"] = clustering.clusters[cluster]["cluster_name"].replace(name, convert[name])

	#Write out distance matrix
	matrix_out = os.path.join(args.outdir, args.prefix + "_distances.txt")
	clustering.write_distance_mat(matrix_out)


	#-------------------------------------------------------------------------------------------------------------#	
	#----------------------------------------- Write all_bindetect file ------------------------------------------#
	#-------------------------------------------------------------------------------------------------------------#

	logger.comment("")
	logger.info("Writing all_bindetect files")

	#Add columns of name / motif_id / prefix
	names = []
	ids = []
	for prefix in info_table.index:
		motif = [motif for motif in motif_list if motif.prefix == prefix]
		names.append(motif[0].name)
		ids.append(motif[0].id)

	info_table.insert(0, "output_prefix", info_table.index)
	info_table.insert(1, "name", names)
	info_table.insert(2, "motif_id", ids)

	#info_table.insert(3, "motif_logo", [os.path.join("motif_logos", os.path.basename(logo_filenames[prefix])) for prefix in info_table["output_prefix"]])	#add relative path to logo
	
	#Add cluster to info_table
	cluster_names = []
	for name in info_table.index:
		for cluster in clustering.clusters:
			if name in clustering.clusters[cluster]["member_names"]:
				cluster_names.append(clustering.clusters[cluster]["cluster_name"])
				
	info_table.insert(3, "cluster", cluster_names)
	
	#Cluster table on motif clusters
	info_table_clustered = info_table.groupby("cluster").mean(numeric_only=True) 	#mean of each column
	info_table_clustered.reset_index(inplace=True)

	#Map correct type
	info_table["total_tfbs"] = info_table["total_tfbs"].map(int)
	for condition in args.cond_names:
		info_table[condition + "_bound"] = info_table[condition + "_bound"].map(int)
	
	#Format comparisons
	for (cond1, cond2) in comparisons:
		base = cond1 + "_" + cond2
		info_table[base + "_change"] = info_table[base + "_change"].round(5)
		info_table[base + "_pvalue"] = info_table[base + "_pvalue"].map("{:.5E}".format, na_action="ignore")

		#Define whether TF is going to be highlighted in plot
		names = info_table["output_prefix"]
		changes = info_table[base + "_change"].astype(float)
		pvalues = info_table[base + "_pvalue"].astype(float)
		pval_min = np.percentile(pvalues[pvalues > 0], 5)	#5% smallest pvalues
		change_min, change_max = np.percentile(changes, [5, 95])	#5% smallest and largest changes
		
		#Add "highlighted" information to info_table
		for i, (change, pvalue) in enumerate(zip(changes, pvalues)):
			if change < change_min or change > change_max or pvalue < pval_min:
				info_table.at[names[i], base + "_highlighted"] = True
			else:
				info_table.at[names[i], base + "_highlighted"] = False
	
	#Write bindetect results tables
	#info_table.insert(0, "TF_name", info_table.index)	 #Set index as first column
	bindetect_out = os.path.join(args.outdir, args.prefix + "_results.txt")
	info_table.to_csv(bindetect_out, sep="\t", index=False, header=True, na_rep="NA")

	#### Write excel ###
	bindetect_excel = os.path.join(args.outdir, args.prefix + "_results.xlsx")

	with pd.ExcelWriter(bindetect_excel, engine='xlsxwriter') as writer:

		#Tables
		info_table.to_excel(writer, index=False, sheet_name="Individual motifs")
		info_table_clustered.to_excel(writer, index=False, sheet_name="Motif clusters")

		for sheet in writer.sheets:
			worksheet = writer.sheets[sheet]
			n_rows = worksheet.dim_rowmax
			n_cols = worksheet.dim_colmax
			worksheet.autofilter(0,0,n_rows,n_cols)


	#-------------------------------------------------------------------------------------------------------------#	
	#------------------------------------------- Make BINDetect plot ---------------------------------------------#	
	#-------------------------------------------------------------------------------------------------------------#	

	if no_conditions > 1:
		logger.info("Creating BINDetect plot(s)")

		#Fill NAs from info_table to enable plotting of log2fcs (NA -> 0 change)
		change_cols = [col for col in info_table.columns if "_change" in col]
		pvalue_cols = [col for col in info_table.columns if "_pvalue" in col]
		info_table[change_cols] = info_table[change_cols].fillna(0)
		info_table[pvalue_cols] = info_table[pvalue_cols].fillna(1)

		#Plotting bindetect per comparison
		for (cond1, cond2) in comparisons:

			logger.info("- {0} / {1} (static plot)".format(cond1, cond2))
			base = cond1 + "_" + cond2

			#Fill motifs with metadata (.change, .pvalue, .logpvalue etc.)
			for motif in motif_list:
				name = motif.prefix
				motif.change = float(info_table.at[name, base + "_change"])	#change for this comparison
				motif.pvalue = float(info_table.at[name, base + "_pvalue"])	#pvalue for this comparison
				motif.logpvalue = -np.log10(motif.pvalue) if motif.pvalue > 0 else -np.log10(1e-308)
				motif.highlighted = info_table.at[name, base + "_highlighted"]

				#Assign each motif to group
				if motif.highlighted == True:
					if motif.change < 0:
						motif.group = cond2 + "_up"
					if motif.change > 0:
						motif.group = cond1 + "_up"
				else:
					motif.group = "n.s."

			#Bindetect plot
			fig = plot_bindetect(motif_list, clustering, [cond1, cond2], args)
			figure_pdf.savefig(fig, bbox_inches='tight')
			plt.close(fig)

			#Interactive BINDetect plot
			logger.info("- {0} / {1} (interactive plot)".format(cond1, cond2))
			html_out = os.path.join(args.outdir, "bindetect_" + base + ".html")
			plot_interactive_bindetect(motif_list, [cond1, cond2], html_out)
			

	#-------------------------------------------------------------------------------------------------------------#	
	#----------------------------- Make heatmap across conditions (for debugging)---------------------------------#	
	#-------------------------------------------------------------------------------------------------------------#	

	if args.debug and len(args.signals) > 1:
		logger.info("Plotting heatmap across conditions for debugging")
		mean_columns = [cond + "_mean_score" for cond in args.cond_names]
		heatmap_table = info_table[mean_columns]
		heatmap_table.index = info_table["output_prefix"]

		#Decide fig size
		rows, cols = heatmap_table.shape
		figsize = (7 + cols, max(10, rows/8.0))
		cm = sns.clustermap(heatmap_table,
							figsize = figsize, 
							z_score = 0,		 	#zscore for rows
							col_cluster = False,	#do not cluster condition columns
							yticklabels = True, 		#show all row annotations
							xticklabels = True,
							cbar_pos = (0, 0, .4, .005),
							dendrogram_ratio = (0.3,0.01),
							cbar_kws = {"orientation": "horizontal", 'label': 'Row z-score'},
							method = "single"
							)
							
		#Adjust width of columns
		#hm = cm.ax_heatmap.get_position()
		#cm.ax_heatmap.set_position([hm.x0, hm.y0, cols * 3 * hm.height / rows, hm.height]) 	#aspect should be equal

		plt.setp(cm.ax_heatmap.get_xticklabels(), fontsize=8, rotation=45, ha="right")
		plt.setp(cm.ax_heatmap.get_yticklabels(), fontsize=5)
		
		cm.ax_col_dendrogram.set_title('Mean scores across conditions', fontsize=20)
		cm.ax_heatmap.set_ylabel("Transcription factor motifs", fontsize=15, rotation=270)
		#cm.ax_heatmap.set_title('Conditions')
		#cm.fig.suptitle('Mean scores across conditions')
		#cm.cax.set_visible(False)

		#Save to output pdf
		plt.tight_layout()
		debug_pdf.savefig(cm.fig, bbox_inches='tight')
		plt.close(cm.fig)	


	#-------------------------------------------------------------------------------------------------------------#
	#-------------------------------------------------- Wrap up---------------------------------------------------#
	#-------------------------------------------------------------------------------------------------------------#
	
	if args.debug:
		debug_pdf.close()

	figure_pdf.close()
	logger.end()


#--------------------------------------------------------------------------------------------------------#
if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser = add_bindetect_arguments(parser)
	args = parser.parse_args()

	if len(sys.argv[1:]) == 0:
		parser.print_help()
		sys.exit()

	run_bindetect(args)
