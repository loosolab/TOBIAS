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
from datetime import datetime
import time
from copy import deepcopy
import logging
import itertools
import pandas as pd

#Machine learning
import sklearn
from sklearn import mixture
import scipy
from scipy.optimize import curve_fit

#Plotting
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import NullFormatter

#Bio-specific packages
import pyBigWig
import pysam

#Internal functions and classes
from tobias.footprinting.bindetect_functions import *
from tobias.utils.utilities import *
from tobias.utils.regions import *
from tobias.utils.sequences import *
from tobias.utils.motifs import *
from tobias.utils.logger import * 

#For warnings from curve_fit
import warnings
from scipy.optimize import OptimizeWarning
warnings.simplefilter("ignore", OptimizeWarning)
warnings.simplefilter("ignore", RuntimeWarning)


#--------------------------------------------------------------------------------------------------------------#

def add_bindetect_arguments(parser):

	parser.formatter_class = lambda prog: argparse.RawDescriptionHelpFormatter(prog, max_help_position=35, width=90)
	description = "BINDetect takes motifs, signals (footprints) and genome as input to estimate bound transcription factor binding sites and differential binding between conditions. "
	description += "The underlying method is a modified motif enrichment test to see which motifs have the largest differences in signal across input conditions. "
	description += "The output is an in-depth overview of global changes as well as the individual binding site signal-differences.\n\n"
	description += "Usage:\nTOBIAS BINDetect --signals <bigwig1> (<bigwig2> (...)) --motifs <motifs.txt> --genome <genome.fasta> --peaks <peaks.bed>\n\n"
	description += "Output files:\n- <outdir>/<prefix>_figures.pdf\n- <outdir>/<prefix>_results.{txt,xlsx}\n- <outdir>/<prefix>_distances.txt\n"
	description += "- <outdir>/<TF>/<TF>_overview.{txt,xlsx} (per motif)\n- <outdir>/<TF>/beds/<TF>_all.bed (per motif)\n"
	description += "- <outdir>/<TF>/beds/<TF>_<condition>_bound.bed (per motif-condition pair)\n- <outdir>/<TF>/beds/<TF>_<condition>_unbound.bed (per motif-condition pair)\n\n"
	parser.description = format_help_description("BINDetect", description)

	parser._action_groups.pop()	#pop -h
	
	required = parser.add_argument_group('Required arguments')
	required.add_argument('--signals', metavar="<bigwig>", help="Signal per condition (.bigwig format)", nargs="*")
	required.add_argument('--peaks', metavar="<bed>", help="Peaks.bed containing open chromatin regions across all conditions")
	required.add_argument('--motifs', metavar="<motifs>", help="Motifs in pfm/jaspar format")
	required.add_argument('--genome', metavar="<fasta>", help="Genome .fasta file")

	optargs = parser.add_argument_group('Optional arguments')
	optargs.add_argument('--cond_names', metavar="<name>", nargs="*", help="Names of conditions fitting to --signals (default: prefix of --signals)")
	optargs.add_argument('--peak_header', metavar="<file>", help="File containing the header of --peaks separated by whitespace or newlines (default: peak columns are named \"_additional_<count>\")")
	#optargs.add_argument('--naming', metavar="<type>", help="Naming convention for TFs ('id', 'name', 'name_id', 'id_name') (default: 'name_id')", choices=["id", "name", "name_id", "id_name"], default="name_id")
	optargs.add_argument('--motif_pvalue', metavar="<float>", type=lambda x: restricted_float(x, 0, 1), help="Set p-value threshold for motif scanning (default: 1e-4)", default=0.0001)
	optargs.add_argument('--bound_pvalue', metavar="<float>", type=lambda x: restricted_float(x, 0, 1), help="Set p-value threshold for bound/unbound split (default: 0.001)", default=0.001)
	optargs.add_argument('--pseudo', type=float, metavar="<float>", help="Pseudocount for calculating log2fcs (default: estimated from data)", default=None)
	optargs.add_argument('--time_series', action='store_true', help="Will only compare signals1<->signals2<->signals3 (...) in order of input, and skip all-against-all comparison.")

	runargs = parser.add_argument_group("Run arguments")
	runargs.add_argument('--outdir', metavar="<directory>", help="Output directory to place TFBS/plots in (default: bindetect_output)", default="bindetect_output")
	optargs.add_argument('--prefix', metavar="<prefix>", help="Prefix for overview files in --outdir folder (default: bindetect)", default="bindetect")
	runargs.add_argument('--cores', metavar="<int>", type=int, help="Number of cores to use for computation (default: 1)", default=1)
	runargs.add_argument('--split', metavar="<int>", type=int, help="Split of multiprocessing jobs (default: 100)", default=100)
	runargs.add_argument('--debug', help=argparse.SUPPRESS, action='store_true')
	
	runargs = add_logger_args(runargs)

	return(parser)


def find_nearest_idx(array, value):
    idx = (np.abs(array - value)).argmin()
    return idx

def norm_fit(x, mean, std, scale):
	return(scale * scipy.stats.norm.pdf(x, mean, std))


#----------------------------------------------------------------------------------------------------------------#
def run_bindetect(args):
	""" Main function to run bindetect algorithm with input files and parameters given in args """

	#Checking input and setting cond_names
	check_required(args, ["signals", "motifs", "genome", "peaks"])
	args.cond_names = [os.path.basename(os.path.splitext(bw)[0]) for bw in args.signals] if args.cond_names is None else args.cond_names
	args.outdir = os.path.abspath(args.outdir)
	args.naming = "name_id"

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

	#Open figure pdf and write overview
	fig_out = os.path.abspath(os.path.join(args.outdir, args.prefix + "_figures.pdf"))
	figure_pdf = PdfPages(fig_out, keep_empty=True)

	plt.figure()
	plt.axis('off')
	plt.text(0.5,0.8, "BINDETECT FIGURES", ha="center", va="center", fontsize=20)

	#output and order
	titles = []
	titles.append("Raw score distributions")
	titles.append("Normalized score distributions")
	if args.debug:
		for (cond1, cond2) in comparisons:
			titles.append("Background log2FCs ({0} / {1})".format(cond1, cond2))	

	for (cond1, cond2) in comparisons:
		titles.append("BINDetect plot ({0} / {1})".format(cond1, cond2))

	plt.text(0.1, 0.6, "\n".join(["Page {0}) {1}".format(i+2, titles[i]) for i in range(len(titles))]) + "\n\n", va="top")
	figure_pdf.savefig(bbox_inches='tight')
	plt.close()

	################# Peaks / GC in peaks ################
	#Read peak and peak_header
	peaks = RegionList().from_bed(args.peaks)
	logger.info("- Found {0} regions in input peaks".format(len(peaks)))
	peaks = peaks.merge()	#merge overlapping peaks
	logger.info("- Merged to {0} regions".format(len(peaks)))

	if len(peaks) == 0:
		logger.error("Input --peaks file is empty!")
		sys.exit()
		
	peak_chroms = peaks.get_chroms()
	peak_columns = len(peaks[0]) #number of columns

	#Header
	if args.peak_header != None:
		content = open(args.peak_header, "r").read()
		args.peak_header_list = content.split()
		logger.debug("Peak header: {0}".format(args.peak_header_list))

		#Check whether peak header fits with number of peak columns
		if len(args.peak_header_list) != peak_columns:
			logger.error("Length of --peak_header ({0}) does not fit number of columns in --peaks ({1}).".format(len(args.peak_header_list), peak_columns))
			sys.exit()
	else:
		args.peak_header_list = ["peak_chr", "peak_start", "peak_end"] + ["additional_" + str(num + 1) for num in range(peak_columns-3)]
	logger.debug("Peak header list: {0}".format(args.peak_header_list))

	##### GC content for motif scanning ######
	fasta_obj = pysam.FastaFile(args.genome)
	fasta_chroms = fasta_obj.references

	if not set(peak_chroms).issubset(fasta_chroms):
		logger.warning("Chromosome(s) found in peaks ({0}) are not found in input FASTA file ({1}). These peaks are skipped.".format(peak_chroms, fasta_chroms))
		peaks.keep_chroms(fasta_chroms)	#peaks are changed in place
	
	#Make chunks of regions for multiprocessing
	peak_chunks = peaks.chunks(args.split)

	logger.info("Estimating GC content from peak sequences") 
	gc_content_pool = pool.starmap(get_gc_content, itertools.product(peak_chunks, [args.genome])) 
	gc_content = np.mean(gc_content_pool)	#fraction
	args.gc = gc_content
	bg = np.array([(1-args.gc)/2.0, args.gc/2.0, args.gc/2.0, (1-args.gc)/2.0])
	logger.info("- GC content estimated at {0:.2f}%".format(gc_content*100))


	################ Get motifs ################
	logger.info("Reading motifs from file") 

	motif_content = open(args.motifs).read()
	converted_content = convert_motif(motif_content, "pfm")
	motif_list = pfm_to_motifs(converted_content) 			#List of OneMotif objects
	no_pfms = len(motif_list)

	logger.info("- Found {0} motifs in file".format(no_pfms))
	logger.debug("Getting motifs ready")
	motif_list.bg = bg

	logger.debug("Getting reverse motifs")
	motif_list.extend([motif.get_reverse() for motif in motif_list])
	logger.spam(motif_list)

	#Set prefixes
	for motif in motif_list:	#now with reverse motifs as well
		motif.set_prefix(args.naming)
		motif.bg = bg

		logger.spam("Getting pssm for motif {0}".format(motif.name))
		motif.get_pssm()
	
	motif_names = list(set([motif.prefix for motif in motif_list]))

	#Get threshold for motifs
	logger.debug("Getting match threshold per motif")
	outlist = pool.starmap(OneMotif.get_threshold, itertools.product(motif_list, [args.motif_pvalue])) 

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
	for TF_names_sub in TF_names_chunks:
		logger.debug("Creating writer queue for {0}".format(TF_names_sub))
		files = [os.path.join(args.outdir, TF, "beds", TF + ".tmp") for TF in TF_names_sub]

		q = manager.Queue()
		qs_list.append(q)

		writer_pool.apply_async(file_writer, args=(q, dict(zip(TF_names_sub,files)), args))	 #, callback = lambda x: finished.append(x) print("Writing time: {0}".format(x)))
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
	logger.stop_logger_queue()	#stop the listening process (wait until all was written)
	
	#--------------------------------------#
	logger.info("Waiting for bedfiles to write")

	#Stop all queues for writing
	logger.debug("Stop all queues by inserting None")
	for q in qs_list:
		q.put((None, None))

	logger.debug("Joining bed_writer queues")
	for i, q in enumerate(qs_list):
		logger.debug("- Queue {0} (size {1})".format(i, q.qsize()))
			
	#Waits until all queues are closed
	writer_pool.join() 

	#-------------------------------------------------------------------------------------------------------------#
	#---------------------------- Process information on background scores and overlaps --------------------------#
	#-------------------------------------------------------------------------------------------------------------#

	logger.info("Merging results from subsets")
	background = merge_dicts([result[0] for result in results])
	TF_overlaps = merge_dicts([result[1] for result in results])
	results = None

	for bigwig in args.cond_names:
		background["signal"][bigwig] = np.array(background["signal"][bigwig])

	logger.comment("")
	logger.info("Estimating score distribution per condition")

	fig = plot_score_distribution([background["signal"][bigwig] for bigwig in args.cond_names], labels=args.cond_names, title="Raw scores per condition")
	figure_pdf.savefig(fig, bbox_inches='tight')
	plt.close()

	logger.info("Normalizing scores")
	list_of_vals = [background["signal"][bigwig] for bigwig in args.cond_names]
	normed, norm_objects = quantile_normalization(list_of_vals)

	args.norm_objects = dict(zip(args.cond_names, norm_objects))
	for bigwig in args.cond_names:
		background["signal"][bigwig] = args.norm_objects[bigwig].normalize(background["signal"][bigwig]) 
	
	fig = plot_score_distribution([background["signal"][bigwig] for bigwig in args.cond_names], labels=args.cond_names, title="Normalized scores per condition")
	figure_pdf.savefig(fig, bbox_inches='tight')
	plt.close()

	###########################################################
	logger.info("Estimating bound/unbound threshold")

	#Prepare scores (remove 0's etc.)
	bg_values = np.array(normed).flatten()
	bg_values = bg_values[np.logical_not(np.isclose(bg_values, 0.0))]	#only non-zero counts
	x_max = np.percentile(bg_values, [99]) 
	bg_values = bg_values[bg_values < x_max]
		
	#Fit mixture of normals
	lowest_bic = np.inf
	for n_components in [2]:	#2 components
		gmm = sklearn.mixture.GaussianMixture(n_components=n_components, random_state=1)
		gmm.fit(np.log(bg_values).reshape(-1, 1))
		
		bic = gmm.bic(np.log(bg_values).reshape(-1,1))
		logger.debug("n_compontents: {0} | bic: {1}".format(n_components, bic))
		if bic < lowest_bic:
			lowest_bic = bic
			best_gmm = gmm
	gmm = best_gmm
	
	#Extract most-right gaussian 
	means = gmm.means_.flatten()
	sds = np.sqrt(gmm.covariances_).flatten()	
	chosen_i = np.argmax(means) 	#Mixture with largest mean

	log_params = scipy.stats.lognorm.fit(bg_values[bg_values < x_max], f0=sds[chosen_i], fscale=np.exp(means[chosen_i]))
	#all_log_params[bigwig] = log_params

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
	mirrored_x = np.concatenate([leftside_x, np.max(leftside_x) + leftside_x]).flatten()
	mirrored_pdf = np.concatenate([leftside_pdf, leftside_pdf[::-1]]).flatten()
	popt, cov = curve_fit(lambda x, std, sc: sc * scipy.stats.norm.pdf(x, mode, std), mirrored_x, mirrored_pdf)
	norm_params = (mode, popt[0])
	logger.debug("Theoretical normal parameters: {0}".format(norm_params))

	#Set threshold for bound/unbound
	threshold = round(scipy.stats.norm.ppf(1-args.bound_pvalue, *norm_params), 5)

	args.thresholds = {bigwig: threshold for bigwig in args.cond_names}
	logger.stats("- Threshold estimated at: {0}".format( threshold))

	#Only plot if args.debug is True
	if args.debug:

		#Plot fit
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

		figure_pdf.savefig(fig)
		plt.close(fig)

	############ Foldchanges between conditions ################
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
			
			norm_params = scipy.stats.norm.fit(log2fcs_fit)

			logger.debug("({0} / {1}) Background log2fc normal distribution: {2}".format(bigwig1, bigwig2, norm_params))
			log2fc_params[(bigwig1, bigwig2)] = norm_params

			#Plot background log2fc to figures
			fig, ax = plt.subplots(1, 1)
			plt.hist(log2fcs, density=True, bins='auto', label="Background log2fc ({0} / {1})".format(bigwig1, bigwig2))

			xvals = np.linspace(plt.xlim()[0], plt.xlim()[1], 100)
			pdf = scipy.stats.norm.pdf(xvals, *log2fc_params[(bigwig1, bigwig2)])
			plt.plot(xvals, pdf, label="Normal distribution fit")
			plt.title("Background log2FCs ({0} / {1})".format(bigwig1, bigwig2))

			plt.xlabel("Log2 fold change")
			plt.ylabel("Density")
			if args.debug:
				figure_pdf.savefig(fig, bbox_inches='tight')
			plt.close()			
			
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
		task_list = [pool.apply_async(process_tfbs, (name, args, log2fc_params)) for name in motif_names]
		monitor_progress(task_list, logger) 	#will not exit before all jobs are done
		results = [task.get() for task in task_list]

	logger.info("Concatenating results from subsets")
	info_table = pd.concat(results)	 	#pandas tables

	pool.terminate()
	pool.join()
	
	#-------------------------------------------------------------------------------------------------------------#	
	#------------------------------------------------ Cluster TFBS -----------------------------------------------#	
	#-------------------------------------------------------------------------------------------------------------#	
	
	clustering = RegionCluster(TF_overlaps)
	clustering.cluster()

	#Convert full ids to alt ids
	convert = {motif.prefix:motif.name for motif in motif_list}
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
	
	#Add cluster to info_table
	cluster_names = []
	for name in info_table.index:
		for cluster in clustering.clusters:
			if name in clustering.clusters[cluster]["member_names"]:
				cluster_names.append(clustering.clusters[cluster]["cluster_name"])
	info_table.insert(3, "cluster", cluster_names)
	
	#Cluster table on motif clusters
	info_table_clustered = info_table.groupby("cluster").mean() #mean of each column
	info_table_clustered.reset_index(inplace=True)

	#Map correct type
	info_table["total_tfbs"] = info_table["total_tfbs"].map(int)
	for condition in args.cond_names:
		info_table[condition + "_bound"] = info_table[condition + "_bound"].map(int)

	#### Write excel ###
	bindetect_excel = os.path.join(args.outdir, args.prefix + "_results.xlsx")
	writer = pd.ExcelWriter(bindetect_excel, engine='xlsxwriter')

	#Tables
	info_table.to_excel(writer, index=False, sheet_name="Individual motifs")	
	info_table_clustered.to_excel(writer, index=False, sheet_name="Motif clusters")
	
	for sheet in writer.sheets:
		worksheet = writer.sheets[sheet]
		n_rows = worksheet.dim_rowmax
		n_cols = worksheet.dim_colmax
		worksheet.autofilter(0,0,n_rows,n_cols)

	writer.save()

	#Format comparisons
	for (cond1, cond2) in comparisons:
		base = cond1 + "_" + cond2
		info_table[base + "_change"] = info_table[base + "_change"].round(5)
		info_table[base + "_pvalue"] = info_table[base + "_pvalue"].map("{:.5E}".format)
	
	#Write bindetect results tables
	#info_table.insert(0, "TF_name", info_table.index)	 #Set index as first column
	bindetect_out = os.path.join(args.outdir, args.prefix + "_results.txt")
	info_table.to_csv(bindetect_out, sep="\t", index=False, header=True, na_rep="NA")


	#-------------------------------------------------------------------------------------------------------------#	
	#------------------------------------------- Make BINDetect plot ---------------------------------------------#	
	#-------------------------------------------------------------------------------------------------------------#	

	if no_conditions > 1:
		logger.info("Creating BINDetect plot(s)")

		#Plotting bindetect per comparison
		for (cond1, cond2) in comparisons:

			logger.info("- {0} / {1}".format(cond1, cond2))
			base = cond1 + "_" + cond2

			#Make copy of motifs and fill in with metadata
			comparison_motifs = motif_list 	#copy.deepcopy(motif_list) - swig pickle error, just overwrite motif_list
			for motif in comparison_motifs:
				name = motif.prefix
				motif.change = float(info_table.at[name, base + "_change"])
				motif.pvalue = float(info_table.at[name, base + "_pvalue"])

			#Bindetect plot
			fig = plot_bindetect(comparison_motifs, clustering, [cond1, cond2], args)
			figure_pdf.savefig(fig, bbox_inches='tight')


	#-------------------------------------------------------------------------------------------------------------#
	#-------------------------------------------------- Wrap up---------------------------------------------------#
	#-------------------------------------------------------------------------------------------------------------#
	
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
