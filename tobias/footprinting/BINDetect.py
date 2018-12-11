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
import logging.handlers
import itertools
import matplotlib.pyplot as plt

from decimal import Decimal
#from sklearn import mixture
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
import scipy

#Bio-specific packages
import pyBigWig
import pysam

#Internal functions and classes
from tobias.footprinting.BINDetect_functions import *
from tobias.utils.utilities import *
from tobias.utils.regions import *
from tobias.utils.sequences import *
from tobias.utils.motifs import *
from tobias.plotting.plot_bindetect import *

np.seterr(divide = 'ignore') 

#--------------------------------------------------------------------------------------------------------------#

def add_bindetect_arguments(parser):

	parser.formatter_class = lambda prog: argparse.RawDescriptionHelpFormatter(prog, max_help_position=35, width=90)
	description = "BINDetect takes motifs, signals (footprints) and genome as input to estimate bound transcription factor binding sites and differential binding between conditions. "
	description += "The underlying method is a modified motif enrichment test to see which motifs have the largest differences in signal across input conditions. "
	description += "The output is an in-depth overview of global changes as well as the individual binding site signal-differences.\n\n"
	description += "Usage:\nTOBIAS BINDetect --signals <bigwig1> (<bigwig2> (...)) --motifs <motifs.txt> --genome <genome.fasta> --peaks <peaks.bed>\n\n"
	description += "Output files:\n- <outdir>/bindetect_figures.pdf\n- <outdir>/bindetect_results.{txt,xlsx}\n- <outdir>/TF_distance_matrix.txt\n"
	description += "- <outdir>/<TF>/<TF>_overview.{txt,xlsx} (per motif)\n- <outdir>/<TF>/beds/<TF>_all.bed (per motif)\n"
	description += "- <outdir>/<TF>/beds/<TF>_<condition>_bound.bed (per motif-condition pair)\n- <outdir>/<TF>/beds/<TF>_<condition>_unbound.bed (per motif-condition pair)\n\n"
	parser.description = format_help_description("BINDetect", description)

	parser._action_groups.pop()	#pop -h
	
	required = parser.add_argument_group('Required arguments')
	required.add_argument('--signals', metavar="<bigwig>", help="Signal per condition (.bigwig format)", nargs="*")
	required.add_argument('--motifs', metavar="<motifs>", help="Motifs in pfm/jaspar format")
	required.add_argument('--genome', metavar="<fasta>", help="Genome .fasta file")
	required.add_argument('--peaks', metavar="<bed>", help="Peaks.bed containing open chromatin regions")

	optargs = parser.add_argument_group('Optional arguments')
	optargs.add_argument('--cond_names', metavar="<name>", nargs="*", help="Names of conditions fitting to --signals (default: prefix of --signals)")
	optargs.add_argument('--peak_header', metavar="<file>", help="File containing the header of --peaks separated by whitespace or newlines (default: peak columns are named \"_additional_<count>\")")
	optargs.add_argument('--naming', metavar="<type>", help="Naming convention for TFs ('id', 'name', 'name_id', 'id_name') (default: 'name_id')", choices=["id", "name", "name_id", "id_name"], default="name_id")
	optargs.add_argument('--motif_pvalue', metavar="<float>", type=lambda x: restricted_float(x, 0, 1), help="Set p-value threshold for motif scanning (default: 1e-4)", default=0.0001)
	optargs.add_argument('--bound_pvalue', metavar="<float>", type=lambda x: restricted_float(x, 0, 1), help="Set p-value threshold for calling bound sites (default: 0.01)", default=0.01)
	optargs.add_argument('--pseudo', type=float, metavar="<float>", help="Pseudocount for calculating log2fcs (default: estimated from data)", default=None)
	optargs.add_argument('--force_overwrite', action="store_true", help="Force overwrite of motif-results that are already existing (default: Only new motif results are calculated)")

	runargs = parser.add_argument_group("Run arguments")
	runargs.add_argument('--outdir', metavar="<directory>", help="Output directory to place TFBS/plots in (default: bindetect_output)", default="bindetect_output")
	runargs.add_argument('--cores', metavar="<int>", type=int, help="Number of cores to use for computation (default: 1)", default=1)
	runargs.add_argument('--split', metavar="<int>", type=int, help="Split of multiprocessing jobs (default: 100)", default=100)
	runargs.add_argument('--verbosity', metavar="<int>", type=int, help="Level of output logging (1 (sparse) / 2 (normal) / 3 (debug)) (default: 2)", choices=[1,2,3], default=2)
	runargs.add_argument('--log', metavar="<file>", help="Full path of logfile (default: log is printed to stdout)")
	runargs.add_argument('--debug', help=argparse.SUPPRESS, action='store_true')
	
	return(parser)


#----------------------------------------------------------------------------------------------------------------#
def run_bindetect(args):

	begin_time = datetime.now()

	#Checking input and setting cond_names
	check_required(args, ["signals", "motifs", "genome", "peaks"])
	args.cond_names = [os.path.basename(os.path.splitext(bw)[0]) for bw in args.signals] if args.cond_names is None else args.cond_names
	args.outdir = os.path.abspath(args.outdir)

	states = ["bound", "unbound"]
	outfiles = [os.path.abspath(os.path.join(args.outdir, "*", "beds", "*_{0}_{1}.bed".format(condition, state))) for (condition, state) in itertools.product(args.cond_names, states)]
	outfiles.append(os.path.abspath(os.path.join(args.outdir, "*", "beds", "*_all.bed")))
	outfiles.append(os.path.abspath(os.path.join(args.outdir, "*", "*_overview.txt")))
	outfiles.append(os.path.abspath(os.path.join(args.outdir, "*", "*_overview.xlsx")))

	outfiles.append(os.path.abspath(os.path.join(args.outdir, "TF_distance_matrix.txt")))
	#outfiles.append(os.path.abspath(os.path.join(args.outdir, "TF_clusters.txt")))
	outfiles.append(os.path.abspath(os.path.join(args.outdir, "bindetect_results.txt")))
	outfiles.append(os.path.abspath(os.path.join(args.outdir, "bindetect_results.xlsx")))
	outfiles.append(os.path.abspath(os.path.join(args.outdir, "bindetect_figures.pdf")))

	
	#----------------------------------------------------------------------------------------------------#
	# Setup logger
	#----------------------------------------------------------------------------------------------------#

	logger = create_logger(args.verbosity, args.log)

	logger.comment("#TOBIAS BINDetect (run started {0})\n".format(begin_time))
	logger.comment("#Command line call: {0}\n".format(" ".join(sys.argv)))

	parser = add_bindetect_arguments(argparse.ArgumentParser())
	logger.comment(arguments_overview(parser, args))

	logger.comment("# ----- Output files -----")
	for outf in outfiles:
		if outf != None:
			logger.comment("# {0}".format(outf))
	logger.comment("\n")

	#-------------------------------------------------------------------------------------------------------------#
	#-------------------------- Pre-processing data: Reading motifs, sequences, peaks ----------------------------#
	#-------------------------------------------------------------------------------------------------------------#

	logger.critical("Processing input data")

	#Check opening/writing of files
	logger.info("Checking reading/writing of files")
	check_files([args.signals, args.motifs, args.genome, args.peaks], action="r")
	check_files(outfiles[-3:], action="w")
	make_directory(args.outdir)

	#Comparisons between conditions
	no_conditions = len(args.signals)
	comparisons = list(itertools.combinations(args.cond_names, 2))

	#Open figure pdf and write overview
	fig_out = os.path.abspath(os.path.join(args.outdir, "bindetect_figures.pdf"))
	figure_pdf = PdfPages(fig_out, keep_empty=True)

	plt.figure()
	plt.axis('off')
	plt.text(0.5,0.8, "BINDETECT FIGURES", ha="center", va="center", fontsize=20)

	#output and order
	titles = []
	for cond in args.cond_names:
		titles.append("Score distribution of {0} scores".format(cond))

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
	logger.info("> Found {0} regions in input peaks".format(len(peaks)))
	peaks = peaks.merge()	#merge overlapping peaks
	logger.info("> Merged to {0} regions".format(len(peaks)))
	
	if args.debug:
		logger.info("Debug on: Peaks reduced to 1000")
		peaks = peaks.subset(10000)
	
	#Make chunks of regions for multiprocessing
	peak_chunks = peaks.chunks(args.split)

	#Header
	if args.peak_header != None:
		content = open(args.peak_header, "r").read()
		args.peak_header_list = content.split()
		logger.debug("Peak header: {0}".format(args.peak_header_list))

		#Check whether peak header fits with number of peak columns
	else:
		args.peak_header_list = None

	
	##### GC content for motif scanning
	fasta_obj = pysam.FastaFile(args.genome)
	logger.info("Estimating GC content from peak sequences") 
	pool = mp.Pool(processes=args.cores)
	gc_content_pool = pool.starmap(get_gc_content, itertools.product(peak_chunks, [args.genome])) 
	gc_content = np.mean(gc_content_pool)	#fraction
	logger.info("> GC content estimated at {0:.2f}%".format(gc_content*100))
	args.gc = gc_content
	bg = np.array([(1-args.gc)/2.0, args.gc/2.0, args.gc/2.0, (1-args.gc)/2.0])


	################ Get motifs ################
	logger.info("Reading motifs from file") 

	motif_content = open(args.motifs).read()
	converted_content = convert_motif(motif_content, "pfm")
	motif_list = pfm_to_motifs(converted_content) 			#List of OneMotif objects
	no_pfms = len(motif_list)
	logger.info("> Found {0} motifs in file".format(no_pfms))

	if args.debug:
		logger.info("Debug on: motifs reduced to 50")
		motif_list = MotifList(motif_list[:50])

	logger.debug("Getting motifs ready")
	motif_list.bg = bg
	motif_names = [motif.name for motif in motif_list]
	logger.debug("Getting reverse motifs")
	motif_list.extend([motif.get_reverse() for motif in motif_list])
	for motif in motif_list:	#now with reverse motifs as well
		motif.set_name(args.naming)
		motif.name = filafy(motif.name)		#remove ()/: etc. which will create problems in filenames
		motif.bg = bg

		logger.debug("Getting pssm for motif {0}".format(motif.name))
		motif.get_pssm()
	
	motif_names = list(set([motif.name for motif in motif_list]))

	#Get threshold for motifs
	logger.debug("Getting match threshold per motif")
	outlist = pool.starmap(OneMotif.get_threshold, itertools.product(motif_list, [args.motif_pvalue])) 
	motif_list = MotifList(outlist)	

	pool.close()
	pool.join()
	

	#-------------------------------------------------------------------------------------------------------------#
	#------------ TF names are known -> test whether they are already calculated -> create subdirs ---------------#
	#-------------------------------------------------------------------------------------------------------------#

	missing_results = []

	#Test whether files for these motifs already exist
	for TF in motif_names:
		overview_file = os.path.join(args.outdir, TF, TF + "_overview.txt") 	#Overview file

		#Check whether file exist
		if not os.path.exists(overview_file):
			logger.debug("Did not find any existing results for {0}".format(TF))
			missing_results.append(TF)
		else:
			logger.debug("Found previous results for {0}".format(TF))

	#Choose which motifs to predict binding of
	if args.force_overwrite == True:
		motif_list_predict = motif_list  # Run bindetect with all input motifs 
		motif_names_predict = list(set([motif.name for motif in motif_list]))
	else:
		#Run bindetect only on missing results
		motif_list_predict = MotifList([motif_obj for motif_obj in motif_list if motif_obj.name in missing_results])
		motif_names_predict = list(set([motif.name for motif in motif_list_predict]))
		logger.info("Prediction run on missing motifs ({0} new motifs). Please use --force_overwrite to rerun for all.".format(len(missing_results)))

	args.new_motifs = motif_names_predict

	logger.info("Creating folder structure for each TF")
	for TF in motif_names:
		make_directory(os.path.join(args.outdir, TF))
		make_directory(os.path.join(args.outdir, TF, "beds"))
		#make_directory(os.path.join(args.outdir, TF, "plots"))

	#-------------------------------------------------------------------------------------------------------------#
	#--------------------- Motif scanning: Find binding sites and match to footprint scores ----------------------#
	#-------------------------------------------------------------------------------------------------------------#

	logger.comment("")
	manager = mp.Manager()

	logger.debug("Setting up listener for log")
	log_q = mp.Manager().Queue()
	listener = mp.Process(target=main_logger_process, args=(log_q, logger))
	listener.start()

	logger.info("Scanning for motifs and matching to signals...")
	worker_cores = max(1, int(args.cores * 0.9))
	writer_cores = max(1, int(args.cores * 0.1))
	logger.debug("Worker cores: {0}".format(worker_cores))
	logger.debug("Writer cores: {0}".format(writer_cores))

	worker_pool = mp.Pool(processes=worker_cores)
	writer_pool = mp.Pool(processes=writer_cores)

	#Create writer queues for bed-file output 
	logger.debug("Setting up writer queues")
	qs_list = []
	qs = {}
	#finished = []
	TF_names_chunks = [motif_names_predict[i::writer_cores] for i in range(writer_cores)]
	for TF_names_sub in TF_names_chunks:
		logger.debug("Creating writer queue for {0}".format(TF_names_sub))
		files = [os.path.join(args.outdir, TF, "beds", TF + "_all.bed") for TF in TF_names_sub]

		q = manager.Queue()
		qs_list.append(q)
		writer_pool.apply_async(file_writer, args=(q, TF_names_sub, files, args)) #, callback = lambda x: finished.append(x) print("Writing time: {0}".format(x)))
		for TF in TF_names_sub:
			qs[TF] = q
	writer_pool.close() #no more jobs applied to writer_pool

	#Start working on data
	if worker_cores == 1:
		logger.info("Running with cores = 1")
		results = []
		for chunk in peak_chunks:
			results.append(scan_and_score(chunk, motif_list, args, log_q, qs))
		
	else: 
		logger.debug("Sending jobs to worker pool")

		task_list = [worker_pool.apply_async(scan_and_score, (chunk, motif_list, args, log_q, qs, )) for chunk in peak_chunks]
		worker_pool.close()
		monitor_progress(task_list, logger)
		worker_pool.join()
		results = [task.get() for task in task_list]
	
	logger.info("Done scanning for TFBS across regions!")

	#---------------------------------------#
	logger.debug("Waiting for listener to finish")
	log_q.put(None)
	while listener.exitcode != 0:
		logger.debug("Listener exitcode is: {0}".format(listener.exitcode))
		time.sleep(1)

	logger.debug("Joining listener")
	listener.join()
	
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


	#--------------------------------------------------------------------------------------#
	#---------------- Process information on background scores and overlaps ---------------#
	#--------------------------------------------------------------------------------------#

	logger.info("Merging results from subsets")
	rand_scores = {}
	TF_overlaps = {}
	for result in results:
		merge_dicts(rand_scores, result[0])
		merge_dicts(TF_overlaps, result[1])
	results = None

	#Mean scores to np array
	logger.info("Estimating score distributions per condition")
	pseudo = []
	args.thresholds = {}
	for bigwig in rand_scores:


		fig, ax = plt.subplots(1, 1)

		rand_scores[bigwig] = np.array(rand_scores[bigwig]) 
		logger.debug("{0} scores for bigwig {1}".format(len(rand_scores[bigwig]), bigwig))

		#Plot histogram
		xmax = np.percentile(rand_scores[bigwig], 99)

		scores = np.copy(rand_scores[bigwig])
		scores = scores[scores > 0]

		to_plot = scores[scores < xmax]
		ax.hist(to_plot, bins=80, density=True, label="Observed score distribution")

		#Fit lognorm to samples
		params = scipy.stats.lognorm.fit(to_plot)
		#print("lognorm params: {0}".format(params))

		xvals = np.linspace(0, xmax, 100)
		probas = scipy.stats.lognorm.pdf(xvals, *params)
		ax.plot(xvals, probas, label="Log-normal fit")

		### Find mode of distribution ###
		xvals = np.linspace(0,max(to_plot),1000)
		vals = scipy.stats.lognorm.pdf(xvals, *params)
		mode = xvals[np.argmax(vals)]
		#print("mode of distribution {0}".format(mode))

		# Estimate theoretical normal 
		a, b = 0.0, mode
		trunc_data = scores[scores <= mode]

		#mirrored
		mirrored = mode + (mode - trunc_data)
		theoretical = np.concatenate((trunc_data, mirrored))
	
		norm_params = scipy.stats.norm.fit(theoretical)
		ax.plot(xvals, scipy.stats.norm.pdf(xvals, *norm_params), color="grey", linestyle="--", label="Theoretical normal")

		#Set threshold for bound/unbound
		threshold = scipy.stats.norm.ppf(1-args.bound_pvalue, *norm_params)

		ax.axvline(threshold, color="black", label="Bound/unbound threshold")
		ymax = plt.ylim()[1]
		ax.text(threshold, ymax, "\n {0:.3f}".format(threshold), va="top")
		args.thresholds[bigwig] = threshold

		#Estimate pseudocount for log2fcs
		pseudo.append(scipy.stats.norm.ppf(0.05, *norm_params))
		
		#Decorate plot
		plt.title("Score distribution of {0} scores".format(bigwig))
		plt.xlabel("Bigwig score")
		plt.ylabel("Density")
		plt.legend()

		#plt.show()
		figure_pdf.savefig()
		plt.close()

	if args.pseudo == None:
		args.pseudo = round(sum(pseudo) / len(pseudo), 5)
		logger.info("Pseudocount estimated at: {0}".format(args.pseudo))

	#Foldchanges between conditions
	log2fc_params = {}
	if len(args.signals) > 1:
		logger.info("Calculating log2 fold changes between conditions")

		for (bigwig1, bigwig2) in comparisons:	#cond1, cond2
			logger.debug("{0} / {1}".format(bigwig1, bigwig2))
			plt.figure()

			log2fc = np.log2(np.true_divide(rand_scores[bigwig1] + args.pseudo, rand_scores[bigwig2] + args.pseudo))
			log2fc = log2fc[np.logical_not(np.isclose(log2fc, 0))]	#remove 0 log2fcs
			n, bins, patches = plt.hist(log2fc, label="log2fc({0} / {1})".format(os.path.basename(bigwig1), os.path.basename(bigwig2)), bins=100, density=True)

			#Curvefit of histogram
			logger.debug("Fitting normal distributions")
			params = list(scipy.stats.norm.fit(log2fc))
	
			logger.debug("({0} / {1} log2fc parameters: {2}".format(bigwig1, bigwig2, log2fc))

			xmin, xmax = plt.xlim()
			x = np.linspace(xmin, xmax, 100)
			p = scipy.stats.norm.pdf(x, *params)

			plt.plot(x, p, linewidth=1, label="Normal distribution fit")
			plt.title("Background log2FCs ({0} / {1})".format(bigwig1, bigwig2))
			plt.xlabel("Log2 fold change")
			plt.ylabel("Density")
			plt.legend()

			figure_pdf.savefig()
			plt.close()

			params.append(len(log2fc))	#add number of sites used to estimate background
			log2fc_params[(bigwig1, bigwig2)] = params

	mean_scores = None

	#-------------------------------------------------------------------------------------------------------------#
	#----------------------------- Read total sites per TF to estimate bound/unbound -----------------------------#
	#-------------------------------------------------------------------------------------------------------------#

	logger.comment("")
	logger.critical("Processing scanned TFBS individually")
	
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
		for name in motif_names_predict:
			logger.info("- {0}".format(name))
			results.append(process_tfbs(name, args, log2fc_params))
	else:
		pool = mp.Pool(processes=args.cores, maxtasksperchild=1)
		task_list = [pool.apply_async(process_tfbs, (name, args, log2fc_params)) for name in motif_names_predict]
		pool.close()	
		monitor_progress(task_list, logger) 	#will not exit before all jobs are done
		pool.terminate()
		pool.join()
		#results = [task.get() for task in task_list]

	#logger.info("Concatenating results from subsets")
	#info_table = pd.concat(results)	 	#pandas tables

	#index_names = info_table.index

	#Check if any TFs are missing, if so: add	
	#index_names = info_table.index
	#for TF in motif_names:
	#	if TF not in index_names:
	#		logger.info("{0} not in index".format(TF))
	#		info_table.append(pd.DataFrame(np.zeros(1,cols), columns=info_columns, index=[TF]))



	#-------------------------------------------------------------------------------------------------------------#	
	#--------------------------- Estimate log2fc effects per TF based on overview file ---------------------------#
	#-------------------------------------------------------------------------------------------------------------#	

	logger.comment("")
	logger.info("Estimating global changes in TF binding")

	pool = mp.Pool(processes=args.cores, maxtasksperchild=1)
	task_list = [pool.apply_async(calculate_global_stats, (name, args, log2fc_params)) for name in motif_names]
	pool.close()	
	monitor_progress(task_list, logger) 	#will not exit before all jobs are done
	pool.terminate()
	pool.join()
	results = [task.get() for task in task_list]


	logger.info("Concatenating results from subsets")
	info_table = pd.concat(results)	 	#pandas tables
	index_names = info_table.index


	#-------------------------------------------------------------------------------------------------------------#	
	#----------------------------------------- Write all_bindetect file ------------------------------------------#
	#-------------------------------------------------------------------------------------------------------------#

	logger.comment("")
	logger.info("Writing all_bindetect files")

	#Condition specific
	info_table["total_tfbs"] = info_table["total_tfbs"].map(int)
	for condition in args.cond_names:
		info_table[condition + "_bound"] = info_table[condition + "_bound"].map(int)
		#info_table[condition + "_threshold"] = info_table[condition + "_threshold"].round(5)

	#### Write excel ###
	bindetect_excel = os.path.join(args.outdir, "bindetect_results.xlsx")
	writer = pd.ExcelWriter(bindetect_excel, engine='xlsxwriter')
	info_table.to_excel(writer)
		
	worksheet = writer.sheets['Sheet1']
	no_rows, no_cols = info_table.shape
	worksheet.autofilter(0,0,no_rows,no_cols)
	writer.save()


	#Format comparisons
	for (cond1, cond2) in comparisons:
		base = cond1 + "_" + cond2
		info_table[base + "_change"] = info_table[base + "_change"].round(5)
		info_table[base + "_pvalue"] = info_table[base + "_pvalue"].map("{:.5E}".format)
	
	#Write bindetect results tables
	bindetect_out = os.path.join(args.outdir, "bindetect_results.txt")
	info_table.insert(0, "TF_name", info_table.index)	 #Set index as first column
	info_table.to_csv(bindetect_out, sep="\t", index=False, header=True, na_rep="NA")


	#-------------------------------------------------------------------------------------------------------------#	
	#------------------------------------------- Make BINDetect plot ---------------------------------------------#	
	#-------------------------------------------------------------------------------------------------------------#	

	if no_conditions > 1:
		logger.info("Creating BINDetect plot(s)")

		#Create distance matrix
		distance_matrix, names = overlap_to_distance(TF_overlaps)
		matrix_out = os.path.join(args.outdir, "TF_distance_matrix.txt")
		np.savetxt(matrix_out, distance_matrix, delimiter="\t", header="\t".join(names), fmt="%.4f")

		#Test index names against names
		index_names = info_table.index
		for name in names:
			if name not in index_names:
				logger.info("{0} not in index".format(name))

		#Plotting bindetect per comparison
		for (cond1, cond2) in comparisons:

			logger.info("- {0} / {1}".format(cond1, cond2))
			
			base = cond1 + "_" + cond2
			changes = [float(info_table.at[name, base + "_change"]) for name in names]
			pvalues = [float(info_table.at[name, base + "_pvalue"]) for name in names]

			#Diffbind plot
			fig = plot_bindetect(names, distance_matrix, changes, pvalues, [cond1, cond2])
			figure_pdf.savefig(fig, bbox_inches='tight')


	#-------------------------------------------------------------------------------------------------------------#	
	#---------------------------------------- Plot all log2fc comparisons ----------------------------------------#	
	#-------------------------------------------------------------------------------------------------------------#	

	"""
	logger.info("Plotting log2fc comparisons per TF")
	
	for i, TF in enumerate(distributions):	#TF name
		logger.info("- {0} ({1} / {2})".format(TF, i+1, len(distributions)))

		fig_out = os.path.join(args.outdir, TF, "plots", TF + "_log2fcs.pdf")
		figure_pdf = PdfPages(fig_out, keep_empty=True)

		for key in distributions[TF]:

			(cond1, cond2) = key
			observed, background = distributions[TF][key]["observed"], distributions[TF][key]["background"]
	
			#Plot distribution comparison
			plt.figure()
			n, bins, patches = plt.hist([observed, background], label=["observed", "background"], color=["red", "black"], alpha=0.5, bins=100, density=True)
			x = np.linspace(bins[0], bins[-1], 100)

			y = scipy.stats.norm.pdf(x, np.mean(observed), np.std(observed))
			plt.plot(x,y, color="red")
			plt.axvline(np.mean(observed), color="red")

			y = scipy.stats.norm.pdf(x, np.mean(background), np.std(background))
			plt.plot(x,y, color="black")
			plt.axvline(np.mean(background), color="black")

			plt.legend()
			plt.xlabel("log2({0} / {1})".format(cond1, cond2))
			plt.ylabel("Density")
			plt.title("{0}".format(TF))
			figure_pdf.savefig()
			plt.close()

		figure_pdf.close()
	"""

	#-------------------------------------------------------------------------------------------------------------#	
	
	figure_pdf.close()

	end_time = datetime.now()
	logger.comment("")
	logger.info("Finished BINDetect run (time elapsed: {0})".format(end_time - begin_time))



#--------------------------------------------------------------------------------------------------------#
if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser = add_bindetect_arguments(parser)
	args = parser.parse_args()

	if len(sys.argv[1:]) == 0:
		parser.print_help()
		sys.exit()

	run_bindetect(args)
