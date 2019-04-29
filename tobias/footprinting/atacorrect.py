#!/usr/bin/env python

"""
ATACorrect.py: Estimates ATAC-seq bias and corrects read counts from .bam and .fasta input

@author: Mette Bentsen
@contact: mette.bentsen (at) mpi-bn.mpg.de
@license: MIT

"""

#--------------------------------------------------------------------------------------------------------#
#----------------------------------------- Import libraries ---------------------------------------------# 
#--------------------------------------------------------------------------------------------------------#

import os
os.environ["MKL_NUM_THREADS"] = "1" 
os.environ["NUMEXPR_NUM_THREADS"] = "1" 
os.environ["OMP_NUM_THREADS"] = "1" 

import sys
import argparse
import numpy as np
import multiprocessing as mp
from datetime import datetime
from copy import deepcopy
import gc

import textwrap
from collections import OrderedDict
import logging
import itertools
from matplotlib.backends.backend_pdf import PdfPages

#Bio-specific packages
import pyBigWig
import pysam

#Internal functions and classes
from tobias.footprinting.atacorrect_functions import *
from tobias.utils.utilities import *
from tobias.utils.regions import *
from tobias.utils.sequences import *
from tobias.utils.ngs import *
from tobias.utils.logger import *

#np.seterr(divide='raise', invalid='raise')

#--------------------------------------------------------------------------------------------------------#
#----------------------------------------- Argument parser ----------------------------------------------# 
#--------------------------------------------------------------------------------------------------------#

def add_atacorrect_arguments(parser):

	parser.formatter_class = lambda prog: argparse.RawDescriptionHelpFormatter(prog, max_help_position=35, width=90)

	description = "ATACorrect corrects the cutsite-signal from ATAC-seq with regard to the underlying sequence preference of Tn5 transposase.\n\n"
	description += "Usage:\nTOBIAS ATACorrect --bam <reads.bam> --genome <genome.fa> --peaks <peaks.bed>\n\n"
	description += "Output files:\n"
	description += "\n".join(["- <outdir>/<prefix>_{0}.bw".format(track) for track in ["uncorrected", "bias", "expected", "corrected"]]) + "\n"
	description += "- <outdir>/<prefix>_atacorrect.pdf"
	parser.description = format_help_description("ATACorrect", description)

	parser._action_groups.pop()	#pop -h

	#Required arguments
	reqargs = parser.add_argument_group('Required arguments')
	reqargs.add_argument('-b', '--bam', metavar="<bam>", help="A .bam-file containing reads to be corrected")
	reqargs.add_argument('-g', '--genome', metavar="<fasta>", help="A .fasta-file containing whole genomic sequence")
	reqargs.add_argument('-p', '--peaks', metavar="<bed>", help="A .bed-file containing ATAC peak regions")

	#Optional arguments
	optargs = parser.add_argument_group('Optional arguments')
	optargs.add_argument('--regions_in', metavar="<bed>", help="Input regions for estimating bias (default: regions not in peaks.bed)")
	optargs.add_argument('--regions_out', metavar="<bed>", help="Output regions (default: peaks.bed)")
	optargs.add_argument('--blacklist', metavar="<bed>", help="Blacklisted regions in .bed-format (default: None)") #file containing blacklisted regions to be excluded from analysis")
	optargs.add_argument('--extend', metavar="<int>", type=int, help="Extend output regions with basepairs upstream/downstream (default: 100)", default=100)
	optargs.add_argument('--split_strands', help="Write out tracks per strand", action="store_true")
	optargs.add_argument('--norm_off', help="Switches off normalization based on number of reads", action='store_true')
	optargs.add_argument('--track_off', metavar="<track>", help="Switch off writing of individual .bigwig-tracks (uncorrected/bias/expected/corrected)", nargs="*", choices=["uncorrected", "bias", "expected", "corrected"], default=[])

	optargs = parser.add_argument_group('Advanced ATACorrect arguments (no need to touch)')
	optargs.add_argument('--k_flank', metavar="<int>", help="Flank +/- of cutsite to estimate bias from (default: 12)", type=int, default=12)
	optargs.add_argument('--read_shift', metavar="<int>", help="Read shift for forward and reverse reads (default: 4 -5)", nargs=2, type=int, default=[4,-5])
	optargs.add_argument('--bg_shift', metavar="<int>", type=int, help="Read shift for estimation of background frequencies (default: 100)", default=100)
	optargs.add_argument('--window', metavar="<int>", help="Window size for calculating expected signal (default: 100)", type=int, default=100)
	optargs.add_argument('--score_mat', metavar="<mat>", help="Type of matrix to use for bias estimation (PWM/DWM) (default: DWM)", choices=["PWM", "DWM"], default="DWM")

	runargs = parser.add_argument_group('Run arguments')
	runargs.add_argument('--prefix', metavar="<prefix>", help="Prefix for output files (default: same as .bam file)")
	runargs.add_argument('--outdir', metavar="<directory>", help="Output directory for files (default: current working directory)", default="")
	runargs.add_argument('--cores', metavar="<int>", type=int, help="Number of cores to use for computation (default: 1)", default=1)
	runargs.add_argument('--split', metavar="<int>", type=int, help="Split of multiprocessing jobs (default: 100)", default=100)
	
	runargs = add_logger_args(runargs)

	return(parser)

#--------------------------------------------------------------------------------------------------------#
#-------------------------------------- Main pipeline function ------------------------------------------# 
#--------------------------------------------------------------------------------------------------------#

def run_atacorrect(args):

	"""
	Function for bias correction of input .bam files
	Calls functions in ATACorrect_functions and several internal classes
	"""

	#Test if required arguments were given:
	if args.bam == None:
		sys.exit("Error: No .bam-file given")
	if args.genome == None:
		sys.exit("Error: No .fasta-file given")
	if args.peaks == None:
		sys.exit("Error: No .peaks-file given")

	#Adjust some parameters depending on input
	args.prefix = os.path.splitext(os.path.basename(args.bam))[0] if args.prefix == None else args.prefix
	args.outdir = os.path.abspath(args.outdir) if args.outdir != None else os.path.abspath(os.getcwd())

	#Set output bigwigs based on input
	tracks = ["uncorrected", "bias", "expected", "corrected"]
	tracks = [track for track in tracks if track not in args.track_off] 	# switch off printing

	if args.split_strands == True:
		strands = ["forward", "reverse"]
	else:
		strands = ["both"]

	output_bws = {}
	for track in tracks:
		output_bws[track] = {}
		for strand in strands:
			elements = [args.prefix, track] if strand == "both" else [args.prefix, track, strand]
			output_bws[track][strand] = {"fn": os.path.join(args.outdir, "{0}.bw".format("_".join(elements)))}

	#Set all output files
	bam_out = os.path.join(args.outdir, args.prefix + "_atacorrect.bam") 
	bigwigs = [output_bws[track][strand]["fn"] for (track, strand) in itertools.product(tracks, strands)]
	figures_f = os.path.join(args.outdir, "{0}_atacorrect.pdf".format(args.prefix))
	
	output_files = bigwigs + [figures_f]
	output_files = list(OrderedDict.fromkeys(output_files)) 	#remove duplicates due to "both" option

	strands = ["forward", "reverse"]

	#----------------------------------------------------------------------------------------------------#
	# Print info on run
	#----------------------------------------------------------------------------------------------------#

	logger = TobiasLogger("ATACorrect", args.verbosity)
	logger.begin()

	parser = add_atacorrect_arguments(argparse.ArgumentParser())
	logger.arguments_overview(parser, args)
	logger.output_files(output_files)

	#----------------------------------------------------------------------------------------------------#
	# Test input file availability for reading 
	#----------------------------------------------------------------------------------------------------#

	logger.info("----- Processing input data -----")

	logger.debug("Testing input file availability")
	check_files([args.bam, args.genome, args.peaks], "r")

	logger.debug("Testing output directory/file writeability")
	make_directory(args.outdir)
	check_files(output_files, "w")

	#Open pdf for figures
	figure_pdf = PdfPages(figures_f, keep_empty=True)

	#----------------------------------------------------------------------------------------------------#
	# Read information in bam/fasta
	#----------------------------------------------------------------------------------------------------#

	logger.info("Reading info from .bam file")
	bamfile = pysam.AlignmentFile(args.bam, "rb")
	if bamfile.has_index() == False:
		logger.warning("No index found for bamfile - creating one via pysam.")
		pysam.index(args.bam)

	bam_references = bamfile.references 	#chromosomes in correct order
	bam_chrom_info = dict(zip(bamfile.references, bamfile.lengths))
	bamfile.close()

	logger.info("Reading info from .fasta file")
	fastafile = pysam.FastaFile(args.genome)
	fasta_chrom_info = dict(zip(fastafile.references, fastafile.lengths))
	fastafile.close()

	#Compare chrom lengths
	chrom_in_common = set(bam_chrom_info.keys()).intersection(fasta_chrom_info.keys())
	for chrom in chrom_in_common:
		bamlen = bam_chrom_info[chrom]
		fastalen = fasta_chrom_info[chrom]
		if bamlen != fastalen:
			logger.warning("(Fastafile)\t{0} has length {1}".format(chrom, fasta_chrom_info[chrom]))
			logger.warning("(Bamfile)\t{0} has length {1}".format(chrom, bam_chrom_info[chrom]))
			sys.exit("Error: .bam and .fasta have different chromosome lengths. Please make sure the genome file is similar to the one used in mapping.")


	#----------------------------------------------------------------------------------------------------#
	# Read regions from bedfiles
	#----------------------------------------------------------------------------------------------------#

	logger.info("Processing input/output regions")

	#Chromosomes included in analysis
	genome_regions = RegionList().from_list([OneRegion([chrom, 0, bam_chrom_info[chrom]]) for chrom in bam_references if not "M" in chrom])		#full genome length
	chrom_in_common = [chrom for chrom in chrom_in_common if "M" not in chrom]
	logger.debug("CHROMS\t{0}".format("; ".join(["{0} ({1})".format(reg.chrom, reg.end) for reg in genome_regions])))
	genome_bp = sum([region.get_length() for region in genome_regions])

	# Process peaks
	peak_regions = RegionList().from_bed(args.peaks)
	peak_regions.merge()
	peak_regions.apply_method(OneRegion.check_boundary, bam_chrom_info, "cut")
	nonpeak_regions = deepcopy(genome_regions).subtract(peak_regions)

	# Process specific input regions if given
	if args.regions_in != None:
		input_regions = RegionList().from_bed(args.regions_in)
		input_regions.merge()
		input_regions.apply_method(OneRegion.check_boundary, bam_chrom_info, "cut")
	else:
		input_regions = nonpeak_regions

	# Process specific output regions
	if args.regions_out != None:
		output_regions = RegionList().from_bed(args.regions_out)
	else:
		output_regions = deepcopy(peak_regions)
	
	output_regions.apply_method(OneRegion.extend_reg, args.extend)
	output_regions.merge()	
	output_regions.apply_method(OneRegion.check_boundary, bam_chrom_info, "cut")

	#Remove blacklisted regions and chromosomes not in common
	blacklist_regions = RegionList().from_bed(args.blacklist) if args.blacklist != None else RegionList([])	 #fill in with regions from args.blacklist
	regions_dict = {"genome": genome_regions, "input_regions":input_regions, "output_regions":output_regions, "peak_regions":peak_regions, "nonpeak_regions":nonpeak_regions, "blacklist_regions": blacklist_regions}
	for sub in ["input_regions", "output_regions", "peak_regions", "nonpeak_regions"]:
		regions_sub = regions_dict[sub]
		regions_sub.subtract(blacklist_regions)
		regions_sub = regions_sub.apply_method(OneRegion.split_region, 50000)

		regions_sub.keep_chroms(chrom_in_common)
		regions_dict[sub] = regions_sub
	
	#write beds to look at in igv
	#input_regions.write_bed(os.path.join(args.outdir, "input_regions.bed"))
	#output_regions.write_bed(os.path.join(args.outdir, "output_regions.bed"))
	#peak_regions.write_bed(os.path.join(args.outdir, "peak_regions.bed"))
	#nonpeak_regions.write_bed(os.path.join(args.outdir, "nonpeak_regions.bed"))

	#Sort according to order in bam_references:
	output_regions.loc_sort(bam_references)
	chrom_order = {bam_references[i]:i for i in range(len(bam_references))}	 #for use later when sorting output 

	#### Statistics about regions ####
	genome_bp = sum([region.get_length() for region in regions_dict["genome"]])
	for key in regions_dict:
		total_bp = sum([region.get_length() for region in regions_dict[key]])
		logger.stats("{0}: {1} regions | {2} bp | {3:.2f}% coverage".format(key, len(regions_dict[key]), total_bp, total_bp/genome_bp*100))

	#Estallish variables for regions to be used
	input_regions = regions_dict["input_regions"]
	output_regions = regions_dict["output_regions"]
	peak_regions = regions_dict["peak_regions"]
	nonpeak_regions = regions_dict["nonpeak_regions"]

	#Exit if no input/output regions were found
	if len(input_regions) == 0 or len(output_regions) == 0 or len(peak_regions) == 0 or len(nonpeak_regions) == 0:
		logger.error("No regions found - exiting!")
		sys.exit()

	#----------------------------------------------------------------------------------------------------#
	# Estimate normalization factors
	#----------------------------------------------------------------------------------------------------#

	#Setup logger queue
	logger.debug("Setting up listener for log")
	logger.start_logger_queue()
	args.log_q = logger.queue

	#----------------------------------------------------------------------------------------------------#

	logger.comment("")
	logger.info("----- Estimating normalization factors -----")

	#If normalization is to be calculated
	if not args.norm_off:

		#Reads in peaks/nonpeaks
		logger.info("Counting reads in peak regions")
		peak_region_chunks = peak_regions.chunks(args.split)
		reads_peaks = sum(run_parallel(count_reads, peak_region_chunks, [args], args.cores, logger))
		logger.comment("")

		logger.info("Counting reads in nonpeak regions")
		nonpeak_region_chunks = nonpeak_regions.chunks(args.split)
		reads_nonpeaks = sum(run_parallel(count_reads, nonpeak_region_chunks, [args], args.cores, logger))

		reads_total = reads_peaks + reads_nonpeaks

		logger.stats("TOTAL_READS\t{0}".format(reads_total))
		logger.stats("PEAK_READS\t{0}".format(reads_peaks))
		logger.stats("NONPEAK_READS\t{0}".format(reads_nonpeaks))

		lib_norm = 10000000/reads_total
		frip = reads_peaks/reads_total
		correct_factor = lib_norm*(1/frip)

		logger.stats("LIB_NORM\t{0:.5f}".format(lib_norm))
		logger.stats("FRiP\t{0:.5f}".format(frip))
	else:
		logger.info("Normalization was switched off")
		correct_factor = 1.0

	logger.stats("CORRECTION_FACTOR:\t{0:.5f}".format(correct_factor))

	#----------------------------------------------------------------------------------------------------#
	# Estimate sequence bias
	#----------------------------------------------------------------------------------------------------#

	logger.comment("")
	logger.info("Started estimation of sequence bias...")

	input_region_chunks = input_regions.chunks(args.split)										#split to 100 chunks (also decides the step of output)
	out_lst = run_parallel(bias_estimation, input_region_chunks, [args], args.cores, logger)	#Output is list of AtacBias objects

	#Join objects
	estimated_bias = out_lst[0]		#initialize object with first output
	for output in out_lst[1:]:
		estimated_bias.join(output)		#bias object contains bias/background SequenceMatrix objects

	#----------------------------------------------------------------------------------------------------#
	# Join estimations from all chunks of regions
	#----------------------------------------------------------------------------------------------------#

	bias_obj = estimated_bias
	bias_obj.correction_factor = correct_factor

	### Bias motif ###
	logger.info("Finalizing bias motif for scoring")
	for strand in strands:
		bias_obj.bias[strand].prepare_mat()
		figure_pdf.savefig(plot_pssm(bias_obj.bias[strand].pssm, "Tn5 insertion bias of reads ({0})".format(strand)))
	
	#----------------------------------------------------------------------------------------------------#
	# Correct read bias and write to bigwig
	#----------------------------------------------------------------------------------------------------#

	logger.comment("")
	logger.info("----- Correcting reads from .bam within output regions -----")

	output_regions.loc_sort(bam_references)		#sort in order of references
	output_regions_chunks = output_regions.chunks(args.split)
	no_tasks = float(len(output_regions_chunks))
	chunk_sizes = [len(chunk) for chunk in output_regions_chunks]
	logger.debug("All regions chunked: {0} ({1})".format(len(output_regions), chunk_sizes))

	### Create key-file linking for bigwigs 
	key2file = {}
	for track in output_bws:
		for strand in output_bws[track]:
			filename = output_bws[track][strand]["fn"]
			key = "{}:{}".format(track, strand)
			key2file[key] = filename

	#Start correction/write cores
	n_bigwig = len(key2file.values())
	args.cores = check_cores(args.cores, logger)
	writer_cores = min(n_bigwig, max(1,int(args.cores*0.1)))	#at most one core per bigwig or 10% of cores (or 1)
	worker_cores = max(1, args.cores - writer_cores) 				
	logger.debug("Worker cores: {0}".format(worker_cores))
	logger.debug("Writer cores: {0}".format(writer_cores))

	worker_pool = mp.Pool(processes=worker_cores)
	writer_pool = mp.Pool(processes=writer_cores)
	manager = mp.Manager()

	#Start bigwig file writers
	header = [(chrom, bam_chrom_info[chrom]) for chrom in bam_references]
	key_chunks = [list(key2file.keys())[i::writer_cores] for i in range(writer_cores)]
	qs_list = []
	qs = {}
	for chunk in key_chunks:
		logger.debug("Creating writer queue for {0}".format(chunk))

		q = manager.Queue()
		qs_list.append(q)

		files = [key2file[key] for key in chunk]
		writer_pool.apply_async(bigwig_writer, args=(q, dict(zip(chunk, files)), header, output_regions, args))	 #, callback = lambda x: finished.append(x) print("Writing time: {0}".format(x)))
		for key in chunk:
			qs[key] = q

	args.qs = qs
	writer_pool.close() #no more jobs applied to writer_pool

	#Start correction
	logger.debug("Starting correction")
	task_list = [worker_pool.apply_async(bias_correction, args=[chunk, args, bias_obj]) for chunk in output_regions_chunks]
	worker_pool.close()
	monitor_progress(task_list, logger, "Correction progress:")	#does not exit until tasks in task_list finished
	results = [task.get() for task in task_list]

	#Get all results 
	pre_bias = results[0][0]	#initialize with first result
	post_bias = results[0][1]	#initialize with first result
	for result in results[1:]:
		pre_bias_chunk = result[0]
		post_bias_chunk = result[1]

		for direction in strands:
			pre_bias[direction].add_counts(pre_bias_chunk[direction])
			post_bias[direction].add_counts(post_bias_chunk[direction])

	#Stop all queues for writing
	logger.debug("Stop all queues by inserting None")
	for q in qs_list:
		q.put((None, None, None))

	logger.debug("Joining bigwig_writer queues")
	
	qsum = sum([q.qsize() for q in qs_list])
	while qsum != 0:
		qsum = sum([q.qsize() for q in qs_list])
		logger.spam("- Queue sizes {0}".format([(key, qs[key].qsize()) for key in qs]))
		time.sleep(0.5)

	#Waits until all queues are closed
	writer_pool.join() 
	worker_pool.terminate()
	worker_pool.join()

	#Stop multiprocessing logger	
	logger.stop_logger_queue()

	#----------------------------------------------------------------------------------------------------#
	# Information and verification of corrected read frequencies
	#----------------------------------------------------------------------------------------------------#		

	logger.comment("")
	logger.info("Verifying bias correction")

	#Calculating variance per base
	for strand in strands:

		#Invert negative counts
		abssum = np.abs(np.sum(post_bias[strand].neg_counts, axis=0))
		post_bias[strand].neg_counts = post_bias[strand].neg_counts + abssum
		
		#Join negative/positive counts
		post_bias[strand].counts += post_bias[strand].neg_counts	#now pos

		pre_bias[strand].prepare_mat()
		post_bias[strand].prepare_mat()

		pre_var = np.mean(np.var(pre_bias[strand].bias_pwm, axis=1)[:4])   #mean of variance per nucleotide
		post_var = np.mean(np.var(post_bias[strand].bias_pwm, axis=1)[:4])
		logger.stats("BIAS\tpre-bias variance {0}:\t{1:.7f}".format(strand, pre_var))
		logger.stats("BIAS\tpost-bias variance {0}:\t{1:.7f}".format(strand, post_var))

		#Plot figure
		fig_title = "Nucleotide frequencies in corrected reads\n({0} strand)".format(strand)
		figure_pdf.savefig(plot_correction(pre_bias[strand].bias_pwm, post_bias[strand].bias_pwm, fig_title))

	
	#----------------------------------------------------------------------------------------------------#
	# Finish up
	#----------------------------------------------------------------------------------------------------#

	figure_pdf.close()
	logger.end()

#--------------------------------------------------------------------------------------------------------#
if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser = add_atacorrect_arguments(parser)
	args = parser.parse_args()

	if len(sys.argv[1:]) == 0:
		parser.print_help()
		sys.exit()

	run_atacorrect(args)
