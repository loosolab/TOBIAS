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
from scipy.optimize import curve_fit
from matplotlib.backends.backend_pdf import PdfPages

#Bio-specific packages
import pyBigWig
import pysam

#Internal functions and classes
from tobias.footprinting.ATACorrect_functions import *
from tobias.utils.utilities import *
from tobias.utils.regions import *
from tobias.utils.sequences import *
from tobias.utils.ngs import *

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
	optargs.add_argument('--split_strands', help="Calculate and correct bias separately per strand", action="store_true")
	optargs.add_argument('--norm_off', help="Switches off normalization based on number of reads", action='store_true')
	optargs.add_argument('--track_off', metavar="<track>", help="Switch off writing of individual .bigwig-tracks (uncorrected/bias/expected/corrected)", nargs="*", choices=["uncorrected", "bias", "expected", "corrected"], default=[])

	optargs = parser.add_argument_group('Advanced ATACorrect arguments (no need to touch)')
	optargs.add_argument('--k_flank', metavar="<int>", help="Flank +/- of cutsite to estimate bias from (default: 12)", type=int, default=12)
	optargs.add_argument('--read_shift', metavar="<int>", help="Read shift for forward and reverse reads (default: 4 -5)", nargs=2, type=int, default=[4,-5])
	optargs.add_argument('--bg_shift', metavar="<int>", type=int, help="Read shift for estimation of background frequencies (default: 100)", default=100)
	optargs.add_argument('--window', metavar="<int>", help="Window size for calculating expected signal (default: 100)", type=int, default=100)
	optargs.add_argument('--score_mat', metavar="<mat>", help="Type of matrix to use for bias estimation (PWM/DWM) (default: DWM)", choices=["PWM", "DWM"], default="DWM")
	#optargs.add_argument('-m', '--method', metavar="", help="Method for correcting cutsites (subtract/log2fc) (default: subtract)", choices=["log2fc", "subtract"], default="subtract")

	runargs = parser.add_argument_group('Run arguments')
	runargs.add_argument('--prefix', metavar="<prefix>", help="Prefix for output files (default: same as .bam file)")
	runargs.add_argument('--outdir', metavar="<directory>", help="Output directory for files (default: current working directory)", default="")
	
	##shared across TOBIAS
	runargs.add_argument('--cores', metavar="<int>", type=int, help="Number of cores to use for computation (default: 1)", default=1)
	runargs.add_argument('--split', metavar="<int>", type=int, help="Split of multiprocessing jobs (default: 100)", default=100)
	runargs.add_argument('--verbosity', metavar="<int>", help="Level of output logging (1 (sparse) / 2 (normal) / 3 (debug)) (default: 2)", choices=[1,2,3], default=2, type=int)
	runargs.add_argument('--log', metavar="<file>", help="Full path of logfile (default: log is printed to stdout)")
	runargs.add_argument('--debug', help=argparse.SUPPRESS, action='store_true')

	return(parser)

#--------------------------------------------------------------------------------------------------------#
#-------------------------------------- Main pipeline function ------------------------------------------# 
#--------------------------------------------------------------------------------------------------------#

def run_atacorrect(args):

	"""
	Function for bias correction
	Calls functions in ATACorrect_functions and several internal classes
	"""

	begin_time = datetime.now()

	#Test if required arguments were given:
	if args.bam == None:
		sys.exit("Error: No .bam-file given")
	if args.genome == None:
		sys.exit("Error: No .fasta-file given")
	if args.peaks == None:
		sys.exit("Error: No .peaks-file given")

	#Adjust input files to full path
	args.bam = os.path.abspath(args.bam)
	args.genome = os.path.abspath(args.genome)
	args.peaks = os.path.abspath(args.peaks) if args.peaks != None else None

	#Adjust some parameters depending on input
	args.prefix = os.path.splitext(os.path.basename(args.bam))[0] if args.prefix == None else args.prefix
	args.outdir = os.path.abspath(args.outdir) if args.outdir != None else os.path.abspath(os.getcwd())
	args.log = os.path.abspath(args.log) if args.log != None else None

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
		for strand in strands + ["both"]:
			elements = [args.prefix, track] if strand == "both" else [args.prefix, track, strand]
			output_bws[track][strand] = {"fn": os.path.join(args.outdir, "{0}.bw".format("_".join(elements)))}

	#Set all output files
	bam_out = os.path.join(args.outdir, args.prefix + "_atacorrect.bam") 
	bigwigs = [output_bws[track][strand]["fn"] for (track, strand) in itertools.product(tracks, strands + ["both"])]
	figures_f = os.path.join(args.outdir, "{0}_atacorrect.pdf".format(args.prefix))
	log_f = args.log
	
	output_files = bigwigs + [figures_f, log_f]
	output_files = list(OrderedDict.fromkeys(output_files)) 	#remove duplicates due to "both" option


	#----------------------------------------------------------------------------------------------------#
	# Print info on run
	#----------------------------------------------------------------------------------------------------#

	logger = create_logger(args.verbosity, args.log)

	logger.comment("#TOBIAS ATACorrect (run started {0})\n".format(begin_time))
	logger.comment("#Command line call: {0}\n".format(" ".join(sys.argv)))

	parser = add_atacorrect_arguments(argparse.ArgumentParser())
	logger.comment(arguments_overview(parser, args))

	logger.comment("# ----- Output files -----")
	for outf in output_files:
		if outf != None:
			logger.comment("# {0}".format(outf))
	logger.comment("\n\n")


	#----------------------------------------------------------------------------------------------------#
	# Test input file availability for reading 
	#----------------------------------------------------------------------------------------------------#

	logger.critical("----- Processing input data -----")

	#Input test
	logger.info("Testing input file availability")
	file_list = [args.bam, args.genome, args.peaks]
	file_list = [file for file in file_list if file != None]				#some files can be None depending on choice
	for path in file_list:
		if not os.path.exists(path):
			logger.info("\nError: {0} does not exist.".format(path))
			sys.exit(1)

	logger.info("Testing output directory/file writeability")
	make_directory(args.outdir)
	if not os.access(args.outdir, os.W_OK):
		logger.info("Error: {0} does not exist or is not writeable.".format(args.outdir))
		sys.exit(1)

	#Output write test
	for path in output_files[:-1]:	#Do not include log-file as this is managed by logger
		if path == None:
			continue
		if os.path.exists(path):
			if not os.access(path, os.W_OK):
				logger.info("Error: {0} could not be opened for writing.".format(path))
				sys.exit(1)

	#Open pdf for figures
	figure_pdf = PdfPages(figures_f, keep_empty=True)


	#----------------------------------------------------------------------------------------------------#
	# Read information in bam/fasta
	#----------------------------------------------------------------------------------------------------#

	logger.info("Reading info from .bam file")
	bamfile = pysam.AlignmentFile(args.bam, "rb")
	if bamfile.has_index() == False:
		logger.info("ERROR: No index found for bamfile")
		sys.exit()

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
			logger.critical("(Fastafile)\t{0} has length {1}".format(chrom, fasta_chrom_info[chrom]))
			logger.critical("(Bamfile)\t{0} has length {1}".format(chrom, bam_chrom_info[chrom]))
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
	regions_dict = {"input_regions":input_regions, "output_regions":output_regions, "peak_regions":peak_regions, "nonpeak_regions":nonpeak_regions}
	for sub in regions_dict:
		regions_sub = regions_dict[sub]
		regions_sub.subtract(blacklist_regions)
		regions_sub = regions_sub.apply_method(OneRegion.split_region, 50000)

		regions_sub.keep_chroms(chrom_in_common)
		regions_dict[sub] = regions_sub
	
	input_regions = regions_dict["input_regions"]
	output_regions = regions_dict["output_regions"]
	peak_regions = regions_dict["peak_regions"]
	nonpeak_regions = regions_dict["nonpeak_regions"]

	
	#write beds to look at in igv
	#input_regions.write_bed(os.path.join(args.outdir, "input_regions.bed"))
	#output_regions.write_bed(os.path.join(args.outdir, "output_regions.bed"))
	#peak_regions.write_bed(os.path.join(args.outdir, "peak_regions.bed"))
	#nonpeak_regions.write_bed(os.path.join(args.outdir, "nonpeak_regions.bed"))

	#Sort according to order in bam_references:
	output_regions.loc_sort(bam_references)
	chrom_order = {bam_references[i]:i for i in range(len(bam_references))}	 #for use later when sorting output 

	#### Statistics about regions ####
	genome_bp = sum([region.get_length() for region in genome_regions])
	blacklist_bp = sum([region.get_length() for region in blacklist_regions])
	peak_bp = sum([region.get_length() for region in peak_regions])
	nonpeak_bp = sum([region.get_length() for region in nonpeak_regions])
	input_bp = sum([region.get_length() for region in input_regions])
	output_bp = sum([region.get_length() for region in output_regions])

	logger.info("INFO\tGENOME\t{0} ({1:.2f}%)".format(genome_bp, genome_bp/genome_bp*100))
	logger.info("INFO\tBLACKLIST_REGIONS\t{0} ({1:.2f}%)".format(blacklist_bp, blacklist_bp/genome_bp*100))
	logger.info("INFO\tPEAK_REGIONS\t{0} ({1:.2f}%)".format(peak_bp, peak_bp/genome_bp*100))
	logger.info("INFO\tNONPEAK_REGIONS\t{0} ({1:.2f}%)".format(nonpeak_bp, nonpeak_bp/genome_bp*100))
	logger.info("INFO\tINPUT_REGIONS\t{0} ({1:.2f}%)".format(input_bp, input_bp/genome_bp*100))
	logger.info("INFO\tOUTPUT_REGIONS\t{0} ({1:.2f}%)".format(output_bp, output_bp/genome_bp*100))

	#----------------------------------------------------------------------------------------------------#
	# Estimate normalization factors
	#----------------------------------------------------------------------------------------------------#

	#Setup logger queue
	logger.debug("Setting up listener for log")
	log_q = mp.Manager().Queue()
	args.log_q = log_q
	listener = mp.Process(target=main_logger_process, args=(log_q, logger))
	listener.start()

	#----------------------------------------------------------------------------------------------------#

	logger.comment("")
	logger.critical("----- Estimating normalization factors -----")

	#If normalization is to be calculated
	if not args.norm_off:

		#Reads in peaks/nonpeaks
		logger.info("Counting reads in peak regions...")
		peak_region_chunks = peak_regions.chunks(args.split)
		reads_peaks = sum(run_parallel(count_reads, peak_region_chunks, [args], args.cores, logger))
		logger.comment("")

		logger.info("Counting reads in nonpeak regions...")
		nonpeak_region_chunks = nonpeak_regions.chunks(args.split)
		reads_nonpeaks = sum(run_parallel(count_reads, nonpeak_region_chunks, [args], args.cores, logger))

		total_reads = reads_peaks + reads_nonpeaks

		logger.comment("")
		logger.info("INFO\tTOTAL_READS\t{0}".format(total_reads))
		logger.info("INFO\tPEAK_READS\t{0}".format(reads_peaks))
		logger.info("INFO\tNONPEAK_READS\t{0}".format(reads_nonpeaks))

		lib_norm = 10000000/total_reads
		noise_norm = reads_nonpeaks/reads_peaks
		correct_factor = lib_norm*noise_norm

		logger.info("INFO\tLIB_NORM\t{0:.5f}".format(lib_norm))
		logger.info("INFO\tNOISE_NORM\t{0:.5f}".format(noise_norm))

	else:
		logger.info("Normalization was switched off")
		correct_factor = 1.0

	logger.info("INFO\tCORRECTION_FACTOR:\t{0:.5f}".format(correct_factor))

	#----------------------------------------------------------------------------------------------------#
	# Estimate sequence bias
	#----------------------------------------------------------------------------------------------------#

	logger.comment("")
	logger.critical("Started estimation of sequence bias...")

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
	logger.critical("Correcting reads from .bam within output regions")

	L = 2 * args.k_flank + 1

	#Getting bigwig files ready
	header = [(chrom, bam_chrom_info[chrom]) for chrom in bam_references]

	for track in output_bws:
		for strand in output_bws[track]:
			filename = output_bws[track][strand]["fn"]
			output_bws[track][strand]["pybw"] = pyBigWig.open(filename, "w")
			output_bws[track][strand]["pybw"].addHeader(header)

	pre_bias = {}
	post_bias = {}
	for direction in strands:
		pre_bias[direction] = SequenceMatrix.create(L, "PWM")
		post_bias[direction] = SequenceMatrix.create(L, "PWM")

	output_regions_chunks = output_regions.chunks(args.split)
	no_tasks = float(len(output_regions_chunks))

	#Start correction
	pool = mp.Pool(processes=args.cores)
	task_list = [pool.apply_async(bias_correction, args=[chunk, args, bias_obj]) for chunk in output_regions_chunks]
	pool.close()

	#Process results as they come in
	write_idx = 0		#index in task_list
	prev_progress = (-1, -1)
	while write_idx < len(task_list):
	
		if task_list[write_idx].ready() == True:

			result = task_list[write_idx].get()
			signals = result[0]
			pre_bias_chunk = result[1][0]
			post_bias_chunk = result[1][1]

			for direction in strands:
				pre_bias[direction].add_counts(pre_bias_chunk[direction])
				post_bias[direction].add_counts(post_bias_chunk[direction])

			#Write tracks to bigwig file
			for region in sorted(signals.keys(), key=lambda tup: (chrom_order[tup[0]], tup[1], tup[2])): #Ensure that positions are written to bigwig in correct order
				
				chrom, reg_start, reg_end = region
				positions = np.arange(reg_start, reg_end)	#genomic 0-based positions

				for track in tracks:	# only prints chosen tracks

					#Join signals if forward/reverse split
					if args.split_strands:
						signals[region][track]["both"] = signals[region][track]["forward"] + signals[region][track]["reverse"]	

					strands_found = signals[region][track]
					for strand in strands_found:

						signal = np.copy(signals[region][track][strand])  #Numpy array of signal

						signal[np.isclose(signal, 0)] = 0	#adjust for weird numpy floating point
						included = signal.nonzero()[0]
						pos = positions[included]
						val = signal[included]

						if len(pos) > 0:
							try:
								output_bws[track][strand]["pybw"].addEntries(chrom, pos, values=val, span=1)
							except:

								logger.info("Error writing region {0}.".format(region))
								print("TRACK: {0}".format(track))
								print("STRAND: {0}".format(strand))
								print("SIGNAL: {0}".format(signals[region][track][strand]))
								sys.exit()

			#Clear memory
			del result
			del signal
			gc.collect()

			#Wait for next chunk to print
			write_idx += 1 #This is also the same as the number of prints done 

		tasks_done = sum([task.ready() for task in task_list])
		if tasks_done != prev_progress[0] or write_idx != prev_progress[1]:
			logger.info("Correction progress: {0:.0f}% | Writing progress: {1:.0f}%".format(tasks_done/no_tasks*100, write_idx/no_tasks*100))
		
		prev_progress = (tasks_done, write_idx)

	#Done computing
	pool.join()
	gc.collect()

	#Done writing to files
	logger.info("Closing bigwig files (this can take a while)")
	for track in output_bws:
		for strand in output_bws[track]:
			output_bws[track][strand]["pybw"].close()



	#---------------------------------------------------#
	
	logger.debug("Waiting for listener to finish")
	log_q.put(None)
	while listener.exitcode != 0:
		logger.debug("Listener exitcode is: {0}".format(listener.exitcode))
		time.sleep(1)

	logger.debug("Joining listener")
	listener.join()
	


	#----------------------------------------------------------------------------------------------------#
	# Information and verification of corrected read frequencies
	#----------------------------------------------------------------------------------------------------#		

	logger.comment("")
	logger.info("Verifying bias correction")

	#Calculating variance per base
	for strand in strands:

		#Join negative/positive counts from post-correction bias
		abssum = np.abs(np.sum(post_bias[strand].neg_counts, axis=0))
		post_bias[strand].neg_counts = abssum + post_bias[strand].neg_counts
		post_bias[strand].counts += post_bias[strand].neg_counts	#now pos

		pre_bias[strand].prepare_mat()
		post_bias[strand].prepare_mat()

		pre_var = np.mean(np.var(pre_bias[strand].bias_pwm, axis=1)[:4])   #mean of variance per nucleotide
		post_var = np.mean(np.var(post_bias[strand].bias_pwm, axis=1)[:4])
		logger.info("INFO\tpre-bias variance {0}:\t{1:.7f}".format(strand, pre_var))
		logger.info("INFO\tpost-bias variance {0}:\t{1:.7f}".format(strand, post_var))

		#Plot figure
		fig_title = "Nucleotide frequencies in corrected reads\n({0} strand)".format(strand)
		figure_pdf.savefig(plot_correction(pre_bias[strand].bias_pwm, post_bias[strand].bias_pwm, fig_title))

	
	#----------------------------------------------------------------------------------------------------#
	# Finish up
	#----------------------------------------------------------------------------------------------------#

	figure_pdf.close()
	end_time = datetime.now()
	logger.comment("")
	logger.critical("Finished ATACorrect run (total time of {0})".format(end_time - begin_time))


#--------------------------------------------------------------------------------------------------------#
if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser = add_atacorrect_arguments(parser)
	args = parser.parse_args()

	if len(sys.argv[1:]) == 0:
		parser.print_help()
		sys.exit()

	run_atacorrect(args)
