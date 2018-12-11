#!/usr/bin/env python

"""
Classes and functions for performing bias estimation, correction and visualization in ATACorrect

@author: Mette Bentsen
@contact: mette.bentsen (at) mpi-bn.mpg.de
@license: MIT

"""

import os
import sys
import gc
import numpy as np
import multiprocessing as mp
import time
from datetime import datetime

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.optimize import curve_fit

#Bio-specific packages
import pysam

#Internal functions and classes
from tobias.utils.sequences import *
from tobias.utils.signals import *
from tobias.utils.ngs import *
from tobias.utils.utilities import * 

#Catch warnings from curve_fit
import warnings
from scipy.optimize import OptimizeWarning
warnings.simplefilter("error", OptimizeWarning)

#--------------------------------------------------------------------------------------------------#
class AtacBias:
	""" Class for storing information about estimated bias """

	def __init__(self, L, stype):

		self.stype = stype
		self.bias = {"forward": SequenceMatrix.create(L, self.stype),
					 "reverse": SequenceMatrix.create(L, self.stype),
					 "both": SequenceMatrix.create(L, self.stype)}
		self.no_reads = 0

	def join(self, obj):
		""" Join counts from AtacBias obj with other AtacBias obj """

		self.bias["forward"].add_counts(obj.bias["forward"])
		self.bias["reverse"].add_counts(obj.bias["reverse"])
		self.bias["both"].add_counts(obj.bias["both"])
		self.no_reads += obj.no_reads


####################################################################################################
####################################################################################################

def run_parallel(FUNC, input_chunks, arguments, n_cores, logger):
	"""
	#FUNC is the function to run
	#input_chunks is the input to loop over
	#arguments are arguments to func
	#logger is a Logging.Logger object
	"""

	no_chunks = len(input_chunks)

	if n_cores > 1:

		#Send jobs to pool
		pool = mp.Pool(processes=n_cores)
		task_list = []
		for input_chunk in input_chunks:
			task_list.append(pool.apply_async(FUNC, args=[input_chunk] + arguments))
		pool.close() 	#done sending jobs to pool

		#Wait for tasks to finish
		count = -1
		finished = sum([task.ready() for task in task_list])
		while finished < no_chunks:
			finished = sum([task.ready() for task in task_list])
			if count != finished:
				logger.info("Progress: {0:.0f}%".format(finished/float(no_chunks)*100))
				count = finished
			else:
				time.sleep(0.5)
		pool.join()

		#Get results from processes
		output_list = [task.get() for task in task_list]

	else:

		output_list = []
		for count, input_chunk in enumerate(input_chunks):
			logger.info("Progress: {0:.0f}%".format(count/float(no_chunks)*100))
			output_list.append(FUNC(input_chunk, *arguments))
	
	return(output_list)

#--------------------------------------------------------------------------------------------------#
def count_reads(regions_list, params):
	""" Count reads from bam within regions (counts position of cutsite to prevent double-counting) """

	bam_f = params.bam
	read_shift = params.read_shift
	bam_obj = pysam.AlignmentFile(bam_f, "rb")

	log_q = params.log_q
	logger = create_mp_logger(params.verbosity, log_q)	#sending all logger calls to log_q

	#Count per region
	read_count = 0
	logger.debug("Started counting region_chunk ({0} -> {1})".format("_".join([str(element) for element in regions_list[0]]), "_".join([str(element) for element in regions_list[-1]])))
	for region in regions_list:
		read_lst = ReadList().from_bam(bam_obj, region)
		
		for read in read_lst:  
			read.get_cutsite(read_shift)
			if read.cutsite > region.start and read.cutsite < region.end:  #only reads within borders
				read_count += 1
				
	logger.debug("Finished counting region_chunk ({0} -> {1})".format("_".join([str(element) for element in regions_list[0]]), "_".join([str(element) for element in regions_list[-1]])))
	bam_obj.close()

	return(read_count)


#--------------------------------------------------------------------------------------------------#
def bias_estimation(regions_list, params):
	""" Estimates bias of insertions within regions """

	#Info on run
	bam_f = params.bam
	fasta_f = params.genome
	k_flank = params.k_flank
	bg_shift = params.bg_shift
	read_shift = params.read_shift
	L = 2 * k_flank + 1
	log_q = params.log_q
	logger = create_mp_logger(params.verbosity, log_q)	#sending all logger calls to log_q

	#Open objects for reading
	bam_obj = pysam.AlignmentFile(bam_f, "rb")
	fasta_obj = pysam.FastaFile(fasta_f)
	chrom_lengths = dict(zip(bam_obj.references, bam_obj.lengths))  #Chromosome boundaries from bam_obj

	bias_obj = AtacBias(L, params.score_mat)

	strands = ["forward", "reverse"] if params.split_strands else ["both"] 

	#Estimate bias at each region
	for region in regions_list:

		read_lst = ReadList().from_bam(bam_obj, region)  #fragment_lst.to_read_list() 				#Fragments to single read list
		for read in read_lst:
			read.get_cutsite(read_shift)

		## Kmer cutting bias ##
		if len(read_lst) > 0:

			#Extract sequence
			extended_region = region.extend_reg(k_flank + bg_shift) 				#Extend to allow full kmers
			extended_region.check_boundary(chrom_lengths, "cut")
			sequence_obj = GenomicSequence(extended_region).from_fasta(fasta_obj)
			
			#Split reads forward/reverse
			for_lst, rev_lst = read_lst.split_strands()
			read_lst_strand = {"both": read_lst, "forward": for_lst, "reverse": rev_lst} 

			for strand in strands:

				#Map reads to positions
				read_per_pos = {}
				for read in read_lst_strand[strand]:
					read_per_pos[read.cutsite] = read_per_pos.get(read.cutsite, []) + [read]

				#Set all reads to forward if "both"
				if strand == "both":
					for cutsite in read_per_pos:
						read_per_pos[cutsite][0].is_reverse = False

				for cutsite in read_per_pos:
					if cutsite > region.start and cutsite < region.end:  	#only reads within borders
						read = read_per_pos[cutsite][0] 					#use first read in list to establish kmer
						no_cut = len(read_per_pos[cutsite]) 	

						read.get_kmer(sequence_obj, k_flank)

						bias_obj.bias[strand].add_sequence(read.kmer, no_cut)
						read.shift_cutsite(bg_shift)
						read.get_kmer(sequence_obj, k_flank)
						bias_obj.bias[strand].add_background(read.kmer, no_cut)

						bias_obj.no_reads += no_cut
	bam_obj.close()
	fasta_obj.close()

	return(bias_obj) 	#object containing information collected on bias 	


#--------------------------------------------------------------------------------------------------#
def relu(x, a, b):
	y = np.maximum(0.0, a*x + b)
	return(y)

def x_thresh(a,b):
	x = (0-b)/a 	#x position of y=0 intersect 
	return(x)

#--------------------------------------------------------------------------------------------------#
def bias_correction(regions_list, params, bias_obj):
	""" Corrects bias in cutsites (from bamfile) using estimated bias """

	bam_f = params.bam
	fasta_f = params.genome
	k_flank = params.k_flank
	read_shift = params.read_shift
	L = 2 * k_flank + 1
	w = params.window
	f = int(w/2.0)

	f_extend = k_flank + f

	strands = ["forward", "reverse"] if params.split_strands else ["both"] 

	pre_bias = {strand: SequenceMatrix.create(L, "PWM") for strand in strands}
	post_bias = {strand: SequenceMatrix.create(L, "PWM") for strand in strands}

	#Open bamfile and fasta
	bam_obj = pysam.AlignmentFile(bam_f, "rb")
	fasta_obj = pysam.FastaFile(fasta_f)
	chrom_lengths = dict(zip(bam_obj.references, bam_obj.lengths))  #Chromosome boundaries from bam_obj
	
	out_signals = {}
	
	#Go through each region
	for region_obj in regions_list:

		region_obj.extend_reg(f_extend)
		region_obj.check_boundary(chrom_lengths, "cut")
		reg_len = region_obj.get_length()	#length including flanking
		reg_key = (region_obj.chrom, region_obj.start+f_extend, region_obj.end-f_extend)	#output region
		out_signals[reg_key] = {"uncorrected":{}, "bias":{}, "expected":{}, "corrected":{}}

		################################
		####### Uncorrected reads ######
		################################

		#Get cutsite positions for each read
		read_lst = ReadList().from_bam(bam_obj, region_obj) 
		for read in read_lst:
			read.get_cutsite(read_shift)

		#Exclude reads with cutsites outside region
		read_lst = ReadList([read for read in read_lst if read.cutsite > region_obj.start and read.cutsite < region_obj.end])
		for_lst, rev_lst = read_lst.split_strands()
		read_lst_strand = {"both":read_lst, "forward": for_lst, "reverse": rev_lst}

		for strand in strands:
			out_signals[reg_key]["uncorrected"][strand] = read_lst_strand[strand].signal(region_obj) 
			out_signals[reg_key]["uncorrected"][strand] *= bias_obj.correction_factor 						#Correction factor
			#out_signals[reg_key]["uncorrected"][strand] = np.log2(out_signals[reg_key]["uncorrected"][strand] + 1)	

		
		del read_lst_strand
		del for_lst

		################################
		###### Estimation of bias ######
		################################

		#Get sequence in this region
		sequence_obj = GenomicSequence(region_obj).from_fasta(fasta_obj)
		
		#Score sequence using forward/reverse motifs
		for strand in strands:
			if strand == "forward" or strand == "both":
				seq = sequence_obj.sequence
				bias = bias_obj.bias[strand].score_sequence(seq)
			elif strand == "reverse":
				seq = sequence_obj.revcomp
				bias = bias_obj.bias[strand].score_sequence(seq)[::-1]  #3'-5'

			out_signals[reg_key]["bias"][strand] = np.nan_to_num(bias) 		#convert any nans to 0
		

		#################################
		###### Correction of reads ######
		#################################
		
		reg_end = reg_len - k_flank
		step = 10
		overlaps = int(params.window / step)
		window_starts = list(range(k_flank, reg_end-params.window, step))
		window_ends = list(range(k_flank+params.window, reg_end, step))
		window_ends[-1] = reg_len
		windows = list(zip(window_starts, window_ends))

		for strand in strands:

			########### Estimate bias threshold ###########
			bias_predictions = np.zeros((overlaps,reg_len))
			row = 0

			for window in windows:

				signal_w = out_signals[reg_key]["uncorrected"][strand][window[0]:window[1]]
				bias_w = out_signals[reg_key]["bias"][strand][window[0]:window[1]]

				signalmax = np.max(signal_w)
				biasmin = np.min(bias_w) 
				biasmax = np.max(bias_w)

				if signalmax > 0:
					try:
						popt, pcov = curve_fit(relu, bias_w, signal_w)
						bias_predict = relu(bias_w, *popt)

					except (OptimizeWarning, RuntimeError):
						cut_positions = np.logical_not(np.isclose(signal_w, 0))
						bias_min = np.min(bias_w[cut_positions])
						bias_predict = bias_w - bias_min
						bias_predict[bias_predict < 0] = 0

					if np.max(bias_predict) > 0:
						bias_predict = bias_predict / np.max(bias_predict)
				else:
					bias_predict = np.zeros(window[1]-window[0])	

				bias_predictions[row, window[0]:window[1]] = bias_predict
				row += 1 if row < overlaps - 1 else 0

			bias_prediction = np.mean(bias_predictions, axis=0)
			bias = bias_prediction 

			######## Calculate expected signal ######
			signal_sum = fast_rolling_math(out_signals[reg_key]["uncorrected"][strand], w, "sum")
			signal_sum[np.isnan(signal_sum)] = 0 	#f-width ends of region

			bias_sum = fast_rolling_math(bias, w, "sum")	#ends of arr are nan
			nulls = np.logical_or(np.isclose(bias_sum, 0), np.isnan(bias_sum))
			bias_sum[nulls] = 1 		# N-regions will give stretches of 0-bias
			bias_probas = bias / bias_sum
			bias_probas[nulls] = 0 		#nan to 0

			out_signals[reg_key]["expected"][strand] = signal_sum * bias_probas 

			######## Correct signal ########
			out_signals[reg_key]["corrected"][strand] = out_signals[reg_key]["uncorrected"][strand] - out_signals[reg_key]["expected"][strand]

			
			##### Rescale back to original sum #####
			pos_corrected = np.abs(out_signals[reg_key]["corrected"][strand])
			corrected_sum = fast_rolling_math(pos_corrected, w, "sum")
			corrected_sum[np.isclose(corrected_sum, 0)] = np.nan

			frac = signal_sum / corrected_sum  	#If corrected sum was smaller than signal -> multiply
			frac[np.isnan(frac)] = 1

			out_signals[reg_key]["corrected"][strand] = out_signals[reg_key]["corrected"][strand] * frac


		#######################################
		########   Verify correction   ########
		#######################################
		
		#Verify correction across all reads
		for direction in strands:
			for idx in range(k_flank,reg_len - k_flank -1): 
				if idx > k_flank and idx < reg_len-k_flank:

					orig = out_signals[reg_key]["uncorrected"][direction][idx]
					correct = out_signals[reg_key]["corrected"][direction][idx]

					if direction == "forward" or direction == "both":
						kmer = sequence_obj.sequence[idx-k_flank:idx+k_flank+1]
					else:
						kmer = sequence_obj.revcomp[reg_len-idx-k_flank-1:reg_len-idx+k_flank]

					#Save kmer for bias correction verification
					pre_bias[direction].add_sequence(kmer, orig)
					post_bias[direction].add_sequence(kmer, correct)

		#Set size back to original
		for track in out_signals[reg_key]:
			for direction in out_signals[reg_key][track]:
				out_signals[reg_key][track][direction] = out_signals[reg_key][track][direction][f_extend:-f_extend]

	bam_obj.close()
	fasta_obj.close()
	
	gc.collect()

	return(out_signals, [pre_bias, post_bias])



####################################################################################################
######################################## Plot functions ############################################
####################################################################################################

colors = {0:"green", 1:"red", 2:"blue", 3:"darkkhaki"}
names = {0:"A", 1:"T", 2:"C", 3:"G"}

def plot_pssm(matrix, title):
	""" Plot pssm in matrix """

	#Make figure
	fig, ax = plt.subplots()
	fig.suptitle(title, fontsize=16, weight="bold")

	#Formatting of x axis
	length = matrix.shape[1]
	flank = int(length/2.0)		
	xvals = np.arange(length)  # each position corresponds to i in mat
	
	#Customize minor tick labels
	xtick_pos = xvals[:-1] + 0.5
	xtick_labels = list(range(-flank, flank))
	ax.xaxis.set_major_locator(ticker.FixedLocator(xvals))
	ax.xaxis.set_major_formatter(ticker.FixedFormatter(xtick_labels))	
	ax.xaxis.set_minor_locator(ticker.FixedLocator(xtick_pos))			#locate minor ticks between major ones (cutsites)
	ax.xaxis.set_minor_formatter(ticker.NullFormatter())

	#Make background grid on major ticks
	plt.grid(color='0.8', which="minor", ls="--", axis="x")

	plt.xlim([0, length-1])
	plt.xlabel('Position from cutsite')
	plt.ylabel('PSSM score')

	######## Plot data #######
	#Plot PSSM / bias motif
	for nuc in range(4):
		plt.plot(xvals, matrix[nuc,:], color=colors[nuc], label=names[nuc])

	#Cutsite-line
	plt.axvline(flank-0.5, linewidth=2, color="black", zorder=100)

	#Finish up
	plt.legend(loc="lower right")
	plt.tight_layout()
	fig.subplots_adjust(top=0.88, hspace=0.5)

	return(fig)

#----------------------------------------------------------------------------------------------------#
def plot_correction(pre_mat, post_mat, title):
	""" Plot comparison of pre-correction and post-correction matrices """

	#Make figure
	fig, ax = plt.subplots()
	fig.suptitle(title, fontsize=16, weight="bold")

	L = pre_mat.shape[1]
	flank = int(L/2.0)		
	xvals = np.arange(L)  # each position corresponds to i in mat
	
	#Customize minor tick labels
	xtick_pos = xvals[:-1] + 0.5
	xtick_labels = list(range(-flank, flank))		#-flank - flank without 0
	ax.xaxis.set_major_locator(ticker.FixedLocator(xvals))
	ax.xaxis.set_major_formatter(ticker.FixedFormatter(xtick_labels))	
	ax.xaxis.set_minor_locator(ticker.FixedLocator(xtick_pos))			#locate minor ticks between major ones (cutsites)
	ax.xaxis.set_minor_formatter(ticker.NullFormatter())
	
	#PWMs for all mats
	pre_pwm = pre_mat
	post_pwm = post_mat

	#Pre correction
	for nuc in range(4):
		yvals = [pre_pwm[nuc, m] for m in range(L)]
		plt.plot(xvals, yvals, linestyle="--", color=colors[nuc], linewidth=1, alpha=0.5)
	
	#Post correction
	for nuc in range(4):
		yvals = [post_pwm[nuc, m] for m in range(L)]
		plt.plot(xvals, yvals, color=colors[nuc], linewidth=2, label=names[nuc])
	
	plt.xlim([0, L-1])
	plt.xlabel('Position from cutsite')
	plt.ylabel('Nucleotide frequency')

	#Set legend
	plt.plot([0],[0], linestyle="--", linewidth=1, color="black", label="pre-correction")
	plt.plot([0],[0], color="black", label="post-correction")
	plt.legend(loc="lower right", prop={'size':6})

	plt.tight_layout()
	fig.subplots_adjust(top=0.88, hspace=0.4)

	return(fig)
