#!/usr/bin/env python

"""
Uilities for argparse, logging, file handling, multiprocessing in TOBIAS

@author: Mette Bentsen
@contact: mette.bentsen (at) mpi-bn.mpg.de
@license: MIT
	
"""

import os
import logging
import sys
import argparse
import collections
import time
import textwrap
import string
import numpy as np
import copy
from difflib import SequenceMatcher

import pyBigWig
from tobias.utils.logger import *


#-------------------------------------------------------------------------------------------#
#----------------------------------- Multiprocessing ---------------------------------------#
#-------------------------------------------------------------------------------------------#
def check_cores(given_cores, logger):
	""" Checks number of available cores and sets the used cores to <= available cores """

	available_cores = mp.cpu_count()
	if given_cores > available_cores:
		logger.warning("Number of available cores is {0} but \'--cores\' is set to {1}. Setting \'--cores\' to {0}.\n".format(available_cores, given_cores))
		return(available_cores)
	else:
		return(given_cores)


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


def file_writer(q, key_file_dict, args):
	""" File-writer per key -> to value file """

	#Open handles for all files (but only once per file!)
	file2handle = {}
	for fil in set(key_file_dict.values()):
		try:
			file2handle[fil] = open(fil, "w")
		except Exception:
			print("Tried opening file {0} in file_writer but something went wrong?".format(fil))
			traceback.print_exc(file=sys.stderr)

	#Assign handles to keys
	handles = {}
	for key in key_file_dict:
		handles[key] = file2handle[key_file_dict[key]]

	#Fetching string content from queue
	while True:
		try:
			(key, content) = q.get()
			if key == None:
				break

			handles[key].write(content)

		except Exception:
			import sys, traceback
			print('Problem in file_writer:', file=sys.stderr)
			traceback.print_exc(file=sys.stderr)
			break

	#Got all regions in queue, close files
	for fil, handle in file2handle.items():
		handle.close()
		
	return(1)	


def bigwig_writer(q, key_file_dict, header, regions, args):
	""" Handle queue to write bigwig, args contain extra info such as verbosity and log_q """

	#todo: check if args.log_q exists
	logger = TobiasLogger("", args.verbosity, args.log_q)	#separate core, needs mp logger through queue
	logger.debug("Opened bigwig writer process for {0}".format(key_file_dict))
	logger.debug("Header: {0}".format(header))

	handles = {}
	for key in key_file_dict:
		logger.debug("Opening file {0} for writing".format(key_file_dict[key]))	
		try:
			handles[key] = pyBigWig.open(key_file_dict[key], "w")
			handles[key].addHeader(header)

		except Exception:
			print("Tried opening file {0} in bigwig_writer but something went wrong?".format(fil))
			traceback.print_exc(file=sys.stderr)

	#Correct order of chromosomes as given in header
	contig_list = [tup[0] for tup in header]
	order_dict = dict(zip(contig_list, range(len(contig_list))))
	
	#Establish order of regions to be writteninput regions
	region_tups = [(region.chrom, region.start, region.end) for region in regions]
	sorted_region_tups = sorted(region_tups, key=lambda tup: (order_dict[tup[0]], tup[1]))			#sort to same order as bigwig header
	no_regions = len(region_tups)

	#Fetching content from queue
	logger.debug("Fetching content from queue")

	i_to_write = {key:0 for key in handles}			#index of next region to write
	ready_to_write = {key:{} for key in handles}	#key:dict; dict is region-tup:signal array 
	while True:

		try:
			(key, region, signal) = q.get()	#key is the bigwig key (e.g. bias:forward), region is a tup of (chr, start, end)
			logger.spam("Received signal {0} for region {1}".format(key, region))
			
			if key == None:	#none is only inserted once all regions have been sent
				for akey in i_to_write:
					if i_to_write[akey] != no_regions - 1:
						logger.error("Wrote {0} regions but there are {1} in total".format(i_to_write[akey], len(region_tups)))
						logger.error("Ready_to_write[{0}]: {1}".format(akey, len(ready_to_write[akey])))
				break

			#Save key:region:signal to ready_to_write
			ready_to_write[key][region] = signal
			
			writing_progress = Progress(no_regions, logger, prefix="Writing progress", round=0)

			#Check if next-to-write region was done
			for key in handles: 

				#Only deal with writing if there are still regions to write for this handle
				if i_to_write[key] != no_regions - 1:
					next_region = sorted_region_tups[i_to_write[key]]	#this is the region to be written next for this key
				
					#If results are in; write wanted entry to bigwig
					while next_region in ready_to_write[key]: 	#When true: Keep writing when the next region is available
						chrom = next_region[0]
						signal = ready_to_write[key][next_region]
						included = signal.nonzero()[0]
						positions = np.arange(next_region[1],next_region[2])		#start-end	(including end)
						pos = positions[included]
						val = signal[included]

						if len(pos) > 0:
							try:
								handles[key].addEntries(chrom, pos, values=val, span=1)
							except:
								logger.error("Error writing key: {0}, region: {1} to bigwig".format(key, next_region))
						logger.spam("Wrote signal {0} from region {1}".format(key, next_region))

						#Check whether this was the last region
						if i_to_write[key] == no_regions - 1: #i_to_write is the last idx in regions; all sites were written
							logger.info("Closing {0} (this might take some time)".format(key_file_dict[key]))
							handles[key].close()
							next_region = None 	#exit the while loop
						else:
							i_to_write[key] += 1
							next_region = sorted_region_tups[i_to_write[key]]	#this is the region to be written next for this key

						#Writing progress
						#progress = sum([i_to_write[key] for key in handles])
						#writing_progress.write(progress)

		except Exception:
			import sys, traceback
			print('Problem in file_writer:', file=sys.stderr)
			traceback.print_exc(file=sys.stderr)
			break

	return(1)


def monitor_progress(task_list, logger, prefix="Progress"):

	prev_done = 0
	no_tasks = len(task_list)
	done = sum([task.ready() for task in task_list])
	logger.info("{0} 0%".format(prefix))
	while done != no_tasks:
		done = sum([task.ready() for task in task_list])
		if done != prev_done:
			#print if done is 
			
			logger.info("{1} {0}%".format(round(done/no_tasks*100.0, max(0,len(str(no_tasks))-1)), prefix))
			prev_done = done
		else:
			time.sleep(0.1)
	
	logger.info("{0} done!".format(prefix))

	return() 	#doesn't return until the while loop exits


#-------------------------------------------------------------------------------------------#
#------------------------------------- Argparser -------------------------------------------#
#-------------------------------------------------------------------------------------------#


def restricted_float(f, f_min, f_max):
    f = float(f)
    if f < f_min or f > f_max:
        raise argparse.ArgumentTypeError("{0} not in range [0.0, 1.0]".format(f))
    return f

def restricted_int(integer, i_min, i_max):
	integer = float(integer)
	if integer < i_min or integer > i_max:
		raise

def format_help_description(name, description, width=90):
	""" Format description of command line tool --help description """

	formatted = "" #initialize

	#Calculate needed whitespace in comparison to header
	header = "TOBIAS ~ {0}".format(name)
	ws = int((width - len(header))/2.0)

	formatted += "_"*width + "\n"*2
	formatted += "{0}{1}{0}\n".format(" "*ws, header)
	formatted += "_"*width + "\n"*2

	for segment in description.split("\n"):
		formatted += "\n".join(textwrap.wrap(segment, width)) + "\n"	#Split description on space 

	if description != "":
		formatted += "\n" + "-"*width + "\n"

	return(formatted)


def check_required(args, required):
	""" Checks required keys in input args """

	for arg in required:
		if getattr(args, arg) == None:
			sys.exit("ERROR: Missing argument --{0}".format(arg))

#-------------------------------------------------------------------------------------------#
#---------------------------------------- Misc ---------------------------------------------#
#-------------------------------------------------------------------------------------------#


class Progress:
	""" Class for writing out progress of processes such as multiprocessing """
	def __init__(self, total_elements, logger, prefix="Progress", round=0):
		
		self.total_elements = total_elements
		self.logger = logger
		self.prefix = prefix
		self.round = round
		self.last_written = -1

	def write(self, progress):
		""" Write out progress if it was not already written """

		percent_to_write = progress / self.total_elements * 100
		if self.round == 0:
			percent_to_write = round(percent_to_write)
		else:
			percent_to_write = round(percent_to_write, self.round)

		#Only write if this level has not already been written
		if percent_to_write != self.last_written:
			self.logger.info("{0}: {1}%".format(self.prefix, percent_to_write))
			self.last_written = percent_to_write


def flatten_list(lst):

    for element in lst:
        if isinstance(element, collections.Iterable) and not isinstance(element, (str, bytes)):
            yield from flatten_list(element)
        else:
            yield element

def check_files(lst_of_files, action="r"):

	flat_lst = flatten_list(lst_of_files)
	for fil in flat_lst:
		if fil != None:
			if action == "r":
				if os.path.exists(fil):
					try: 
						with open(fil) as f:
							pass
					except:
						sys.exit("ERROR: Could not open file \"{0}\" for reading".format(fil))
				else:
					sys.exit("ERROR: File \"{0}\" does not exists".format(fil))

			elif action == "w":
				if os.path.exists(fil):
					try: 
						with open(fil, "w") as f:
							pass
					except:
						sys.exit("ERROR: Could not open file \"{0}\" for writing. Please check that you do not have the file open and that you have write permission.".format(fil))
				
def make_directory(directory):
	if not os.path.isfile(directory) and not os.path.exists(directory):
		os.makedirs(directory)
	

def merge_dicts(dicts):
	""" Merge recursive keys and values for list of dicts into one dict. Values are added numerically / lists are extended / numpy arrays are added"""

	def merger(dct, dct_to_add):
		""" Merge recursive keys and values of dct_to_add into dct. Values are added numerically / lists are extended / numpy arrays are added
		No return - dct is changed in place """

		for k, v in dct_to_add.items():

				#If k is a dict, go one level down
				if (k in dct and isinstance(dct[k], dict)):
					merger(dct[k], dct_to_add[k])
				else:
					if not k in dct:
						dct[k] = dct_to_add[k]
					else:
						dct[k] += dct_to_add[k]
	
	#Initialize with the first dict in list
	out_dict = copy.deepcopy(dicts[0])
	for dct in dicts[1:]:
		merger(out_dict, dct)

	return(out_dict)

def filafy(astring):
	""" Make string into accepted filename """

	valid_chars = "-_.%s%s" % (string.ascii_letters, string.digits)
	filename = ''.join(char for char in astring if char in valid_chars)
	return(filename)


def get_closest(value, arr): 
	"""Find element the element in arr which is closest to value """

	idx = (np.abs(arr-value)).argmin()
	return(arr[idx])


#-------------------------------------------------------------------------------------------#
#-------------------------------------- Matching -------------------------------------------#
#-------------------------------------------------------------------------------------------#

def common_prefix(strings):
	""" Find the longest string that is a prefix of all the strings. Used in PlotChanges to find TF names from list """
	if not strings:
		return ''
	prefix = strings[0]
	for s in strings:
		if len(s) < len(prefix):
			prefix = prefix[:len(s)]
		if not prefix:
			return ''
		for i in range(len(prefix)):
			if prefix[i] != s[i]:
				prefix = prefix[:i]
				break
	return prefix


def match_lists(lofl): # list of lists
	""" Find matches between list1 and list2 (output will be the length of list1 with one or more matches per element)."""

	#Remove common prefixes/suffixes per list
	prefixes = []
	suffixes = []
	for i, lst in enumerate(lofl):
		if len(lst) > 1:	#make sure to only compare between more than one element (otherwise it sets the whole element as prefix)
			prefix = common_prefix(lst)
			suffix = common_prefix([element[::-1] for element in lst])[::-1]
		else:
			prefix = ""
			suffix = ""

		lofl[i] = [element.replace(prefix, "").replace(suffix, "") for element in lst]
		prefixes.append(prefix)
		suffixes.append(suffix)

	#Initialize
	matches = [[[] for _ in lofl[0]] for _ in range(len(lofl)-1)]	#list of lists (len(lofl) * len(lofl[0]) * list for collecting results)

	#Find best matches to each element in col:
	for col in range(len(lofl)-1):

		lst1 = lofl[0]  	#[row for row in lofl[0]]	#Always match to first list 
		lst2 = lofl[col+1] 	# [row for row in lofl[col+1]] 	#compare to column

		for i, element1 in enumerate(lst1):

			local_match_scores = []
			global_match_scores = []
			for element2 in lst2:
				match_tups = SequenceMatcher(None, element1.lower(), element2.lower()).get_matching_blocks()
				match_sum = sum([tup.size for tup in match_tups if tup.size > 1])

				local_match_scores.append(match_sum / float(len(element1)))		#local score shows how well element1 fits locally in element2
				global_match_scores.append(match_sum / float(len(element1) + len(element2)))	#takes into account length of both element1 and 2

			#Best local match scores
			best_score = max(local_match_scores)
			best_idx = [idx for idx, score in enumerate(local_match_scores) if score == best_score]

			#Sort best_idx by global match scores
			best_global_match_scores = [global_match_scores[idx] for idx in best_idx]
			best_idx = [idx for _,idx in sorted(zip(best_global_match_scores, best_idx), key=lambda pair: pair[0], reverse=True)]	#largest global match scores first
			best_elements = [lst2[idx] for idx in best_idx]		#elements in lst2 that fit to element1 
			matches[col][i].extend(best_elements)

	#Map back to full names
	for col in range(1, len(lofl)): 	# number of lists to get matches from
		for row in range(len(lofl[0])): 	#elements in lofl[0]
			matches[col-1][row] = [prefixes[col] + match + suffixes[col] for match in matches[col-1][row]]

	return(matches)
