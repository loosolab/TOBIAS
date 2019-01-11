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
from difflib import SequenceMatcher


#-------------------------------------------------------------------------------------------#
#----------------------------------- Logger stuff ------------------------------------------#
#-------------------------------------------------------------------------------------------#

class TOBIASFormatter(logging.Formatter):

	default_fmt = logging.Formatter("%(asctime)s (%(processName)s)\t%(message)s", "%Y-%m-%d %H:%M:%S")
	comment_fmt = logging.Formatter("%(message)s")

	#def __init__(self, fmt="(%(processName)s) %(asctime)s\t\t%(message)s"):
	#	logging.Formatter.__init__(self, fmt)

	def format(self, record):
		format_orig = self._fmt

		#Comments
		if record.levelno > logging.CRITICAL:
			return self.comment_fmt.format(record)
		else:
			return self.default_fmt.format(record)	


def create_logger(verbosity=2, log_f=None):
	
	#verbose_levels = {1:logging.CRITICAL, 2:logging.ERROR, 3:logging.WARNING, 4:logging.INFO, 5:logging.DEBUG}
	verbose_levels = {1:logging.CRITICAL, 2:logging.INFO, 3:logging.DEBUG}
	logger = logging.getLogger("TOBIAS")	#Set logger specifically to not get other module logs (which write to root)
	logger.setLevel(verbose_levels[verbosity])

	formatter = TOBIASFormatter()

	#Log file stream
	if log_f != None:
		log = logging.FileHandler(log_f, "w")
		log.setLevel(verbose_levels[verbosity])
		log.setFormatter(formatter)
		logger.addHandler(log)

	#Stdout stream
	else:
		con = logging.StreamHandler(sys.stdout)		#console output
		con.setLevel(verbose_levels[verbosity])
		con.setFormatter(formatter)
		logger.addHandler(con)


	#Create custom level for comments (always shown)
	comment_level = logging.CRITICAL+10
	logging.addLevelName(comment_level, "comment") #log_levels[lvl])
	setattr(logger, 'comment', lambda *args: logger.log(comment_level, *args))

	return(logger)


def create_mp_logger(verbosity, queue):

	verbose_levels = {1:logging.CRITICAL, 2:logging.INFO, 3:logging.DEBUG}
	
	h = logging.handlers.QueueHandler(queue)  	# Just the one handler needed
	logger = logging.getLogger("Worker")
	logger.handlers = []
	logger.addHandler(h)
	logger.setLevel(verbose_levels[verbosity])

	return(logger)


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


#-------------------------------------------------------------------------------------------#
#----------------------------------- Multiprocessing ---------------------------------------#
#-------------------------------------------------------------------------------------------#

def file_writer(q, key_file_dict, args):
	""" File-writer per key -> to value file """

	#time_spent_writing = 0
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

			#begin_time = time.time()
			#print("writing content of {0} chars for TF {1}".format(len(content), TF))
			handles[key].write(content)
			#end_time = time.time()
			#time_spent_writing += end_time - begin_time
			#except Queue.Empty:

		except Exception:
			import sys, traceback
			print('Problem in file_writer:', file=sys.stderr)
			traceback.print_exc(file=sys.stderr)
			break

	#Got all regions in queue, close files
	for fil, handle in file2handle.items():
		handle.close()
		
	return(1)	


def monitor_progress(task_list, logger):

	prev_done = 0
	no_tasks = len(task_list)
	done = sum([task.ready() for task in task_list])
	logger.info("Progress 0%")
	while done != no_tasks:
		done = sum([task.ready() for task in task_list])
		if done != prev_done:
			logger.info("Progress {0}%".format(round(done/no_tasks*100.0, max(0,len(str(no_tasks))-1))))
			prev_done = done
		else:
			time.sleep(0.1)
	logger.info("All tasks completed!")

	return() 	#doesn't return until the while loop exits

#-------------------------------------------------------------------------------------------#
#------------------------------------- Argparser -------------------------------------------#
#-------------------------------------------------------------------------------------------#

def restricted_float(f, f_min, f_max):
    f = float(f)
    if f < f_min or f > f_max:
        raise argparse.ArgumentTypeError("{0} not in range [0.0, 1.0]".format(f))
    return f

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



def arguments_overview(parser, args):
	""" Return string of arguments and options to print to stdout/log"""

	content = ""
	content += "# ----- Input parameters -----\n"
	for group in parser._action_groups:
			group_actions = group._group_actions
			if len(group_actions) > 0:
				#content += "# ----- {0} -----\n".format(group.title)
				for option in group_actions:
					name = option.dest
					attr = getattr(args, name, None)
					content += "# {0}:\t{1}\n".format(name, attr)
				#content += "\n"
	return(content)


#-------------------------------------------------------------------------------------------#
#---------------------------------------- Misc ---------------------------------------------#
#-------------------------------------------------------------------------------------------#

def check_required(args, required):

	for arg in required:
		if getattr(args, arg) == None:
			sys.exit("ERROR: Missing argument --{0}".format(arg))


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
	


def merge_dicts(dct, merge_dct):
	""" Merge recursive keys and values of merge_dct into dct. Values are added numerically 
		No return - dct is changed in place """

	for k, v in merge_dct.items():

		#If k is a dict, go one level down
		if (k in dct and isinstance(dct[k], dict)): # and isinstance(merge_dct[k], collections.Mapping)):
			merge_dicts(dct[k], merge_dct[k])
		else:
			if not k in dct:
				dct[k] = merge_dct[k]
			else:
				dct[k] += merge_dct[k]


def filafy(astring): 	#Make name into accepted filename
	
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
	""" Find the longest string that is a prefix of all the strings """
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
	""" Find matches between list1 and list2 (output will be the length of list1 with one or more matches per element) """

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
