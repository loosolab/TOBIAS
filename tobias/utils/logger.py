#!/usr/bin/env python

"""
Class for dealing with logger-function used across all TOBIAS tools

@author: Mette Bentsen
@contact: mette.bentsen (at) mpi-bn.mpg.de
@license: MIT

"""

import sys
import os
from datetime import datetime
import logging
import logging.handlers
import multiprocessing as mp
import time

def add_logger_args(args):
	""" Function for adding TOBIAS-wide verbosity to command-line parsers """
	args.add_argument('--verbosity', metavar="<int>", help="Level of output logging (0: silent, 1: errors/warnings, 2: info, 3: stats, 4: debug, 5: spam) (default: 3)", choices=[0,1,2,3,4,5], default=3, type=int)
	return(args)

class TobiasLogger(logging.Logger):
	""" TobiasLogger is an instance of a logging.Logger with special functions for formatting and creating automatic logging """

	logger_levels = {
						0: 0,
						1: logging.WARNING,							#also includes errors
						2: logging.INFO, 							#info 
						3: int((logging.INFO+logging.DEBUG)/2),		#statistics
						4: logging.DEBUG,							#debugging info
						5: logging.DEBUG - 5						#spam-level debug
					}

	def __init__(self, tool_name="TOBIAS", level=3, queue=None):

		self.tool_name = tool_name		#name of tool within TOBIAS
		logging.Logger.__init__(self, self.tool_name)

		if level == 0:
			self.disabled = True

		####### Setup custom levels #######
		#Create custom level for comments (Same level as errors/warnings)
		comment_level = TobiasLogger.logger_levels[1] + 1
		logging.addLevelName(comment_level, "comment") #log_levels[lvl])
		setattr(self, 'comment', lambda *args: self.log(comment_level, *args))

		#Create custom level for stats (between info and debug)
		stats_level = TobiasLogger.logger_levels[3]
		logging.addLevelName(stats_level, "STATS") #log_levels[lvl])
		setattr(self, 'stats', lambda *args: self.log(stats_level, *args))
		
		#Create custom level for spamming debug messages
		spam_level = TobiasLogger.logger_levels[5]
		logging.addLevelName(spam_level, "SPAM") #log_levels[lvl])
		setattr(self, 'spam', lambda *args: self.log(spam_level, *args))

		#Set level
		self.level = TobiasLogger.logger_levels[level]
		
		######### Setup formatter ########

		#Setup formatting
		self.formatter = TOBIASFormatter()
		self.setLevel(self.level)
		
		########## Setup streaming #########
		##Log file stream
		#if log_f != None:
		#	log = logging.FileHandler(log_f, "w")
		#	log.setLevel(self.level)
		#	log.setFormatter(self.formatter)
		#	self.addHandler(log)

		if queue == None:
			#Stdout stream
			con = logging.StreamHandler(sys.stdout)		#console output
			con.setLevel(self.level)
			con.setFormatter(self.formatter)
			self.addHandler(con)
		else:
			h = logging.handlers.QueueHandler(queue)  	# Just the one handler needed
			self.handlers = []
			self.addHandler(h)

		#Lastly, initialize time
		self.begin_time = datetime.now()
		self.end_time = None
		self.total_time = None


	def begin(self):
		""" Begin logging by writing comments about the current run """
		from tobias import __version__ as TOBIAS_VERSION

		#Print info on run
		self.comment("# TOBIAS {0} {1} (run started {2})".format(TOBIAS_VERSION, self.tool_name, self.begin_time))
		self.comment("# Working directory: {0}".format(os.getcwd()))
		self.comment("# Command line call: {0}\n".format(" ".join(sys.argv)))

	def stop(self):
		""" Stop without printing status """
		
		self.end_time = datetime.now()
		self.total_time = self.end_time - self.begin_time
		
	def end(self):
		""" End logging - write out how long it took """

		self.end_time = datetime.now()
		self.total_time = self.end_time - self.begin_time

		self.comment("")	#creates empty line; only for pretty output
		self.info("Finished {0} run (total time elapsed: {1})".format(self.tool_name, self.total_time))


	def start_logger_queue(self):
		""" start process for listening and handling through the main logger queue """

		self.debug("Starting logger queue for multiprocessing")
		self.queue = mp.Manager().Queue()
		self.listener = mp.Process(target=self.main_logger_process)
		self.listener.start()
		

	def stop_logger_queue(self):
		""" Stop process for listening """

		self.debug("Waiting for listener to finish")
		self.queue.put(None)
		while self.listener.exitcode != 0:
			self.debug("Listener exitcode is: {0}. Waiting for exitcode = 0.".format(self.listener.exitcode))
			time.sleep(0.1)

		self.debug("Joining listener")
		self.listener.join()


	def main_logger_process(self):

		self.debug("Started main logger process")
		while True:
			try:
				record = self.queue.get()
				if record is None:
					break
				self.handle(record) 	#this logger is coming from the main process

			except Exception:
				import sys, traceback
				print('Problem in main logger process:', file=sys.stderr)
				traceback.print_exc(file=sys.stderr)
				break

		return(1)


	def arguments_overview(self, parser, args):
		""" Creates string of arguments and options to print using logger """

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
		self.comment(content + "\n")


	def output_files(self, outfiles):
		""" Print out list of output files"""

		self.comment("# ----- Output files -----")
		for outf in outfiles:
			if outf != None:
				self.comment("# {0}".format(outf))
		self.comment("\n")



class TOBIASFormatter(logging.Formatter):
	""" Formatter class used in TobiasLogger """
	default_fmt = logging.Formatter("%(asctime)s (%(process)d) [%(levelname)s]\t%(message)s", "%Y-%m-%d %H:%M:%S")
	comment_fmt = logging.Formatter("%(message)s")

	def format(self, record):
		format_orig = self._fmt

		#Comments
		if record.levelname == "comment":
			return self.comment_fmt.format(record)
		elif record.levelno != 0:
			return self.default_fmt.format(record)	
		else:
			return
