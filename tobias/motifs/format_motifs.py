#!/usr/bin/env python

"""
Utility to format (convert/join/split) motifs between pfm/jaspar/meme formats

@author: Mette Bentsen
@contact: mette.bentsen (at) mpi-bn.mpg.de
@license: MIT
"""

import argparse
import sys
import os
import re
import textwrap

#Internal functions/classes
from tobias.utils.motifs import *
from tobias.utils.utilities import *
from tobias.utils.logger import *

#--------------------------------------------------------------------------------------------------------#
def add_formatmotifs_arguments(parser):


	parser.formatter_class = lambda prog: argparse.RawDescriptionHelpFormatter(prog, max_help_position=40, width=90)
	description = ""
	parser.description = format_help_description("FormatMotifs", description) 
	
	parser._action_groups.pop()	#pop -h

	#Required arguments
	required = parser.add_argument_group('Required arguments')
	required.add_argument('--input', metavar="", nargs="*", help="One or more input motif files")			
	required.add_argument('--format', metavar="", help="Desired motif output format (pfm, jaspar, meme) (default: \"jaspar\")", choices=["pfm", "jaspar", "meme"], default="jaspar")
	required.add_argument('--task', metavar="", help="Which task to perform on motif files (join/split) (default: join)", choices=["join", "split"], default="join")
	required.add_argument('--filter', metavar="", help="File containing list of motif names/ids to filter on. Only motifs fitting entries in filter will be output.")
	required.add_argument('--output', metavar="", help="If task == join, output is the joined output file; if task == split, output is a directory")
	
	additional = parser.add_argument_group('Additional arguments')
	additional = add_logger_args(additional)

	return(parser)


#--------------------------------------------------------------------------------------------------------#
def run_formatmotifs(args):

	check_required(args, ["input", "output"])	#Check input arguments
	motif_files = expand_dirs(args.input)		#Expand any dirs in input
	check_files(motif_files) 					#Check if files exist

	# Create logger and write argument overview
	logger = TobiasLogger("FormatMotifs", args.verbosity)
	logger.begin()

	parser = add_formatmotifs_arguments(argparse.ArgumentParser())
	logger.arguments_overview(parser, args)
	logger.output_files([args.output])

	####### Getting ready #######
	if args.task == "split":
		logger.info("Making directory {0} if not existing".format(args.output))
		make_directory(args.output)		#Check and create output directory

	### Read motifs from files ###
	logger.info("Reading input files...")
	motif_list = MotifList()
	converted_content = ""
	for f in motif_files:
		logger.info("- {0}".format(f))
		motif_list.extend(MotifList().from_file(f))

	logger.info("Read {} motifs\n".format(len(motif_list)))

	### Filter motif list ###
	if args.filter != None:

		#Read filter
		pass


	#### Write out results ####
	if args.task == "split":
		logger.info("Writing individual files to directory {0}".format(args.output))

		for motif in motif_list:
			motif_string = MotifList([motif]).as_string(args.format)

			#Open file and write
			out_path = os.path.join(args.output, motif.id + "." + args.format)
			logger.info("- {0}".format(out_path))
			f_out = open(out_path, "w")
			f_out.write(motif_string)
			f_out.close()
	
	elif args.task == "join":
		logger.info("Writing converted motifs to file {0}".format(args.output))

		f_out = open(args.output, "w")
		motif_string = motif_list.as_string(args.format)
		f_out.write(motif_string)
		f_out.close()

	logger.end()

#--------------------------------------------------------------------------------------------------------#

if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser = add_formatmotifs_arguments(parser)
	args = parser.parse_args()

	if len(sys.argv[1:]) == 0:
		parser.print_help()
		sys.exit()
		
	run_formatmotifs(args)
