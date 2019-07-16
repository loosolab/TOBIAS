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
	check_files(motif_files + [args.filter]) 	#Check if files exist

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
		logger.debug("- {0}".format(f))
		motif_list.extend(MotifList().from_file(f))

	logger.info("Read {} motifs\n".format(len(motif_list)))

	#Sort out duplicate motifs
	all_motif_ids = [motif.id for motif in motif_list]
	unique_motif_ids = set(all_motif_ids)
	if len(all_motif_ids) != len(unique_motif_ids):
		logger.info("Found duplicate motif ids in file - choosing first motif with unique id.")
		motif_list = MotifList([motif_list[all_motif_ids.index(motifid)] for motifid in unique_motif_ids])
		logger.info("Reduced to {0} unique motif ids".format(len(motif_list)))

	### Filter motif list ###
	if args.filter != None:

		#Read filter
		logger.info("Reading entries in {0}".format(args.filter))
		entries = open(args.filter, "r").read().split()
		logger.info("Read {0} unique filter values".format(len(set(entries))))

		#Match to input motifs #print(entries)
		logger.info("Matching motifs to filter")
		used_filters = []
		filtered_list = MotifList()
		for input_motif in motif_list:
			found_in_filter = 0
			i = -1
			while found_in_filter == 0 and i < len(entries) - 1:
				i += 1 
				if entries[i].lower() in input_motif.name.lower() or entries[i].lower() in input_motif.id.lower():
					filtered_list.append(input_motif)
					logger.debug("Selected motif {0} ({1}) due to filter value {2}".format(input_motif.name, input_motif.id, entries[i]))
					found_in_filter = 1
					used_filters.append(entries[i])

			if found_in_filter == 0:
				logger.debug("Could not find any match to motif {0} ({1}) in filter".format(input_motif.name, input_motif.id))

		logger.info("Filtered number of motifs from {0} to {1}".format(len(motif_list), len(filtered_list)))
		motif_list = filtered_list

		logger.debug("Filters not used: {0}".format(list(set(entries) - set(used_filters))))

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
