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
from tobias.parsers import add_formatmotifs_arguments
from tobias.utils.motifs import OneMotif, MotifList
from tobias.utils.utilities import * #check_required, check_files, expand_dirs
from tobias.utils.logger import TobiasLogger

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
