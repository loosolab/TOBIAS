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
	required.add_argument('--output', metavar="", help="If task == join, output is the joined output file; if task == split, output is a directory")
	
	additional = parser.add_argument_group('Additional arguments')
	additional = add_logger_args(additional)

	return(parser)


#--------------------------------------------------------------------------------------------------------#
def run_formatmotifs(args):

	check_required(args, ["input", "output"])	#Check input arguments
	check_files(args.input) 					#Check if files exist

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
	else:
		logger.info("Opening file {0} for writing".format(args.output))
		out_f = open(args.output, "w")	#open file

	### Estimate format of input files
	logger.info("Reading input files...")
	converted_content = ""
	for f in args.input:
		content = open(f).read()
		motif_format = get_motif_format(content)
		logger.info("- {0}: {1}".format(f, motif_format))

		#Convert to output format
		converted_content += convert_motif(content, args.format)

	logger.comment("")

	#### Write out results ####
	#todo: if meme is input, estimate meme header from input file; else set it to standard
	meme_header = "MEME version 4\n\n"
	meme_header += "ALPHABET=ACGT\n\n"
	meme_header += "strands: + -\n\n"
	meme_header += "Background letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n\n"

	if args.task == "split":
		logger.info("Writing individual files to directory {0}".format(args.output))

		#Split on ">"
		if args.format == "jaspar" or args.format == "pfm":
			splitted_content = [ ">" + content for content in converted_content.split(">") if content != ""]
		elif args.format == "meme":
			splitted_content = [content for content in converted_content.split("MOTIF") if content != ""]
			splitted_content = [content if content.startswith("MOTIF") else "MOTIF" + content for content in splitted_content]
			
		#Write out
		for motif in splitted_content:

			#Get name from headerline
			elements = motif.replace("MOTIF", "").replace(">", "").split()
			motifid, name = elements[0], elements[1]


			#Open file and write
			out_path = os.path.join(args.output, motifid + "." + args.format)
			logger.info("- {0}".format(out_path))
			f_out = open(out_path, "w")

			if args.format == "meme":
				f_out.write(meme_header + motif) #make sure every motif has a header
			else:
				f_out.write(motif)

			f_out.close()
	
	elif args.task == "join":
		logger.info("Writing converted motifs to file {0}".format(args.output))

		if args.format == "meme":
			converted_content = meme_header + converted_content

		out_f.write(converted_content + "\n")

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
