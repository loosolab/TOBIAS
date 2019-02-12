#!/usr/bin/env python

"""
Small utility for mering pdfs

@author: Mette Bentsen
@contact: mette.bentsen (at) mpi-bn.mpg.de
@license: MIT
"""

from PyPDF2 import PdfFileMerger, PdfFileReader
import argparse
import sys
import os

#Internal functions
from tobias.utils.utilities import *


#--------------------------------------------------------------------------------------------------------#
def add_mergepdf_arguments(parser):

	parser.formatter_class = lambda prog: argparse.RawDescriptionHelpFormatter(prog, max_help_position=40, width=90)
	description = "Merge single PDF-files to one file"
	parser.description = format_help_description("MergePDF", description)
	parser._action_groups.pop()	#pop -h

	reqargs = parser.add_argument_group('Required arguments')
	reqargs.add_argument('--input', metavar="", nargs="*", help="PDF files to join")
	reqargs.add_argument('--output', metavar="", help="Path to output file (default: ./merged.pdf)", default="merged.pdf")

	return(parser)


def run_mergepdf(args):

	check_required(args, ["input", "output"])
	print("Number of input files: {0}".format(len(args.input)))

	#Preliminary checks
	print("Checking read/write status")
	check_files(args.input, action="r")
	check_files([args.output], action="w")

	#Join pdfs
	print("Starting to merge PDFs")
	merger = PdfFileMerger(strict=False)
	for pdf in args.input:
		if os.stat(pdf).st_size != 0:	#only join files containing plots
			merger.append(PdfFileReader(pdf))
	
	print("Writing merged file: {0}".format(args.output))
	merger.write(args.output)

	print("PDFs merged successfully!")


#--------------------------------------------------------------------------------------------------------#
if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser = add_mergepdf_arguments(parser)
	args = parser.parse_args()

	if len(sys.argv[1:]) == 0:
		parser.print_help()
		sys.exit()

	run_mergepdf(args)
