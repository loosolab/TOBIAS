#!/usr/bin/env python

"""
TOBIAS top-level parser

@author: Mette Bentsen
@contact: mette.bentsen (at) mpi-bn.mpg.de
@license: MIT
"""

import sys
import argparse
from argparse import SUPPRESS
import textwrap
import importlib

#Import extra dependencies
try:
	import MOODS
except:
	sys.exit("ERROR: Package MOODS is not installed and is needed by TOBIAS. You can install it via conda using:\n"
			"$ conda install moods -c bioconda\n\n"
			"Or directly from source:\n"
			"$ wget https://github.com/jhkorhonen/MOODS/releases/download/v1.9.3/MOODS-python-1.9.3.tar.gz\n"
			"$ tar xzvf MOODS-python-1.9.3.tar.gz\n"
			"$ cd  MOODS-python-1.9.3\n"
			"$ python setup.py install"
			)

#Import parsers from tobias
from tobias.parsers import *
from tobias import __version__ as TOBIAS_VERSION

#Ignore gimmemotifs plot warning
import warnings
import matplotlib
import matplotlib.cbook
warnings.filterwarnings("ignore", category=matplotlib.cbook.mplDeprecation)


def main():

	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

	all_parser_info = {"Tools for footprinting analysis":
							{
							"ATACorrect":{"help":"Correct reads with regards to Tn5 sequence bias", "add_arguments": add_atacorrect_arguments, "function": "tobias.tools.atacorrect.run_atacorrect"},
							"ScoreBigwig":{"help":"Calculate scores such as footprints from cutsites", "add_arguments": add_scorebigwig_arguments, "function": "tobias.tools.score_bigwig.run_scorebigwig", "replaces":"FootprintScores"},
							"BINDetect":{"help":"Detect TF binding from footprints and motifs", "add_arguments": add_bindetect_arguments, "function": "tobias.tools.bindetect.run_bindetect"},
							},

						"Tools for working with motifs/TFBS":
							{
							"TFBScan": {"help":"Identify positions of TFBS given sequence and motifs", "add_arguments": add_tfbscan_arguments, "function": "tobias.tools.tfbscan.run_tfbscan"},
							"FormatMotifs": {"help": "Utility to deal with motif files", "add_arguments": add_formatmotifs_arguments, "function": "tobias.tools.format_motifs.run_formatmotifs"},
							"ClusterMotifs": {"help": "Cluster motifs by similarity", "add_arguments": add_motifclust_arguments, "function": "tobias.tools.motif_clust.run_motifclust", "space":"\t"},
							"ScoreBed": {"help":"Score .bed-file with signal from .bigwig-file(s)", "add_arguments": add_scorebed_arguments, "function": "tobias.tools.score_bed.run_scorebed"},
							},

						"Visualization tools":
							{
							"PlotAggregate": {"help": "Aggregate of .bigwig-signal across TF binding sites", "add_arguments": add_aggregate_arguments, "function": "tobias.tools.plot_aggregate.run_aggregate", "space":"\t"},
							"PlotHeatmap": {"help": "Heatmap of .bigwig-signal across TF binding sites", "add_arguments": add_heatmap_arguments, "function": "tobias.tools.plot_heatmap.run_heatmap"},
							"PlotChanges": {"help": "Plot changes in TF binding across multiple conditions (from BINDetect output)", "add_arguments": add_plotchanges_arguments, "function": "tobias.tools.plot_changes.run_plotchanges"},
							"PlotTracks": {"help": "Plot genomic tracks using the svist4get package", "add_arguments": add_tracks_arguments, "function": "tobias.tools.plot_tracks.run_tracks"}
							},

						"Miscellaneous tools":
							{
							"DownloadData": {"help": "Download test data for the TOBIAS tools", "add_arguments": add_downloaddata_arguments, "function": "tobias.tools.download_data.run_downloaddata"},
							"MergePDF": {"help": "Merge pdf files to one", "add_arguments": add_mergepdf_arguments, "function":"tobias.tools.merge_pdfs.run_mergepdf"},
							"MaxPos": {"help": "Get .bed-positions of highest bigwig signal within .bed-regions", "add_arguments": add_maxpos_arguments, "function": "tobias.tools.maxpos.run_maxpos"},
							"SubsampleBam": {"help": "Subsample a .bam-file using samtools", "add_arguments": add_subsample_arguments, "function": "tobias.tools.subsample_bam.run_subsampling"},
							"CreateNetwork": {"help": "Create TF-gene network from annotated TFBS", "add_arguments": add_network_arguments, "function": "tobias.tools.create_network.run_network", "space":"\t"},
							"Log2Table": {"help": "Convert logs from PlotAggregate to tab-delimitered tables of footprint stats", "add_arguments": add_log2table_arguments, "function": "tobias.tools.log2table.run_log2table"},
							"FilterFragments": {"help": "Tool for filtering fragments from a .bam-file based on the overlap of reads with .bed-regions", "add_arguments": add_filterfragments_arguments, "function": "tobias.tools.filter_fragments.run_filterfragments", "space":"\t"}
							}
						}

	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#					

	parser = argparse.ArgumentParser("TOBIAS", usage=SUPPRESS)
	parser._action_groups.pop()
	parser.formatter_class = lambda prog: argparse.RawDescriptionHelpFormatter(prog, max_help_position=25, width=90)
	parser.description = textwrap.dedent('''
										 ______________________________________________________________________________
										|                                                                              |
										|                                ~ T O B I A S ~                               |
										|                   Transcription factor Occupancy prediction                  |
										|                       By Investigation of ATAC-seq Signal                    |
										|______________________________________________________________________________|
									
										Usage: TOBIAS <TOOLNAME> [arguments]

										''')

	subparsers = parser.add_subparsers(title=None, metavar="")
	
	#Add all tools to parser
	all_tool_parsers = {}				
	for group in all_parser_info:
		parser.description += group + ":\n"

		info = all_parser_info[group]
		for tool in info:
			parser.description += "   {0}{1}{2}\n".format(tool, info[tool].get("space", "\t\t"), info[tool]["help"])
			subparser = subparsers.add_parser(tool, usage=SUPPRESS)
			subparser = info[tool]["add_arguments"](subparser)
			subparser.set_defaults(module=info[tool]["function"])
			all_tool_parsers[tool.lower()] = subparser

			#Add version to subparser
			subparser.add_argument("--version", action='version', version=TOBIAS_VERSION)
			subparser = add_underscore_options(subparser)

			#Add parser for old tool names
			if "replaces" in info[tool]:
				replace_tool = info[tool]["replaces"]
				subparser = subparsers.add_parser(replace_tool, usage=SUPPRESS)
				subparser = info[tool]["add_arguments"](subparser)
				subparser.set_defaults(module=info[tool]["function"])
				all_tool_parsers[replace_tool.lower()] = subparser
			
		parser.description += "\n"

	parser.description += "For help on each tool, please run: TOBIAS <TOOLNAME> --help\n"
	
	#Add version number to upper TOBIAS parser 
	parser.description += "For version number: TOBIAS --version"
	parser.add_argument("--version", action='version', version=TOBIAS_VERSION)

	#If no args, print help for top-level TOBIAS
	if len(sys.argv[1:]) == 0:
		parser.print_help()
		sys.exit()

	#if args are pointing to specific tool, and no other arguments are given, print help for this
	if sys.argv[1].lower() in all_tool_parsers and len(sys.argv) == 2:
		chosen_tool = sys.argv[1]
		if chosen_tool != "DownloadData":	#Downloaddata can be run without options
			all_tool_parsers[chosen_tool.lower()].print_help()
			sys.exit()
	
	args = parser.parse_args()

	#Depending on subparser chosen, load main script entry and run
	function_str = args.module
	mod_name, func_name = function_str.rsplit(".", 1)
	module = importlib.import_module(mod_name)	#load specific module
	func = getattr(module, func_name)
	args.func = func

	#Run specified function with arguments
	args.func(args)		

	
if __name__ == "__main__":
    main()