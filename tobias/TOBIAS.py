#!/usr/bin/env python

"""
TOBIAS top-level parser

@author: Mette Bentsen
@contact: mette.bentsen (at) mpi-bn.mpg.de
@license: MIT
"""

import sys

#Import extra dependencies
try:
	import MOODS
except:
	sys.exit("ERROR: Package MOODS is not installed and is needed by TOBIAS. You can install it via conds using:\n"
			"$ conda install moods -c bioconda\n\n"
			"Or directly from source:\n"
			"$ wget https://github.com/jhkorhonen/MOODS/releases/download/v1.9.3/MOODS-python-1.9.3.tar.gz\n"
			"$ tar xzvf MOODS-python-1.9.3.tar.gz\n"
			"$ cd  MOODS-python-1.9.3\n"
			"$ python setup.py install"
			)

#Import general 
import argparse
from argparse import SUPPRESS
import textwrap

from tobias.footprinting.atacorrect import *
from tobias.footprinting.scorebigwig import *
from tobias.footprinting.bindetect import *

from tobias.plotting.plot_aggregate import *
from tobias.plotting.plot_heatmap import *
from tobias.plotting.plot_changes import *

from tobias.motifs.tfbscan import * 
from tobias.motifs.format_motifs import * 
#from tobias.motifs.cluster_tfbs import *
from tobias.motifs.score_bed import *

from tobias.misc.subsample_bam import *
from tobias.misc.merge_pdfs import *
from tobias.misc.maxpos import *
#from tobias.misc.create_network import *
from tobias.misc.log2table import *

from tobias import __version__ as TOBIAS_VERSION

def main():

	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

	all_parser_info = {"Tools for footprinting analysis":
							{
							"ATACorrect":{"help":"Correct reads with regards to Tn5 sequence bias", "add_arguments": add_atacorrect_arguments, "function":run_atacorrect},
							"ScoreBigwig":{"help":"Calculate scores such as footprints from cutsites", "add_arguments": add_scorebigwig_arguments, "function":run_scorebigwig, "replaces":"FootprintScores"},
							"BINDetect":{"help":"Detect TF binding from footprints and motifs", "add_arguments": add_bindetect_arguments, "function":run_bindetect},
							},

						"Tools for working with motifs/TFBS":
							{
							"TFBScan": {"help":"Identify positions of TFBS given sequence and motifs", "add_arguments": add_tfbscan_arguments, "function": run_tfbscan},
							"FormatMotifs": {"help": "Utility to deal with motif files", "add_arguments": add_formatmotifs_arguments, "function": run_formatmotifs},
							#"ClusterTFBS": {"help": "Cluster TFs based on overlap of sites", "add_arguments": add_clustering_arguments, "function": run_clustering},
							"ScoreBed": {"help":"Score .bed-file with signal from .bigwig-file(s)", "add_arguments": add_scorebed_arguments, "function": run_scorebed},
							},

						"Visualization tools":
							{
							"PlotAggregate": {"help": "Aggregate of .bigwig-signal across TF binding sites", "add_arguments": add_aggregate_arguments, "function": run_aggregate, "space":"\t"},
							"PlotHeatmap": {"help": "Heatmap of .bigwig-signal across TF binding sites", "add_arguments": add_heatmap_arguments, "function": run_heatmap},
							"PlotChanges": {"help": "Plot changes in TF binding across multiple conditions (from BINDetect output)", "add_arguments": add_plotchanges_arguments, "function": run_plotchanges}
							},

						"Miscellaneous tools":
							{
							"MergePDF": {"help": "Merge pdf files to one", "add_arguments":add_mergepdf_arguments, "function":run_mergepdf},
							"MaxPos": {"help": "Get .bed-positions of highest bigwig signal within .bed-regions", "add_arguments": add_maxpos_arguments, "function": run_maxpos},
							"SubsampleBam": {"help": "Subsample a .bam-file using samtools", "add_arguments": add_subsample_arguments, "function": run_subsampling},
							#"CreateNetwork": {"help": "Create TF-gene network from annotated TFBS", "add_arguments": add_network_arguments, "function": run_network, "space":"\t"},
							"Log2Table": {"help": "Convert logs from PlotAggregate to tab-delimitered tables of footprint stats", "add_arguments": add_log2table_arguments, "function": run_log2table}
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
			subparser.set_defaults(func=info[tool]["function"])
			all_tool_parsers[tool.lower()] = subparser

			#Add version to subparser
			subparser.add_argument("--version", action='version', version=TOBIAS_VERSION)

			#Add parser for old tool names
			if "replaces" in info[tool]:
				replace_tool = info[tool]["replaces"]
				subparser = subparsers.add_parser(replace_tool, usage=SUPPRESS)
				subparser = info[tool]["add_arguments"](subparser)
				subparser.set_defaults(func=info[tool]["function"])
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
		all_tool_parsers[chosen_tool.lower()].print_help()
		sys.exit()
	
	args = parser.parse_args()
	args.func(args)		#run specified function with arguments

if __name__ == "__main__":
    main()