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

from tobias.footprinting.ATACorrect import *
from tobias.footprinting.footprint_scores import *
from tobias.footprinting.BINDetect import *

from tobias.plotting.plot_aggregate import *
from tobias.plotting.plot_heatmap import *
from tobias.plotting.plot_bindetect import *
from tobias.plotting.plot_changes import *

from tobias.motifs.tfbscan import * 
from tobias.motifs.format_motifs import * 

from tobias.utils.subsample_bam import *
from tobias.utils.merge_pdfs import *
from tobias.utils.score_bed import *

def main():
	parser = argparse.ArgumentParser("TOBIAS", usage=SUPPRESS)
	parser._action_groups.pop()
	parser.formatter_class = lambda prog: argparse.RawDescriptionHelpFormatter(prog, max_help_position=25, width=90)
	parser.description = textwrap.dedent('''
										 ______________________________________________________________________________
										|                                                                              |
										|                                ~ T O B I A S ~                               |
										|                    Transcription factor Occupancy prediction                 |
										|                       By Investigation of ATAC-seq Signal                    |
										|______________________________________________________________________________|
									
										Usage: TOBIAS <TOOLNAME> [arguments]

										''')

	subparsers = parser.add_subparsers(title=None, metavar="")
	all_tool_parsers = {}

	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Tools for footprinting ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

	parser.description += "Tools for footprinting analysis:\n"

	name, hlp = "ATACorrect", "Correct reads with regards to Tn5 sequence bias"
	parser.description += "   {0}\t\t{1}\n".format(name, hlp)
	atacorrect_parser = subparsers.add_parser(name, usage=SUPPRESS)
	atacorrect_parser = add_atacorrect_arguments(atacorrect_parser) 	#add atacorrect arguments to the atacorrect subparser
	atacorrect_parser.set_defaults(func=run_atacorrect)
	all_tool_parsers[name.lower()] = atacorrect_parser

	name, hlp = "FootprintScores", "Calculate footprint scores from cutsites"
	parser.description += "   {0}\t{1}\n".format(name, hlp)
	footprint_parser = subparsers.add_parser(name, usage=SUPPRESS)
	footprint_parser = add_footprint_arguments(footprint_parser)
	footprint_parser.set_defaults(func=run_footprinting)
	all_tool_parsers[name.lower()] = footprint_parser

	name, hlp = "BINDetect", "Detects TF binding from footprints"
	parser.description += "   {0}\t\t{1}\n".format(name, hlp)
	detect_parser = subparsers.add_parser(name, usage=SUPPRESS)
	detect_parser = add_bindetect_arguments(detect_parser)
	detect_parser.set_defaults(func=run_bindetect)
	all_tool_parsers[name.lower()] = detect_parser


	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~ Tools for working with TFBS ~~~~~~~~~~~~~~~~~~~~~~~~~#
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

	parser.description += "\nTools for working with motifs/TFBS:\n"
	
	name, hlp = "TFBScan", "Identify positions of TFBS given sequence and motifs"
	parser.description += "   {0}\t\t{1}\n".format(name, hlp)
	tfbscan_parser = subparsers.add_parser(name, usage=SUPPRESS)
	tfbscan_parser = add_tfbscan_arguments(tfbscan_parser)
	tfbscan_parser.set_defaults(func=run_tfbscan)
	all_tool_parsers[name.lower()] = tfbscan_parser
	
	name, hlp = "FormatMotifs", "Utility to deal with motif files"
	parser.description += "   {0}\t\t{1}\n".format(name, hlp)
	formatmotifs_parser = subparsers.add_parser(name, usage=SUPPRESS) 
	formatmotifs_parser = add_formatmotifs_arguments(formatmotifs_parser)
	formatmotifs_parser.set_defaults(func=run_formatmotifs)
	all_tool_parsers[name.lower()] = formatmotifs_parser

	"""
	name, hlp = "ClusterTF", "Cluster TFs based on overlap of sites"
	parser.description += "   {0}\t\t{1}\n".format(name, hlp)
	clustering_parser = subparsers.add_parser(name, usage=SUPPRESS)
	clustering_parser = add_clustering_arguments(clustering_parser)
	clustering_parser.set_defaults(func=run_clustering)
	all_tool_parsers[name] = clustering_parser
	"""
	
	name, hlp = "ScoreBed", "Score .bed-file with signal from .bigwig-file(s)"
	parser.description += "   {0}\t\t{1}\n".format(name, hlp)
	scorebed_parser = subparsers.add_parser(name, usage=SUPPRESS) 
	scorebed_parser = add_scorebed_arguments(scorebed_parser)
	scorebed_parser.set_defaults(func=run_scorebed)
	all_tool_parsers[name.lower()] = scorebed_parser



	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Tools for plotting ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

	parser.description += "\nVisualization tools:\n"

	name, hlp = "PlotAggregate", "Aggregate of .bigwig-signal across TF binding sites"
	parser.description += "   {0}\t{1}\n".format(name, hlp)
	aggregate_parser = subparsers.add_parser(name, usage=SUPPRESS)
	aggregate_parser = add_aggregate_arguments(aggregate_parser)
	aggregate_parser.set_defaults(func=run_aggregate)
	all_tool_parsers[name.lower()] = aggregate_parser

	name, hlp = "PlotHeatmap", "Heatmap of .bigwig-signal across TF binding sites"
	parser.description += "   {0}\t\t{1}\n".format(name, hlp)
	heatmap_parser = subparsers.add_parser(name, usage=SUPPRESS)
	heatmap_parser = add_heatmap_arguments(heatmap_parser)
	heatmap_parser.set_defaults(func=run_heatmap)
	all_tool_parsers[name.lower()] = heatmap_parser
	
	name, hlp = "PlotBINDetect", "Plotting function from BINDetect (to re-plot output)"
	parser.description += "   {0}\t{1}\n".format(name, hlp)
	diffplot_parser = subparsers.add_parser(name, usage=SUPPRESS)
	diffplot_parser = add_diffplot_arguments(diffplot_parser)
	diffplot_parser.set_defaults(func=run_diffplot)
	all_tool_parsers[name.lower()] = diffplot_parser

	name, hlp = "PlotChanges", "Plot changes in TF binding across multiple conditions (from BINDetect output)"
	parser.description += "   {0}\t\t{1}\n".format(name, hlp)
	changeplot_parser = subparsers.add_parser(name, usage=SUPPRESS)
	changeplot_parser = add_plotchanges_arguments(changeplot_parser)
	changeplot_parser.set_defaults(func=run_plotchanges)
	all_tool_parsers[name.lower()] = changeplot_parser



	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Misc tools ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

	parser.description += "\nMiscellaneous tools:\n"

	name, hlp = "MergePDF", "Merge pdf files to one"
	parser.description += "   {0}\t\t{1}\n".format(name, hlp)
	mergepdf_parser = subparsers.add_parser(name, usage=SUPPRESS)
	mergepdf_parser = add_mergepdf_arguments(mergepdf_parser)
	mergepdf_parser.set_defaults(func=run_mergepdf)
	all_tool_parsers[name.lower()] = mergepdf_parser

	name, hlp = "SubsampleBam", "Subsample a .bam-file using samtools"
	parser.description += "   {0}\t\t{1}\n".format(name, hlp)
	subsample_parser = subparsers.add_parser(name, usage=SUPPRESS)
	subsample_parser = add_subsample_arguments(subsample_parser)
	subsample_parser.set_defaults(func=run_subsampling)
	all_tool_parsers[name.lower()] = subsample_parser

	parser.description += "\n"

	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

	parser.description += "For help on each tool, please run: TOBIAS <TOOLNAME> --help"

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