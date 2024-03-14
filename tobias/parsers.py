#!/usr/bin/env python

"""
TOBIAS top-level parser

@author: Mette Bentsen
@contact: mette.bentsen (at) mpi-bn.mpg.de
@license: MIT
"""

import argparse
from tobias.utils.utilities import format_help_description, restricted_float, add_underscore_options
from tobias.utils.logger import add_logger_args

#--------------------------------------------------------------------------------------------------------#
def add_atacorrect_arguments(parser):

	parser.formatter_class = lambda prog: argparse.RawDescriptionHelpFormatter(prog, max_help_position=35, width=90)

	description = "ATACorrect corrects the cutsite-signal from ATAC-seq with regard to the underlying sequence preference of Tn5 transposase.\n\n"
	description += "Usage:\nTOBIAS ATACorrect --bam <reads.bam> --genome <genome.fa> --peaks <peaks.bed>\n\n"
	description += "Output files:\n"
	description += "\n".join(["- <outdir>/<prefix>_{0}.bw".format(track) for track in ["uncorrected", "bias", "expected", "corrected"]]) + "\n"
	description += "- <outdir>/<prefix>_atacorrect.pdf"
	parser.description = format_help_description("ATACorrect", description)

	parser._action_groups.pop()	#pop -h

	#Required arguments
	reqargs = parser.add_argument_group('Required arguments')
	reqargs.add_argument('-b', '--bam', metavar="<bam>", help="A .bam-file containing reads to be corrected")
	reqargs.add_argument('-g', '--genome', metavar="<fasta>", help="A .fasta-file containing whole genomic sequence")
	reqargs.add_argument('-p', '--peaks', metavar="<bed>", help="A .bed-file containing ATAC peak regions")

	#Optional arguments
	optargs = parser.add_argument_group('Optional arguments')
	optargs.add_argument('--regions-in', metavar="<bed>", help="Input regions for estimating bias (default: regions not in peaks.bed)")
	optargs.add_argument('--regions-out', metavar="<bed>", help="Output regions (default: peaks.bed)")
	optargs.add_argument('--blacklist', metavar="<bed>", help="Blacklisted regions in .bed-format (default: None)") #file containing blacklisted regions to be excluded from analysis")
	optargs.add_argument('--extend', metavar="<int>", type=int, help="Extend output regions with basepairs upstream/downstream (default: 100)", default=100)
	optargs.add_argument('--split-strands', help="Write out tracks per strand", action="store_true")
	optargs.add_argument('--norm-off', help="Switches off normalization based on number of reads", action='store_true')
	optargs.add_argument('--track-off', metavar="<track>", help="Switch off writing of individual .bigwig-tracks (uncorrected/bias/expected/corrected)", nargs="*", choices=["uncorrected", "bias", "expected", "corrected"], default=[])
	optargs.add_argument('--drop-chroms', metavar="<chrom>", help="Drop any chromosomes in the list from the correction. The default is to drop the mitochrondrial chromosome. Default: ['chrM', 'chrMT', 'M', 'MT', 'Mito']", nargs="*", default=['chrM', 'chrMT', 'M', "MT", "Mito"])

	optargs = parser.add_argument_group('Advanced ATACorrect arguments (no need to touch)')
	optargs.add_argument('--k_flank', metavar="<int>", help="Flank +/- of cutsite to estimate bias from (default: 12)", type=int, default=12)
	optargs.add_argument('--read_shift', metavar="<int>", help="Read shift for forward and reverse reads (default: 4 -5)", nargs=2, type=int, default=[4,-5])
	optargs.add_argument('--bg_shift', metavar="<int>", type=int, help="Read shift for estimation of background frequencies (default: 100)", default=100)
	optargs.add_argument('--window', metavar="<int>", help="Window size for calculating expected signal (default: 100)", type=int, default=100)
	optargs.add_argument('--score_mat', metavar="<mat>", help="Type of matrix to use for bias estimation (PWM/DWM) (default: DWM)", choices=["PWM", "DWM"], default="DWM")
	optargs.add_argument('--bias-pkl', metavar="<obj>", help="Path to a pre-calculated AtacBias.pkl-object, as output from a previous ATACorrect run (default: None). Can be used to bypass the internal bias estimation.", default=None)

	runargs = parser.add_argument_group('Run arguments')
	runargs.add_argument('--prefix', metavar="<prefix>", help="Prefix for output files (default: same as .bam file)")
	runargs.add_argument('--outdir', metavar="<directory>", help="Output directory for files (default: current working directory)", default="")
	runargs.add_argument('--cores', metavar="<int>", type=int, help="Number of cores to use for computation (default: 1)", default=1)
	runargs.add_argument('--split', metavar="<int>", type=int, help="Split of multiprocessing jobs (default: 100)", default=100)
	
	runargs = add_logger_args(runargs)

	return(parser)

#--------------------------------------------------------------------------------------------------------#
def add_scorebigwig_arguments(parser):

	parser.formatter_class = lambda prog: argparse.RawDescriptionHelpFormatter(prog, max_help_position=40, width=90)
	description = "ScoreBigwig calculates scores (such as footprint-scores) from bigwig files (such as ATAC-seq cutsites calculated using the ATACorrect tool). " 
	description += "NOTE: ScoreBigwig replaces the previous FootprintScores command.\n\n"
	description += "Usage: ScoreBigwig --signal <cutsites.bw> --regions <regions.bed> --output <output.bw>\n\n"
	description += "Output:\n- <output.bw>"
	parser.description = format_help_description("ScoreBigwig", description)
	
	parser._action_groups.pop()	#pop -h

	#Required arguments
	required = parser.add_argument_group('Required arguments')
	required.add_argument('-s', '--signal', metavar="<bigwig>", help="A .bw file of ATAC-seq cutsite signal")
	required.add_argument('-o', '--output', metavar="<bigwig>", help="Full path to output bigwig")			
	required.add_argument('-r', '--regions', metavar="<bed>", help="Genomic regions to run footprinting within")

	optargs = parser.add_argument_group('Optional arguments')
	optargs.add_argument('--score', metavar="<score>", choices=["footprint", "sum", "mean", "none"], help="Type of scoring to perform on cutsites (footprint/sum/mean/none) (default: footprint)", default="footprint")
	optargs.add_argument('--absolute', action='store_true', help="Convert bigwig signal to absolute values before calculating score")
	optargs.add_argument('--extend', metavar="<int>", type=int, help="Extend input regions with bp (default: 100)", default=100)
	optargs.add_argument('--smooth', metavar="<int>", type=int, help="Smooth output signal by mean in <bp> windows (default: no smoothing)", default=1)
	optargs.add_argument('--min-limit', metavar="<float>", type=float, help="Limit input bigwig value range (default: no lower limit)") 		#default none
	optargs.add_argument('--max-limit', metavar="<float>", type=float, help="Limit input bigwig value range (default: no upper limit)") 		#default none

	footprintargs = parser.add_argument_group('Parameters for score == footprint')
	footprintargs.add_argument('--fp-min', metavar="<int>", type=int, help="Minimum footprint width (default: 20)", default=20)
	footprintargs.add_argument('--fp-max', metavar="<int>", type=int, help="Maximum footprint width (default: 50)", default=50)
	footprintargs.add_argument('--flank-min', metavar="<int>", type=int, help="Minimum range of flanking regions (default: 10)", default=10)
	footprintargs.add_argument('--flank-max', metavar="<int>", type=int, help="Maximum range of flanking regions (default: 30)", default=30)
	
	sumargs = parser.add_argument_group('Parameters for score == sum')
	sumargs.add_argument('--window', metavar="<int>", type=int, help="The window for calculation of sum (default: 100)", default=100)

	runargs = parser.add_argument_group('Run arguments')
	runargs.add_argument('--cores', metavar="<int>", type=int, help="Number of cores to use for computation (default: 1)", default=1)
	runargs.add_argument('--split', metavar="<int>", type=int, help="Split of multiprocessing jobs (default: 100)", default=100)
	runargs = add_logger_args(runargs)

	return(parser)

#--------------------------------------------------------------------------------------------------------#
def add_bindetect_arguments(parser):

	parser.formatter_class = lambda prog: argparse.RawDescriptionHelpFormatter(prog, max_help_position=35, width=90)
	description = "BINDetect takes motifs, signals (footprints) and genome as input to estimate bound transcription factor binding sites and differential binding between conditions. "
	description += "The underlying method is a modified motif enrichment test to see which motifs have the largest differences in signal across input conditions. "
	description += "The output is an in-depth overview of global changes as well as the individual binding site signal-differences.\n\n"
	description += "Usage:\nTOBIAS BINDetect --signals <bigwig1> (<bigwig2> (...)) --motifs <motifs.txt> --genome <genome.fasta> --peaks <peaks.bed>\n\n"
	description += "Output files:\n- <outdir>/<prefix>_figures.pdf\n- <outdir>/<prefix>_results.{txt,xlsx}\n- <outdir>/<prefix>_distances.txt\n"
	description += "- <outdir>/<TF>/<TF>_overview.{txt,xlsx} (per motif)\n- <outdir>/<TF>/beds/<TF>_all.bed (per motif)\n"
	description += "- <outdir>/<TF>/beds/<TF>_<condition>_bound.bed (per motif-condition pair)\n- <outdir>/<TF>/beds/<TF>_<condition>_unbound.bed (per motif-condition pair)\n\n"
	parser.description = format_help_description("BINDetect", description)

	parser._action_groups.pop()	#pop -h
	
	required = parser.add_argument_group('Required arguments')
	required.add_argument('--signals', metavar="<bigwig>", help="Signal per condition (.bigwig format)", nargs="*")
	required.add_argument('--peaks', metavar="<bed>", help="Peaks.bed containing open chromatin regions across all conditions")
	required.add_argument('--motifs', metavar="<motifs>", help="Motif file(s) in pfm/jaspar/meme/transfac format", nargs="*")
	required.add_argument('--genome', metavar="<fasta>", help="Genome .fasta file")

	optargs = parser.add_argument_group('Optional arguments')
	optargs.add_argument('--cond-names', metavar="<name>", nargs="*", help="Names of conditions fitting to --signals (default: prefix of --signals)")
	optargs.add_argument('--peak-header', metavar="<file>", help="File containing the header of --peaks separated by whitespace or newlines (default: peak columns are named \"_additional_<count>\")")
	optargs.add_argument('--naming', metavar="<string>", help="Naming convention for TF output files ('id', 'name', 'name_id', 'id_name') (default: 'name_id')", choices=["id", "name", "name_id", "id_name"], default="name_id")
	optargs.add_argument('--motif-pvalue', metavar="<float>", type=lambda x: restricted_float(x, 0, 1), help="Set p-value threshold for motif scanning (default: 1e-4)", default=0.0001)
	optargs.add_argument('--bound-pvalue', metavar="<float>", type=lambda x: restricted_float(x, 0, 1), help="Set p-value threshold for bound/unbound split (default: 0.001)", default=0.001)
	optargs.add_argument('--cluster-threshold', metavar="<float>", type=lambda x: restricted_float(x, 0, 1), help="Set the clustering threshold. Motifs below this threshold will be assigned to one cluster (default: 0.5)", default=0.5)
	#optargs.add_argument('--volcano-diff-thresh', metavar="<float>", help="", default=0.2)	#not yet implemented
	#optargs.add_argument('--volcano-p-thresh', metavar="<float>", help="", default=0.05)	#not yet implemented

	optargs.add_argument('--pseudo', type=float, metavar="<float>", help="Pseudocount for calculating log2fcs (default: estimated from data)", default=None)
	optargs.add_argument('--time-series', action='store_true', help="Will only compare signals1<->signals2<->signals3 (...) in order of input, and skip all-against-all comparison.")
	optargs.add_argument('--skip-excel', action='store_true', help="Skip creation of excel files - for large datasets, this will speed up BINDetect considerably")
	optargs.add_argument('--output-peaks', metavar="<bed>", help="""Gives the possibility to set the output peak set differently than the input --peaks.
													 				This will limit all analysis to the regions in --output-peaks. 
																	NOTE: --peaks must still be set to the full peak set!""")
	optargs.add_argument('--norm-off', action='store_true', help="Turn off normalization of footprint scores across conditions")

	runargs = parser.add_argument_group("Run arguments")
	runargs.add_argument('--outdir', metavar="<directory>", help="Output directory to place TFBS/plots in (default: bindetect_output)", default="bindetect_output")
	optargs.add_argument('--prefix', metavar="<prefix>", help="Prefix for overview files in --outdir folder (default: bindetect)", default="bindetect")
	runargs.add_argument('--cores', metavar="<int>", type=int, help="Number of cores to use for computation (default: 1)", default=1)
	runargs.add_argument('--split', metavar="<int>", type=int, help="Split of multiprocessing jobs (default: 100)", default=100)
	runargs.add_argument('--debug', action='store_true', help="Creates an additional '_debug.pdf'-file with debug plots")	#creates extra output for debugging
	
	runargs = add_logger_args(runargs)

	return(parser)

#--------------------------------------------------------------------------------------------------------#
def add_tfbscan_arguments(parser):

	parser.formatter_class = lambda prog: argparse.RawDescriptionHelpFormatter(prog, max_help_position=35, width=90)
	description = "Find positions of Transcription Factor Binding Sites (TFBS) in FASTA sequences by scanning with motifs.\n\n" 
	description += "Usage:\nTOBIAS TFBScan --motifs <motifs.txt> --fasta <genome.fa> \n\n"
	description += "By setting --outdir, the output files are:\n- <outdir>/<TF1>.bed\n- <outdir>/<TF2>.bed\n- (...)\n\n"
	description += "By setting --outfile, all TFBS are written to one file (with motif specified in the 4th column of the .bed)."
	parser.description = format_help_description("TFBScan", description)

	parser._action_groups.pop()	#pop -h

	required_arguments = parser.add_argument_group('Required arguments')
	required_arguments.add_argument('-m', '--motifs', metavar="", help='File containing motifs in either MEME, PFM, JASPAR or TRANSFAC format')
	required_arguments.add_argument('-f', '--fasta', metavar="", help='A fasta file of sequences to use for scanning motifs') 	# whole genome file or regions of interest in FASTA format to be scanned with motifs')

	#all other arguments are optional
	optional_arguments = parser.add_argument_group('Optional arguments')
	optional_arguments.add_argument('-r', '--regions', metavar="", help='Subset scanning to regions of interest')
	optional_arguments.add_argument('--outdir', metavar="", help='Output directory for TFBS sites in one file per motif (default: ./tfbscan_output/). NOTE: Select either --outdir or --outfile.', default=None)
	optional_arguments.add_argument('--outfile', metavar="", help='Output file for TFBS sites joined in one bed-file (default: not set). NOTE: Select either --outdir or --outfile.', default=None)

	optional_arguments.add_argument('--naming', metavar="", help="Naming convention for bed-ids and output files ('id', 'name', 'name_id', 'id_name') (default: 'name_id')", choices=["id", "name", "name_id", "id_name"], default="name_id")
	optional_arguments.add_argument('--gc', metavar="", type=lambda x: restricted_float(x,0,1), help='Set the gc content for background regions (default: will be estimated from fasta)')
	optional_arguments.add_argument('--pvalue', metavar="", type=lambda x: restricted_float(x,0,1), help='Set p-value for motif matches (default: 0.0001)', default=0.0001)
	optional_arguments.add_argument('--keep-overlaps', action='store_true', help='Keep overlaps of same motifs (default: overlaps are resolved by keeping best-scoring site)')
	optional_arguments.add_argument('--add-region-columns', action='store_true', help="Add extra information columns (starting from 4th column) from --regions to the output .bed-file(s) (default: off)")

	RUN = parser.add_argument_group('Run arguments')
	RUN.add_argument('--split', metavar="<int>", type=int, help="Split of multiprocessing jobs (default: 100)", default=100)
	RUN.add_argument('--cores', metavar="", type=int, help='Number of cores to use (default: 1)', default=1)
	RUN.add_argument('--debug', action="store_true", help=argparse.SUPPRESS)
	RUN = add_logger_args(optional_arguments)

	return(parser)

#--------------------------------------------------------------------------------------------------------#
def add_formatmotifs_arguments(parser):

	parser.formatter_class = lambda prog: argparse.RawDescriptionHelpFormatter(prog, max_help_position=40, width=90)
	description = ""
	parser.description = format_help_description("FormatMotifs", description) 
	
	parser._action_groups.pop()	#pop -h

	#Required arguments
	required = parser.add_argument_group('Required arguments')
	required.add_argument('--input', metavar="", nargs="*", help="One or more input motif files (required)")			
	required.add_argument('--output', metavar="", help="If task == join, output is the joined output file; if task == split, output is a directory (required)")
	
	additional = parser.add_argument_group('Additional arguments')
	additional.add_argument('--format', metavar="", help="Desired motif output format (pfm, jaspar, meme) (default: \"jaspar\")", choices=["pfm", "jaspar", "meme"], default="jaspar")
	additional.add_argument('--task', metavar="", help="Which task to perform on motif files (join/split) (default: join)", choices=["join", "split"], default="join")
	additional.add_argument('--filter', metavar="", help="File containing list of motif names/ids to filter on. Only motifs fitting entries in filter will be output.")
	additional = add_logger_args(additional)

	return(parser)

#--------------------------------------------------------------------------------------------------------#	
def add_scorebed_arguments(parser):

	parser.formatter_class = lambda prog: argparse.RawDescriptionHelpFormatter(prog, max_help_position=40, width=90)
	description = "ScoreBed is a utility to score .bed-file regions with values from a .bigwig-file. The output is a .bed-file with the bigwig value(s) as extra column(s). Options --position and --math can be used to adjust scoring scheme."
	parser.description = format_help_description("ScoreBed", description)

	parser._action_groups.pop()	#pop -h
	
	#Required arguments
	required = parser.add_argument_group('Required arguments')
	required.add_argument('--bed', metavar="", help="Sites to score (.bed file)")
	required.add_argument('--bigwigs', metavar="", nargs="*",  help="Scores to assign to regions in .bed (.bw file(s))")
	
	#Optional arguments
	optional = parser.add_argument_group('Optional arguments')
	optional.add_argument('--output', metavar="", help="Path to output .bed-file (default: scored sites are written to stdout)") 
	optional.add_argument('--subset', metavar="", help="Subset scoring to .bed regions and set all other sites to --null value (default: all sites in input file will be scored)")
	optional.add_argument('--null', metavar="", help="If --subset is given, which score/label to add to non-scored regions (default: 0)", default="0", type=float)
	optional.add_argument('--position', metavar="", help="Position in sites to score (start/mid/end/full) (default: full)", choices=["mid", "start", "end", "full"], default="full")
	optional.add_argument('--math', metavar="", help="If position == full, choose math to perform on signal (min/max/mean/sum) (default: mean)", choices=["min", "max", "mean", "sum"], default="mean")
	optional = add_logger_args(optional)
	#optional.add_argument('--buffer', metavar="", help="Lines to buffer before writing (default: 10000)", type=int, default=10000)

	return(parser)

#--------------------------------------------------------------------------------------------------------#		
def add_aggregate_arguments(parser):

	parser.formatter_class = lambda prog: argparse.RawDescriptionHelpFormatter(prog, max_help_position=40, width=90)
	description = ""
	parser.description = format_help_description("PlotAggregate", description)

	parser._action_groups.pop()	#pop -h

	IO = parser.add_argument_group('Input / output arguments')
	IO.add_argument('--TFBS', metavar="<bed>", nargs="*", help="TFBS sites (*required)") 						#default is None
	IO.add_argument('--signals', metavar="<bigwig>", nargs="*", help="Signals in bigwig format (*required)")	#default is None
	IO.add_argument('--regions', metavar="<bed>", nargs="*", help="Regions to overlap with TFBS (optional)", default=[])
	IO.add_argument('--whitelist', metavar="<bed>", nargs="*", help="Only plot sites overlapping whitelist (optional)", default=[])
	IO.add_argument('--blacklist', metavar="<bed>", nargs="*", help="Exclude sites overlapping blacklist (optional)", default=[])
	IO.add_argument('--output', metavar="", help="Path to output plot (default: TOBIAS_aggregate.pdf)", default="TOBIAS_aggregate.pdf")
	IO.add_argument('--output-txt', metavar="", help="Path to output file for aggregates in .txt-format (default: None)")

	PLOT = parser.add_argument_group('Plot arguments')
	PLOT.add_argument('--title', metavar="", help="Title of plot (default: \"Aggregated signals\")", default="Aggregated signals")
	PLOT.add_argument('--flank', metavar="", help="Flanking basepairs (+/-) to show in plot (counted from middle of the TFBS) (default: 60)", default=60, type=int)
	PLOT.add_argument('--TFBS-labels', metavar="", help="Labels used for each TFBS file (default: prefix of each --TFBS)", nargs="*")
	PLOT.add_argument('--signal-labels', metavar="", help="Labels used for each signal file (default: prefix of each --signals)", nargs="*")
	PLOT.add_argument('--region-labels', metavar="", help="Labels used for each regions file (default: prefix of each --regions)", nargs="*")
	PLOT.add_argument('--share-y', metavar="", help="Share y-axis range across plots (none/signals/sites/both). Use \"--share-y signals\" if bigwig signals have similar ranges. Use \"--share_y sites\" if sites per bigwig are comparable, but bigwigs themselves aren't comparable (default: none)", choices=["none", "signals", "sites", "both"], default="none")
	
	#Signals / regions
	PLOT.add_argument('--normalize', action='store_true', help="Normalize the aggregate signal(s) to be between 0-1 (default: the true range of values is shown)")
	PLOT.add_argument('--negate', action='store_true', help="Negate overlap with regions")
	PLOT.add_argument('--smooth', metavar="<int>", type=int, help="Smooth output signal by taking the mean of <smooth> bp windows (default: 1 (no smooth)", default=1)
	PLOT.add_argument('--log-transform', help="Log transform the signals before aggregation", action="store_true")
	PLOT.add_argument('--plot-boundaries', help="Plot TFBS boundaries (Note: estimated from first region in each --TFBS)", action='store_true')
	PLOT.add_argument('--signal-on-x', help="Show signals on x-axis and TFBSs on y-axis (default: signal is on y-axis)", action='store_true')
	PLOT.add_argument('--remove-outliers', metavar="<float>", help="Value between 0-1 indicating the percentile of regions to include, e.g. 0.99 to remove the sites with 1%% highest values (default: 1)", type=lambda x: restricted_float(x, 0, 1), default=1)

	RUN = parser.add_argument_group("Run arguments")
	RUN = add_logger_args(RUN)
	
	return(parser)

#--------------------------------------------------------------------------------------------------------#		
def add_plotchanges_arguments(parser):

	parser.formatter_class = lambda prog: argparse.RawDescriptionHelpFormatter(prog, max_help_position=35, width=90)
	description = "PlotChanges is a utility to plot the changes in TF binding across multiple conditions as predicted by TOBIAS BINDetect.\n\n"
	description += "Example usage:\n$ echo CTCF GATA > TFS.txt\n$ TOBIAS PlotChanges --bindetect <bindetect_results.txt> --TFS TFS.txt\n\n"

	parser.description = format_help_description("PlotChanges", description)

	parser._action_groups.pop()	#pop -h

	required_arguments = parser.add_argument_group('Required arguments')
	required_arguments.add_argument('--bindetect', metavar="", help='Bindetect_results.txt file from BINDetect run')
	
	#All other arguments are optional
	optional_arguments = parser.add_argument_group('Optional arguments')
	optional_arguments.add_argument('--TFS', metavar="", help='Text file containing names of TFs to show in plot (one per line)') 
	optional_arguments.add_argument('--output', metavar="", help='Output file for plot (default: bindetect_changes.pdf)', default="bindetect_changes.pdf")
	optional_arguments.add_argument('--conditions', metavar="", help="Ordered list of conditions to show (default: conditions are ordered as within the bindetect file)", nargs="*")
	optional_arguments = add_logger_args(optional_arguments)
	
	return(parser)

#--------------------------------------------------------------------------------------------------------#
def add_heatmap_arguments(parser):

	parser.formatter_class = lambda prog: argparse.RawDescriptionHelpFormatter(prog, max_help_position=40, width=90)
	description = "PlotHeatmap plots a heatmap of signals from bigwig(s) (each row is one site) as well as the aggregate signal across all sites."
	parser.description = format_help_description("PlotHeatmap", description)
	
	parser._action_groups.pop()	#pop -h
	
	IO = parser.add_argument_group('Input / output arguments')
	IO.add_argument('--TFBS', metavar="", nargs="*", action='append', help="TFBS sites per column (*required)")	#if more than one, set to next column
	IO.add_argument('--signals', metavar="", nargs="*", help="Signals in bigwig format (*required)")
	IO.add_argument('--output',  metavar="", help="Output filename (default: TOBIAS_heatmap.pdf)", default="TOBIAS_heatmap.pdf")

	PLOT = parser.add_argument_group('Plot arguments')
	PLOT.add_argument('--plot-boundaries', help="Plot TFBS boundaries", action='store_true')
	PLOT.add_argument('--share-colorbar', help="Share colorbar across all bigwigs (default: estimate colorbar per bigwig)", action='store_true')
	PLOT.add_argument('--flank', metavar="", help="", type=int, default=75)
	
	PLOT.add_argument('--title', metavar="", default="TOBIAS heatmap")
	PLOT.add_argument('--TFBS-labels', metavar="", nargs="*", action='append', help="Labels of TFBS (default: basename of --TFBS)")
	PLOT.add_argument('--signal-labels', metavar="", nargs="*", help="Labels of signals (default: basename of --signals)")
	PLOT.add_argument('--cmap', metavar="", help="Name of colormap to use for heatmaps (default: seismic)", default="seismic")

	PLOT.add_argument('--show-columns', nargs="*", metavar="", type=int, help="Show scores from TFBS column besides heatmap. Set to 0-based python coordinates (for example -1 for last column) (default: None)", default=[])
	PLOT.add_argument('--sort-by', metavar="", help="Columns in .bed to sort heatmap by (default: input .beds are not sorted)", type=int)

	RUN = parser.add_argument_group('Run arguments')
	RUN = add_logger_args(RUN)

	return(parser)


#--------------------------------------------------------------------------------------------------------#
def add_tracks_arguments(parser):

	parser.formatter_class = lambda prog: argparse.RawDescriptionHelpFormatter(prog, max_help_position=40, width=90)
	description = "Plot genomic tracks (such as cutsite or footprint scores) in specific genomic regions.\n"
	description += "This function is a wrapper for the svist4get package (Artyom A. Egorov, \"svist4get: a simple visualization tool for genomic tracks from sequencing experiments\", BMC Bioinformatics, Volume 20, 2019, 113)"
	description += " which allows for automatic creation of multiple plots by using bigwigs/bedfiles."
	parser.description = format_help_description("PlotTracks", description)
	parser._action_groups.pop()		#pop -h

	IO = parser.add_argument_group('Input / output arguments')
	IO.add_argument('--bigwigs', metavar="",  action='append', nargs="*", help="One or more bigwigs to show. Note: All bigwigs within one \"--bigwigs\" argument will be normalized to each other. " + 
																				"It is possible to give multiple \"--bigwigs\" arguments, which will be normalized independently per group (required)", default=[])
	IO.add_argument('--regions', metavar="", help="Genomic regions to plot (required)")
	IO.add_argument('--sites', metavar="", help="Genomic sites to show in plot (optional)")
	IO.add_argument('--highlight', metavar="", help="Regions to highlight in plot (optional)")
	IO.add_argument('--gtf', metavar="", help="GTF file containing genes to show (optional)")

	IO.add_argument('--width', metavar="", help="Width of plot in cm (default 15)", type=float, default=15)
	#IO.add_argument('--height', metavar="")
	IO.add_argument('--colors', metavar="", nargs="*", help="List of specific colors to use for plotting tracks", default=None)
	IO.add_argument('--labels', metavar="", nargs="*", help="Labels for tracks (default: prefix of bigwig)")
	IO.add_argument('--max-transcripts', metavar="", type=int, help="Set a limit on number of transcripts per gene shown in plot (default: 1)", default=1)

	IO.add_argument('--outdir', metavar="", help="Output folder (default: plottracks_output)", default="plottracks_output")
	IO = add_logger_args(IO)
	
	return(parser)

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

#--------------------------------------------------------------------------------------------------------#
def add_maxpos_arguments(parser):

	parser.formatter_class = lambda prog: argparse.RawDescriptionHelpFormatter(prog, max_help_position=40, width=90)
	description = "MaxPos identifies the position of maximum signal (from bigwig) within a given set of .bed-regions. Used to identify peak of signals such as accessibility or footprint scores.\n\n"
	description += "Usage:\nTOBIAS MaxPos --bed <regions.bed> --bigwig <signal.bw>\n\n"
	description += "The output is a 4-column .bed-file (default output is stdout, but you can use --output to set a specific output file).\n"
	parser.description = format_help_description("MaxPos", description)

	parser._action_groups.pop()	#pop -h
	
	#Required arguments
	required = parser.add_argument_group('Required arguments')
	required.add_argument('--bed', metavar="", help="Regions to search for max position within")
	required.add_argument('--bigwig', metavar="", help="Scores used to identify maximum value")

	#Optional arguments
	optional = parser.add_argument_group('Optional arguments')
	optional.add_argument('--output', metavar="", help="Path to output .bed-file (default: scored sites are written to stdout)") 
	optional.add_argument('--invert', help="Find minimum position instead of maximum position", action='store_true', default=False)

	return(parser)

#--------------------------------------------------------------------------------------------------------#
def add_subsample_arguments(parser):

	parser.formatter_class=lambda prog: argparse.RawDescriptionHelpFormatter(prog, max_help_position=40, width=100)
	description = "Subsample a .bam-file into random subsets (of percentage size --start to --end with step --step).\n"
	description = "NOTE: Requires 'samtools' command available at runtime."
	parser.description = format_help_description("SubsampleBam", description)

	parser._action_groups.pop()	#pop -h

	#Required arguments
	args = parser.add_argument_group('Input arguments')
	args.add_argument('--bam', metavar="", help="Path to .bam-file")
	args.add_argument('--no_rand', metavar="", type=int, help="Number of randomizations (per step) (default: 3)", default=3)
	args.add_argument('--start', metavar="", type=int, help="Start of percent subsample (default: 0)", default=0)
	args.add_argument('--end', metavar="", type=int, help="End of percent subsample (default: 100)", default=100)
	args.add_argument('--step', metavar="", type=int, help="Step between --start and --end (default: 10)", default=10)
	args.add_argument('--cores', metavar="", type=int, help="Cores for multiprocessing (default: 1)", default=1)
	args.add_argument('--outdir', metavar="", help="Output directory (default: subsamplebam_output)", default="subsamplebam_output")
	args.add_argument('--prefix', metavar="", help="Prefix for output files (default: prefix of .bam)")
	args.add_argument('--force', action="store_true", help="Force creation of subsampled .bam-files (default: only create if not existing)")
	args = add_logger_args(args)

	return(parser)

#--------------------------------------------------------------------------------------------------------#
def add_network_arguments(parser):

	parser.formatter_class = lambda prog: argparse.RawDescriptionHelpFormatter(prog, max_help_position=40, width=90)
	description = "Creates a TF-TF gene regulation network from annotated transcription factor binding sites" 
	parser.description = format_help_description("CreateNetwork", description)

	parser._action_groups.pop()	#pop -h
	
	#Required arguments
	required = parser.add_argument_group('Required arguments')
	required.add_argument('--TFBS', metavar="", help="File(s) containing TFBS (with annotation) to create network from", nargs="*")
	required.add_argument('--origin', metavar="", help="File containing mapping of TF <-> origin gene")

	#Optional arguments
	optional = parser.add_argument_group('Optional arguments')
	optional.add_argument('--start', metavar="", help="Name of node to start in (default: all nodes)")
	optional.add_argument('--max-len', metavar="", help="Maximum number of nodes in paths through graph (default: 4)", type=int, default=4)
	#optional.add_argument('--unique', action='store_true', help="Only include edges once (default: edges can occur multiple times in case of multiple binding sites)")
	
	runargs = parser.add_argument_group("Run arguments")
	runargs.add_argument('--outdir', metavar="", help="Path to output directory (default: createnetwork_output)", default="createnetwork_output") 
	runargs = add_logger_args(runargs)

	return(parser)

#--------------------------------------------------------------------------------------------------------#
def add_log2table_arguments(parser):

	parser.formatter_class = lambda prog: argparse.RawDescriptionHelpFormatter(prog, max_help_position=40, width=90)
	description = "Log2Table creates tables of footprint depth (FPD) and aggregate correlations from the PlotAggregate logfiles." 
	parser.description = format_help_description("Log2Table", description)

	parser._action_groups.pop()	#pop -h
	
	#Required arguments
	required = parser.add_argument_group('Required arguments')
	required.add_argument('--logfiles', nargs="*", metavar="", help="Logfiles from PlotAggregate")
	required.add_argument('--outdir', metavar="", help="Output directory for tables (default: current dir)", default=".")
	required.add_argument('--prefix', metavar="", help="Prefix of output files", default="aggregate")

	return(parser)


#----------------------------------------------------------------------------------------#
def add_filterfragments_arguments(parser):

	parser.formatter_class = lambda prog: argparse.RawDescriptionHelpFormatter(prog, max_help_position=40, width=90)
	description = "FilterFragments can filter out fragments from a .bam-file of paired-end reads based on the overlap with regions in the given .bed-file."
	#description += "Usage: "
	parser.description = format_help_description("FilterFragments", description)

	parser._action_groups.pop()	#pop -h

	IO = parser.add_argument_group('Input / output arguments')
	IO.add_argument('--bam', metavar="", help=".bam-file to filter")
	IO.add_argument('--regions', metavar="", help=".bed-file containing regions to filter fragments from")
	IO.add_argument('--mode', help="Mode 1: Remove fragment if both reads overlap .bed-regions. Mode 2: Remove whole fragment if one read is overlapping .bed-regions (default: 1)", choices=[1,2], default=1)
	IO.add_argument('--output', metavar="", help="Path to the filtered .bam-file (default: <prefix of --bam>_filtered.bam)")
	IO.add_argument('--threads', metavar="", type=int, help="Number of threads used for decompressing/compressing bam (default: 10)", default=10)

	IO = add_logger_args(IO)

	return(parser)

#--------------------------------------------------------------------------------------------------------#
def add_motifclust_arguments(parser):
	"""Parsing arguments using argparse

	Returns
	-------
	ArgumentParser
		an object containing all parameters given.
	"""

	parser.formatter_class = lambda prog: argparse.RawDescriptionHelpFormatter(prog, max_help_position=40, width=90)
	description = "Cluster motifs based on similarity and create one consensus motif per cluster.\n\n"
	description += "Usage:\nTOBIAS ClusterMotifs --motifs <motifs.jaspar>\n\n"
	description += "NOTE: This tool requres the python package 'gimmemotifs'."
	parser.description = format_help_description("ClusterMotifs", description) 

	parser._action_groups.pop() #pop -h

	required = parser.add_argument_group('Required arguments')
	required.add_argument("-m", "--motifs", required=True, help="One or more motif files to compare and cluster", nargs="*", metavar="")
	
	optional = parser.add_argument_group('Optional arguments')
	optional.add_argument("-t", "--threshold", metavar="", help="Clustering threshold (Default = 0.3)", type=float, default=0.3)  
	optional.add_argument('--dist_method', metavar="", help="Method for calculating similarity between motifs (default: pcc)", choices=["pcc", "seqcor", "wic", "akl"], default="pcc") 
	optional.add_argument('--clust_method', metavar="", help="Method for clustering (See: https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html)", default="average", choices=["single","complete","average","weighted","centroid","median","ward"])
	optional.add_argument("-a", "--cons_format", metavar="", choices= ['meme', 'pfm', 'jaspar'], help="Format of consensus motif file [‘meme’, 'pfm', 'jaspar'] (Default: jaspar)", default="jaspar")
	optional.add_argument("-p", "--prefix", metavar="", help="Output prefix (default: 'motif_comparison')", default="motif_comparison")
	optional.add_argument("-o", "--outdir", metavar="", help="Output directory (default: 'clustermotifs_output')", default="clustermotifs_output")
	optional = add_logger_args(optional)

	visualisation = parser.add_argument_group('Visualisation arguments')
	visualisation.add_argument("-e", "--type", dest="type", choices= ['png', 'pdf', 'jpg'], help="Plot file type [png, pdf, jpg] (Default: pdf)", default="pdf")
	#visualisation.add_argument("-x", "--width", dest="width", help="Width of Heatmap (Default: autoscaling)", type=int)
	#visualisation.add_argument("-y", "--height", dest="height", help="Height of Heatmap (Default: autoscaling)", type=int)
	visualisation.add_argument("-d", "--dpi", dest="dpi", help="Dpi for plots (Default: 100)", type=int, default=100)
	#visualisation.add_argument("-ncc", "--no_col_clust", dest="ncc", help="No column clustering", action="store_true")
	#visualisation.add_argument("-nrc", "--no_row_clust", dest="nrc", help="No row clustering", action="store_true")
	visualisation.add_argument("-c", "--color_palette", dest="color", help="Color palette (All possible paletts: https://python-graph-gallery.com/197-available-color-palettes-with-matplotlib/. Add '_r' to reverse palette.)", default="YlOrRd_r")

	return parser

#--------------------------------------------------------------------------------------------------------#

def add_downloaddata_arguments(parser):

	parser.formatter_class = lambda prog: argparse.RawDescriptionHelpFormatter(prog, max_help_position=40, width=90)
	description = "Download test data for the TOBIAS commandline examples"
	parser.description = format_help_description("DownloadData", description) 

	parser._action_groups.pop()	#pop -h

	arguments = parser.add_argument_group('Arguments')
	arguments.add_argument('--endpoint', metavar="", help="Link to the s3 server (default: The loosolab s3 server)", default="https://s3.mpi-bn.mpg.de")
	arguments.add_argument('--bucket', metavar="", help="Name of bucket to download (default: data-tobias-2020)", default="data-tobias-2020")
	#arguments.add_argument('--target', metavar="", help="Name of directory to save files to (default: name of bucket)")
	arguments.add_argument('--patterns', metavar="", help="List of patterns for files to download e.g. '*.txt' (default: *)", default="*")
	arguments.add_argument('--username', metavar="", help="Username for endpoint (default: None set)")
	arguments.add_argument("--key", metavar="", help="Access key for endpoint (default: None set)")
	arguments.add_argument("--yaml", metavar="", help="Set the endpoint/bucket/access information through a config file in .yaml format (NOTE: overwrites commandline arguments)")
	arguments.add_argument("--force", action="store_true", help="Force download of already existing files (default: warn and skip)")

	arguments = add_logger_args(arguments)

	return parser
