#!/usr/bin/env python

"""
PlotTracks: Plot genomic tracks such as footprint scores in specific regions. 
Wraps the commandline tool svist4get to automatically plot scores/sites for each region in the input '--regions' .bed-file.

@author: Mette Bentsen
@contact: mette.bentsen (at) mpi-bn.mpg.de
@license: MIT
"""

import os
import sys
import argparse
import numpy as np

from PyPDF2 import PdfFileMerger, PdfFileReader
import matplotlib
from matplotlib import textpath
from matplotlib.font_manager import findfont, FontProperties
import matplotlib.colors as mcolors
from shutil import which

#Bio-stuff
import pyBigWig
import pybedtools as pb
import svist4get as sv4g

from tobias.parsers import add_tracks_arguments
from tobias.utils.utilities import *

#--------------------------------------------------------------------------------------------------------#
def svist4get_defaults():
	""" Config parameters used to visualize TOBIAS tracks """

	c = {}	#initalize config dict

	#General options
	c["png_dpi"] = 600
	c["revcomp_transform"] = 0
	c["hide_introns"] = 0
	c["show_title"] = 1
	c["font_size_title"] = 7.5
	c["show_aa_seq_track"] = 1
	c["show_genomic_axis"] = 1
	c["show_genomic_axis_tics"] = 1 
	c["show_nt_seq_track"] = 0
	c["transcript_label"] = "name"

	c["triplet_code"] = ""
	c["output_filename"] = "svist4get_out.pdf"
	c["num_of_digits_after_the_dot"] = 2

	#Get fonts
	data_dir = os.path.join(os.path.dirname(sv4g.__file__), "svist4get_data")	#data dir of the svist4get module
	if os.path.exists(data_dir):
		c["mono_font"] = os.path.join(data_dir, "fonts", "iosevka-regular.ttf")
		c["regular_font"] = os.path.join(data_dir, "fonts", "Lato-Regular.ttf")
	else:
		c["mono_font"] = findfont(FontProperties(family=['monospace']))
		c["regular_font"] = findfont(FontProperties(family=['sans-serif']))

	#Page parameters
	c["page_width"] = 8
	c["png_dpi"] = 600
	c["gap"] = 0.15
	c["margin"] = 0.5

	#Vertical grid
	c["vgrid_step"] = 0.5
	c["c_vgrid_alpha"] = 0		#grid is invisible
	c["c_vgrid"] = "blue"

	#Highlight settings
	c["c_fill_hframe"] = "grey"
	c["c_fill_hframe_alpha"] = 0.4
	c["c_stroke_hframe"] = "grey"
	c["c_stroke_hframe_alpha"] = 0
	c["c_text_hframe"] = "black"
	c["hframe_font_size"] = 5

	#Bedgraph settings
	c["bedgraph_track_height"] = 0.75
	c["bedgraph_column_min_width"] = 0.0004 	#force 1nt resolution
	c["bedgraph_bar"] = "none"					#force 1nt resolution
	c["bedgraph_tics_font_size"] = 5    
	c["bedgraph_label_font_size"] = 7  
	c["c_bedgraph_label_alpha"] = 0.65
	c["c_bedgraph_tracks"] = ["purple","green","blue","pink","yellow","orange","brown","red"] 	#default, can be overwritten
	c["c_bedgraph_alpha"] = 1
	c['bedgraph_label_position'] = "right"
	c["bedgraph_axis_tics"] = "auto"
	c["bedgraph_axis_tics_step"] = 0   
	c["bedgraph_upper_limit"] = "auto"
	c["bedgraph_lower_limit"] = "auto"

	#Nucleotide sequence track 
	c["c_A"] = "blue"
	c["c_A_alpha"] = 0.8
	c["c_G"] = "orange"
	c["c_G_alpha"] = 0.8
	c["c_C"] = "yellow"
	c["c_C_alpha"] = 0.8
	c["c_T"] = "purple"
	c["c_T_alpha"] = 0.8

	#Genomic axis tics
	c["gaxis_label_font_size"] = 7
	c["gaxis_tics_font_size"] = 6
	c["gaxis_tics_step"] = 0

	#Transcript structure
	c["stroke_width_CDS"] = 0.14
	c["stroke_width_exon"] = 0.2
	c["transcript_id_font_size"] = 5.5   
	c["arrow_height"] = 0.05
	c["c_marks"] = "deepred"
	c["transcript_label_style"] = "name"

	#;[GENOMIC REGIONS TRACK SETTINGS]
	c["regions_line_width"] = 0.1
	c["font_size_regions_label"] = 6
	c["c_regions"] = "black"
	c["c_regions_alpha"] = 0.7

	return(c)

def write_out_config(config, outfile):
	""" Config is a dictionary containing values to write to outfile in 'key = value' format """
	with open(outfile, "w") as f:
		for key in config:
			if isinstance(config[key], tuple) or isinstance(config[key], list):
				f.write("{0} = {1}\n".format(key, ",".join([str(element) for element in config[key]])))
			else: 
				f.write("{0} = {1}\n".format(key, config[key]))

def is_executable(tool):
	"""Checks whether 'tool' is available on the system PATH and can be executed."""

	return which(tool) is not None #True if tool exists

#--------------------------------------------------------------------------------------------------------#
def run_tracks(args):

	#Check if input files exist and create outdir
	check_required(args, ["bigwigs", "regions"])
	check_files([args.bigwigs, args.regions, args.sites, args.highlight, args.gtf], "r")
	make_directory(args.outdir)

	#Setup logger
	logger = TobiasLogger("PlotTracks", args.verbosity)
	logger.begin()

	parser = add_tracks_arguments(argparse.ArgumentParser())
	logger.arguments_overview(parser, args)

	############# Check that dependencies are available #############

	if is_executable("gs") == False:
		logger.error("gs (Ghostscript) is not available on PATH but is needed for PlotTracks. Please install gs to continue.")
		sys.exit(1)

	################# Setup custom config file ####################

	#Create custom palette
	palette_dict = mcolors.CSS4_COLORS
	palette_str = "\n".join("{0} = {1}".format(tup[0], tup[1]) for tup in palette_dict.items())   #tup is color_name, color_hash
	palette_path = os.path.abspath(os.path.join(args.outdir, "palette.txt"))
	f = open(palette_path, "w")
	f.write(palette_str)
	f.close()

	#Setup config
	c = svist4get_defaults()
	c["palette"] = palette_path
	c["triplet_code"] = "empty.txt"
	c["page_width"] = args.width

	#If sites/highlight/gtf are None; set to empty file
	emptyfile = os.path.abspath(os.path.join(args.outdir, "empty.txt"))
	open(emptyfile, "w").close()
	
	args.regions = emptyfile if args.regions is None else args.regions
	args.sites = emptyfile if args.sites is None else args.sites
	args.highlight = emptyfile if args.highlight is None else args.highlight
	args.gtf = emptyfile if args.gtf is None else args.gtf

	#Overwrite config with commandline parameters
	c["c_bedgraph_tracks"] = args.colors if args.colors is not None else c["c_bedgraph_tracks"]
	
	#Check whether colors are in color dict:
	for color in c["c_bedgraph_tracks"]:
		if color not in palette_dict:
			logger.error("Color {0} is unknown to matplotlib - please choose another color.".format(color))
			sys.exit()

	#Duplicate colors as needed
	n_bigwigs = len(sum(args.bigwigs, []))
	c["c_bedgraph_tracks"] = c["c_bedgraph_tracks"] * int(np.ceil(n_bigwigs/len(c["c_bedgraph_tracks"]))) + [c["c_bedgraph_tracks"][0]]
	c["c_bedgraph_tracks"] += c[["c_bedgraph_tracks"][0]]	#add extra color to ensure that colors are being read as a list

	#Write out custom config
	custom_config = os.path.join(args.outdir, "custom_config.cfg")
	write_out_config(c, custom_config)


	################# Get input data ready ####################

	### Unpack bigwigs per group ###
	bigwigs = {}
	bigwig_names = []
	for group in range(len(args.bigwigs)):
		for bw in args.bigwigs[group]:

			name = os.path.splitext(os.path.basename(bw))[0]
			bigwig_names.append(name)
			bigwigs[name] = {"path": bw,
							 "group":group,
							 "pybw":pyBigWig.open(bw)} 

	#Load sites/highlight
	pb_sites = pb.BedTool(args.sites)
	pb_highlight = pb.BedTool(args.highlight)

	#### Go through each region and plot tracks ####
	output_plots = []
	with open(args.regions) as f:
		for line in f:
			if len(line.strip()) == 0:	#empty line
				continue

			columns = line.rstrip().split()
			prefix = "-".join(columns[:3])
			logger.info("Plotting region: {0}".format(columns[:3]))

			chrom, start, end = columns[0], int(columns[1]), int(columns[2])
			
			f_prefix = os.path.join(args.outdir, prefix)
			make_directory(os.path.join(args.outdir, prefix))

			#Extract bigwig signals for this region
			max_vals = {}
			min_vals = {}
			for bw in bigwig_names:
				logger.debug("- Reading signal from {0}".format(bw))
				intervals = bigwigs[bw]["pybw"].intervals(chrom, start, end)
				
				if intervals is not None:
					bedgraph = "\n".join([chrom + "\t" + "\t".join([str(element) for element in interval]) for interval in intervals])
					max_vals[bw] = np.ceil(max([tup[-1] for tup in intervals]) * 10) / 10
					min_vals[bw] = np.floor(min([tup[-1] for tup in intervals]) * 10) / 10
				else:
					bedgraph = ""
					min_vals[bw] = 0
					max_vals[bw] = 0

				f = open(os.path.join(args.outdir, prefix, bw + ".bedgraph"), "w")
				f.write(bedgraph)
				f.close()

			#Establish max per bigwig group
			logger.debug("- Establish min/max per bigwig group")
			group_max = {}
			group_min = {}
			for bw in bigwig_names:
				group = bigwigs[bw]["group"]

				group_max[group] = max([group_max.get(group, max_vals[bw]), max_vals[bw]])
				group_max[group] = round(group_max[group] + group_max[group]*0.1, 2)	#add 10% to prevent tracks from being over the borders

				group_min[group] = round(min([group_min.get(group, min_vals[bw]), min_vals[bw]]), 2)

			logger.debug("  > Group maximums: {0}".format(group_max))
			logger.debug("  > Group minimums: {0}".format(group_min))

			#Extract sites to highlight
			logger.debug("- Extract sites to show/highlight")
			pb_region = pb.BedTool(line, from_string=True)
			sites_in_region = pb_sites.intersect(pb_region, u=True).sort().saveas(os.path.join(f_prefix, "sites.bed"))
			highlight_in_region = pb_highlight.intersect(pb_region, u=True).saveas(os.path.join(f_prefix, "highlight.bed"))	
			#logger.debug("Sites in region: {0}".format(sites_in_region))

			#----------- Setup svist4get -----------# 
			#Call commandline svist4get
			logger.debug("Setting up svist4get call")

			pa = sv4g.manager.Parameters()
			pa.initialize(os.path.join(args.outdir, "custom_config.cfg"))

			#Setup parameters specific for this region
			pa.config['window'] = [chrom, start, end]
			pa.config["output_filename"] = os.path.join(args.outdir, prefix, prefix) 	#without suffix
			pa.config["gtf_file"] = args.gtf

			pa.config["bedgraph"] = [os.path.join(f_prefix, bw + ".bedgraph") for bw in bigwig_names]
			pa.config['bedgraph_label'] = bigwig_names

			pa.config['genomic_intervals'] = [["{0}-{1}".format(site[1], site[2])] for site in sites_in_region]
			pa.config['genomic_intervals_label'] = [site[3] for site in sites_in_region]
			pa.config['highlight_frame'] = [[highlight[1], highlight[2], ""] for highlight in highlight_in_region]
			pa.config['image_title'] = ""

			pa.config['bedgraph_upper_limit'] = [group_max[bigwigs[bw]["group"]] for bw in bigwig_names]
			pa.config['bedgraph_lower_limit'] = [group_min[bigwigs[bw]["group"]] for bw in bigwig_names]

			# gtf data extraction
			gtf = sv4g.data_processing.Gtf_helper(pa.config['gtf_file'])
			transcripts = gtf.extract_transcripts_from_widnow(*pa.config['window'])
			data_from_gtf = (gtf.extract_data_about_transcripts(transcripts))
			
			#Add gene_id to the annotation dicts with no annotation; hack for error
			for anno in data_from_gtf['for_transcript']:
				anno["gene_id"] = "" if "gene_id" not in anno else anno["gene_id"]
				anno["exons"] = "{0}-{1}".format(anno["start_of_transcript"], anno["end_of_transcript"]) if anno["exons"] == "" else anno["exons"]
				anno["CDS"] = "{0}-{1}".format(anno["start_of_transcript"], anno["end_of_transcript"]) if anno["CDS"] == "" else anno["CDS"]

			#Set maximum on transcripts (per gene) to minimize showing transcript 1-100 for the same gene
			final_data_from_gtf = data_from_gtf.copy()
			final_data_from_gtf["for_transcript"] = []
			transcripts_per_gene = {}
			for transcript in data_from_gtf["for_transcript"]:
				transcripts_per_gene[transcript["gene_id"]] = transcripts_per_gene.get(transcript["gene_id"], 0) + 1
				if transcripts_per_gene[transcript["gene_id"]] <= args.max_transcripts:
					final_data_from_gtf["for_transcript"].append(transcript)

			logger.debug("data_from_gtf: {0}".format(final_data_from_gtf))
			pa.add_gtf_data(final_data_from_gtf)	#after subset using args.max_transcripts

			#Write out specific config for future use
			region_config = os.path.join(f_prefix, "config.cfg")
			write_out_config(pa.config, region_config)

			#Collect tracks
			tracks = []

			tracks += sv4g.manager.Title_tracks_maker(pa).create_tracks()
			tracks += sv4g.manager.Axis_tics_tracks_maker(pa).create_tracks()
			tracks += sv4g.manager.Transcript_struct_tracks_maker(pa).create_tracks()
			tracks += sv4g.manager.Bedgraph_tracks_maker(pa).create_tracks()

			#Handle plotting of genomic intervals such as TFBS
			n_intervals = len(pa.config['genomic_intervals'])
			logger.debug("n genomic intervals: {0}".format(n_intervals))
			if len(pa.config['genomic_intervals']) > 0:
				
				genomic_intervals = pa.config['genomic_intervals'][:]
				genomic_intervals_label = pa.config['genomic_intervals_label'][:]

				bp_per_width = (end - start)/pa.config["page_width"]	#in cm
				width_per_point = 2.1/50	#2.5cm per 50 points rough estimation
				
				label_bp = []
				extended_intervals = []
				for interval, label in zip(pa.config['genomic_intervals'], pa.config['genomic_intervals_label']):

					t = matplotlib.textpath.TextPath((0,0), label, size=pa.config["font_size_regions_label"])
					bb = t.get_extents()
					text_w_points = bb.width

					#Q: how many bp does this amount to in plot?
					text_bp = int(text_w_points * width_per_point * bp_per_width)
					label_bp.append(text_bp)
					
					#Extend each site with slop
					interval_start, interval_end = (int(pos) for pos in interval[0].split("-"))
					slop_bp = int((text_bp - (interval_end - interval_start))/2)
					extended_intervals.append((interval_start-slop_bp, interval_end+slop_bp))

				#Find maximum number of overlapping transcripts
				positions = {}
				for (interval_start, interval_end) in extended_intervals:
					for i in range(interval_start, interval_end):
						positions[i] = positions.get(i, 0) + 1
				n_interval_tracks = max(positions.values())
				logger.debug("Number of interval tracks: {0}".format(n_interval_tracks))

				#Create tracks of separate genomic intervals
				for track_n in range(n_interval_tracks):
					selection = list(range(track_n, n_intervals, n_interval_tracks))

					pa.config['genomic_intervals'] = [genomic_intervals[i] for i in selection]
					pa.config['genomic_intervals_label'] = [genomic_intervals_label[i] for i in selection]
					tracks += sv4g.manager.Genomic_intervals_tracks_maker(pa).create_tracks()	

			#Draw to pdf
			sv4g.manager.Image(tracks, pa).draw()
			sv4g.methods.pdf_page_to_png(pa)	#pdf to png
			output_plots.append(pa.config["output_filename"])

	#Done plotting all regions; join to one pdf
	merger = PdfFileMerger(strict=False)
	pdf_filenames = [prefix + ".pdf" for prefix in output_plots]
	for pdf in pdf_filenames:
		if os.stat(pdf).st_size != 0:	#only join files containing plots
			merger.append(PdfFileReader(pdf))
	merger.write(os.path.join(args.outdir, "all_plots.pdf"))

	#End PlotTracks
	logger.end()


#--------------------------------------------------------------------------------------------------------#
if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser = add_tracks_arguments(parser)
	args = parser.parse_args()

	if len(sys.argv[1:]) == 0:
		parser.print_help()
		sys.exit()
		
	run_tracks(args)
