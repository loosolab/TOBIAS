"""
TFBScan.py produces the data to be used to join the footprint and motif information across genome

@author: Anastasiia Petrova and Mette Bentsen
@contact: anastasiia.petrova(at)mpi-bn.mpg.de and mette.bentsen(at)mpi-bn.mpg.de

"""

import argparse
import sys
import os
import re
from datetime import datetime
import multiprocessing as mp
import itertools

import pysam

from tobias.utils.utilities import *
from tobias.utils.regions import *
from tobias.utils.sequences import * 
from tobias.utils.motifs import *


def restricted_float(f, f_min, f_max):
    f = float(f)
    if f < f_min or f > f_max:
        raise argparse.ArgumentTypeError("{0} not in range [0.0, 1.0]".format(f))
    return f


#----------------------------------------------------------------------------------------------------------#
def add_tfbscan_arguments(parser):

	parser.formatter_class = lambda prog: argparse.RawDescriptionHelpFormatter(prog, max_help_position=35, width=90)
	description = "Find positions of Transcription Factor Binding Sites (TFBS) in FASTA sequences by scanning with motifs.\n\n" 
	description += "Usage:\nTOBIAS TFBScan --motifs <motifs.txt> --fasta <genome.fa> \n\n"
	description += "Output files:\n- <outdir>/<TF1>.bed\n- <outdir>/<TF2>.bed\n- (...)\n\n"
	parser.description = format_help_description("TFBScan", description)

	parser._action_groups.pop()	#pop -h

	required_arguments = parser.add_argument_group('Required arguments')
	required_arguments.add_argument('-m', '--motifs', metavar="", help='File containing motifs in either MEME, PFM or JASPAR format')
	required_arguments.add_argument('-f', '--fasta', metavar="", help='A fasta file of sequences to use for scanning motifs') 	# whole genome file or regions of interest in FASTA format to be scanned with motifs')

	#all other arguments are optional
	optional_arguments = parser.add_argument_group('Optional arguments')
	optional_arguments.add_argument('-o', '--outdir', metavar="", help='Output directory (default: ./motifscan_output/)', default="motifscan_output")
	optional_arguments.add_argument('-r', '--regions', metavar="", help='Subset scanning to regions of interest')
	optional_arguments.add_argument('--naming', metavar="", help="Naming convention for bed-ids and output files ('id', 'name', 'name_id', 'id_name') (default: 'name_id')", choices=["id", "name", "name_id", "id_name"], default="name_id")
	optional_arguments.add_argument('--gc', metavar="", type=lambda x: restricted_float(x,0,1), help='Set the gc content for background regions (default: will be estimated from fasta)')
	optional_arguments.add_argument('--pvalue', metavar="", type=lambda x: restricted_float(x,0,1), help='Set p-value for motif matches (default: 0.0001)', default=0.0001)
	##optional_arguments.add_argument('--tool', metavar="", choices=['FIMO', 'MOODS'], help='Choose the tool to perform motif scanning with (default: MOODS)', default='MOODS')
	optional_arguments.add_argument('--keep_overlaps', action='store_true', help='Keep overlaps of same motifs (default: overlaps are resolved by keeping best-scoring site)')
	
	optional_arguments.add_argument('--split', metavar="<int>", type=int, help="Split of multiprocessing jobs (default: 100)", default=100)
	optional_arguments.add_argument('--cores', metavar="", type=int, help='Number of cores to use (default: 1)', default=1)
	optional_arguments.add_argument('--log', metavar="", help="Path to logfile (default: writes to stdout)")
	optional_arguments.add_argument('--debug', action="store_true", help=argparse.SUPPRESS)


	return(parser)

#----------------------------------------------------------------------------------------------------------#

def motif_scanning(regions, args, motifs_obj):
	""" motifs_obj is a MotifList object """

	motifs_obj.setup_moods_scanner()	#setup scanner in object
	fasta_obj = pysam.FastaFile(args.fasta)
	qs = args.qs

	#Scan motifs_obj against sequences per region
	all_TFBS = {name:RegionList() for name in [motif.name for motif in motifs_obj]}

	for region in regions:

		seq = fasta_obj.fetch(region.chrom, region.start, region.end)
		region_TFBS = motifs_obj.scan_sequence(seq, region)		#Scan sequence

		#Split to single TFs
		for TFBS in region_TFBS:
			all_TFBS[TFBS.name].append(TFBS)

	#Write to queue
	for name in all_TFBS:
		bed_content = all_TFBS[name].as_bed()	#string 
		qs[name].put((name, bed_content))		#send to writers

	return(1)


def process_TFBS(file, args):
	""" Process TFBS bedfile - remove duplicates, overlaps and sort """

	outfile = file.replace(".tmp", ".bed")

	#Read sites to regionList
	TFBS = RegionList().from_bed(file)

	#Remove overlaps
	if args.keep_overlaps == False:
		TFBS = TFBS.resolve_overlaps()		#automatically removes duplicates as well
	else:
		TFBS = TFBS.remove_duplicates()		#

	#If chrom looks like subset -> convert to true coordinates 
	for site in TFBS:

		m = re.match("(.+)\:([0-9]+)\-([0-9]+)", site.chrom)
		if m:
			reg_chrom, reg_start, reg_end = m.group(1), m.group(2), m.group(3)

			site.chrom = reg_chrom
			site.start = int(reg_start) + site.start
			site.end = int(reg_start) + site.end
			site.update()

	#Write out to bed
	TFBS.loc_sort()
	TFBS.write_bed(outfile)
	
	#Remove the tmp file
	if not args.debug:
		try:
			os.remove(file)
		except:
			print("Error removing file {0}".format(file))

	return(1)


#----------------------------------------------------------------------------------------------------------#
def run_tfbscan(args):

	begin_time = datetime.now()

	###### Check input arguments ######
	check_required(args, ["motifs", "fasta"])				#Check input arguments
	check_files([args.motifs, args.fasta, args.regions]) 	#Check if files exist
	make_directory(args.outdir)							#Check and create output directory

	###### Create logger and write argument overview ######
	logger = create_logger(2, args.log) 

	logger.comment("#TOBIAS TFBScan (run started {0})\n".format(begin_time))
	logger.comment("#Command line call: {0}\n".format(" ".join(sys.argv)))
	
	parser = add_tfbscan_arguments(argparse.ArgumentParser())
	logger.comment(arguments_overview(parser, args))	

	######## Read sequences from file and estimate background gc ########
	logger.info("Reading sequences from fasta")

	fastafile = pysam.FastaFile(args.fasta)
	fasta_chrom_info = dict(zip(fastafile.references, fastafile.lengths))
	fastafile.close()
	logger.info("- Found {0} sequences".format(len(fasta_chrom_info)))
	
	#Create regions available in fasta 	
	logger.info("Setting up regions")
	fasta_regions = RegionList([OneRegion([header, 0, fasta_chrom_info[header]]) for header in fasta_chrom_info])

	#If subset, setup regions
	if args.regions:
		regions = RegionList().from_bed(args.regions)

	else:	#set up regions from fasta references
		regions = fasta_regions

		#Subset regions to maximum length and extend to overlap at junctions
		regions = regions.apply_method(OneRegion.split_region, 1000000)
		regions = regions.apply_method(OneRegion.extend_reg, 50)
		regions = regions.apply_method(OneRegion.check_boundary, fasta_chrom_info, "cut")
		#logger.info("Split regions: {0} ".format(len(regions)))

	logger.info("- Total of {0} regions (after splitting)".format(len(regions)))

	#Background gc
	if args.gc == None:
		logger.info("Estimating GC content from fasta")
		args.gc = get_gc_content(regions, args.fasta)
	
	bg = np.array([(1-args.gc)/2.0, args.gc/2.0, args.gc/2.0, (1-args.gc)/2.0])

	#Split regions
	region_chunks = regions.chunks(args.split)
	

	#################### Read motifs from file ####################
	logger.info("Reading motifs from file")

	motif_content = open(args.motifs).read()
	converted_content = convert_motif(motif_content, "pfm")
	motif_list = pfm_to_motifs(converted_content) 			#List of OneMotif objects

	logger.info("- Found {0} motifs".format(len(motif_list)))
	
	logger.debug("Getting motifs ready")
	motif_list.bg = bg
	motif_names = [motif.name for motif in motif_list]
	motif_list.extend([motif.get_reverse() for motif in motif_list])
	for motif in motif_list:	#now with reverse motifs as well
		motif.set_name(args.naming)
		motif.name = filafy(motif.name)	#remove ()/: etc. which will create problems in filenames
		motif.bg = bg
		motif.get_pssm()
	
	motif_names = list(set([motif.name for motif in motif_list]))


	pool = mp.Pool(processes=args.cores)
	outlist = pool.starmap(OneMotif.get_threshold, itertools.product(motif_list, [args.pvalue])) 
	motif_list = MotifList(outlist)	
	pool.close()
	pool.join()


	#################### Find TFBS in regions #####################

	logger.comment("")
	logger.critical("Scanning for TFBS with all motifs")

	manager = mp.Manager()

	writer_cores = max(1,int(args.cores*0.1))
	worker_cores = max(1,int(args.cores*0.9))

	worker_pool = mp.Pool(processes=worker_cores, maxtasksperchild=1)
	writer_pool = mp.Pool(processes=writer_cores)

	#Setup bed-writers 
	temp_files = []
	qs = {}
	TF_names_chunks = [motif_names[i::writer_cores] for i in range(writer_cores)]
	for TF_names_sub in TF_names_chunks:
		logger.debug("Creating writer queue for {0}".format(TF_names_sub))
		files = [os.path.join(args.outdir, TF + ".tmp") for TF in TF_names_sub]
		temp_files.extend(files)

		q = manager.Queue()
		#qs_list.append(q)
		writer_pool.apply_async(file_writer, args=(q, TF_names_sub, files, args)) 	#, callback = lambda x: finished.append(x) print("Writing time: {0}".format(x)))
		for TF in TF_names_sub:
			qs[TF] = q
	writer_pool.close() #no more jobs applied to writer_pool
	args.qs = qs #qs is a dict

	#Setup scanners pool
	input_arguments = [(chunk, args, motif_list) for chunk in region_chunks]
	task_list = [worker_pool.apply_async(motif_scanning, (chunk, args, motif_list, )) for chunk in region_chunks]
	monitor_progress(task_list, logger)
	results = [task.get() for task in task_list]
	logger.info("Done!")

	#Wait for files to write
	for TF in qs:
		qs[TF].put((None, None))
	writer_pool.join()	

	#Process each file output and write out
	logger.comment("")
	logger.critical("Processing results for each TF")
	task_list = [worker_pool.apply_async(process_TFBS, (file, args)) for file in temp_files]
	worker_pool.close()
	monitor_progress(task_list, logger)
	worker_pool.terminate()
	results = [task.get() for task in task_list]
	logger.info("Done!")

	logger.debug("Joining multiprocessing pools")
	worker_pool.join()	
	writer_pool.join()

	end_time = datetime.now()
	logger.comment("")
	logger.critical("Finished TFBScan run (total time of {0})".format(end_time - begin_time))


#----------------------------------------------------------------------------------------------------------#
if __name__ == "__main__":

	parser = argparse.ArgumentParser()
	parser = add_tfbscan_arguments(parser)
	args = parser.parse_args()

	if len(sys.argv[1:]) == 0:
		parser.print_help()
		sys.exit()

	run_tfbscan(args)

