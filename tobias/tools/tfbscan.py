#!/usr/bin/env python

"""
TFBScan.py scans for positions of transcription factor binding sites across the genome

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

import pysam	#for reading fastafile


from tobias.parsers import add_tfbscan_arguments
from tobias.utils.utilities import *
from tobias.utils.regions import *
from tobias.utils.sequences import get_gc_content
from tobias.utils.motifs import *
from tobias.utils.logger import *

#----------------------------------------------------------------------------------------------------------#
def motif_scanning(regions, args, motifs_obj):
	""" motifs_obj is a MotifList object """

	fasta_obj = pysam.FastaFile(args.fasta)
	qs = args.qs

	#Scan motifs_obj against sequences per region
	all_TFBS = {name: RegionList() for name in [motif.prefix for motif in motifs_obj]}

	for region in regions:

		seq = fasta_obj.fetch(region.chrom, region.start, region.end)
		region_TFBS = motifs_obj.scan_sequence(seq, region)		#Scan sequence, returns list of OneRegion TFBS
		
		#Add region columns if chosen
		if args.add_region_columns:
			for TFBS in region_TFBS:
				TFBS.extend(region[3:])

		#Check format of region chromosome and convert sites if needed
		m = re.match(r"(.+)\:([0-9]+)\-([0-9]+)\s+.+", region.chrom)
		if m:
			reg_chrom, reg_start, reg_end = m.group(1), m.group(2), m.group(3)
			for TFBS in region_TFBS:
				TFBS.chrom = region.chrom
				TFBS.start = int(reg_start) + TFBS.start
				TFBS.end = int(reg_start) + TFBS.end

		#Split to single TFs
		for TFBS in region_TFBS:
			all_TFBS[TFBS.name].append(TFBS)

	#Resolve overlaps
	if args.keep_overlaps == False:	#default
		for name in all_TFBS:
			all_TFBS[name] = all_TFBS[name].resolve_overlaps()

	#Write to queue
	for name in all_TFBS:
		bed_content = all_TFBS[name].as_bed()	#string 
		qs[name].put((name, bed_content))		#send to writers

	return(1)


def process_TFBS(infile, args):
	""" Process TFBS bedfile - remove duplicates and sort """

	if args.outdir != None:
		outfile = infile.replace(".tmp", ".bed")	 #single files, the names are controlled as motif names
	elif args.outfile != None:
		outfile = args.outfile

	#Sort and print unique lines
	os.system("sort -k1,1 -k2,2n {0} | uniq > {1}".format(infile, outfile))

	#Remove the tmp file
	if not args.debug:
		try:
			os.remove(infile)
		except:
			print("Error removing file {0}".format(infile))

	return(1)


#----------------------------------------------------------------------------------------------------------#
def run_tfbscan(args):

	###### Check input arguments ######
	check_required(args, ["motifs", "fasta"])				#Check input arguments
	check_files([args.motifs, args.fasta, args.regions]) 	#Check if files exist

	##Test input
	if args.outdir != None and args.outfile != None:								#Error - both set
		sys.exit("ERROR: Please choose either --outdir or --outfile")
	elif ((args.outdir == None or args.outdir != None) and args.outfile == None): 	#Separate files	
		args.outdir = "tfbscan_output/" if args.outdir == None else args.outdir
		make_directory(args.outdir) #Check and create output directory
	elif args.outdir == None and args.outfile != None: 								#Joined file
		check_files([args.outfile], "w")


	###### Create logger and write argument overview ######
	logger = TobiasLogger("TFBScan", args.verbosity)
	logger.begin()
	parser = add_tfbscan_arguments(argparse.ArgumentParser())
	
	logger.arguments_overview(parser, args)
	
	if args.outfile != None:
		logger.output_files([args.outfile])

	######## Read sequences from file and estimate background gc ########
	
	logger.info("Handling input files")
	logger.info("Reading sequences from fasta")

	fastafile = pysam.FastaFile(args.fasta)
	fasta_chrom_info = dict(zip(fastafile.references, fastafile.lengths))
	fastafile.close()
	logger.stats("- Found {0} sequences in fasta".format(len(fasta_chrom_info)))
	
	#Create regions available in fasta 	
	logger.info("Setting up regions")
	fasta_regions = RegionList([OneRegion([header, 0, fasta_chrom_info[header]]) for header in fasta_chrom_info])

	#If subset, setup regions
	if args.regions:
		regions = RegionList().from_bed(args.regions)

	else:	#set up regions from fasta references
		regions = fasta_regions
		regions = regions.apply_method(OneRegion.split_region, 1000000)	
		regions = regions.apply_method(OneRegion.extend_reg, 50)		#extend to overlap at junctions	

	#Clip regions at chromosome boundaries
	regions = regions.apply_method(OneRegion.check_boundary, fasta_chrom_info, "cut")
	if len(regions) == 0:
		logger.error("No regions found.")
		sys.exit()
	logger.info("- Total of {0} regions (after splitting)".format(len(regions)))
	
	#Background gc
	if args.gc == None:
		logger.info("Estimating GC content from fasta (set --gc to skip this step)")
		args.gc = get_gc_content(regions, args.fasta)
		logger.info("- GC content: {0}".format(round(args.gc, 5)))
	
	bg = np.array([(1-args.gc)/2.0, args.gc/2.0, args.gc/2.0, (1-args.gc)/2.0])

	#Split regions
	region_chunks = regions.chunks(args.split)
	

	#################### Read motifs from file ####################

	logger.info("Reading motifs from file")

	motif_list = MotifList().from_file(args.motifs)
	logger.stats("- Found {0} motifs".format(len(motif_list)))
	
	logger.debug("Getting motifs ready")
	motif_list.bg = bg
	for motif in motif_list:
		motif.set_prefix(args.naming)
		motif.bg = bg
		motif.get_pssm()
	
	motif_names = list(set([motif.prefix for motif in motif_list]))

	#Calculate scanning-threshold for each motif
	pool = mp.Pool(processes=args.cores)
	outlist = pool.starmap(OneMotif.get_threshold, itertools.product(motif_list, [args.pvalue])) 
	motif_list = MotifList(outlist)	

	pool.close()
	pool.join()


	#################### Find TFBS in regions #####################

	logger.comment("")
	logger.info("Scanning for TFBS with all motifs")

	manager = mp.Manager()

	if args.outdir != None:
		writer_cores = max(1,int(args.cores*0.1))
		worker_cores = max(1,args.cores - writer_cores)

	elif args.outfile != None:	#Write to one file
		writer_cores = 1
		worker_cores = max(1,args.cores - writer_cores)

	#Setup pools
	logger.debug("Writer cores: {0}".format(writer_cores))
	logger.debug("Worker cores: {0}".format(worker_cores))
	worker_pool = mp.Pool(processes=worker_cores, maxtasksperchild=1)
	writer_pool = mp.Pool(processes=writer_cores)

	#Setup bed-writers based on --outdir or --outfile
	temp_files = []
	qs = {}
	TF_names_chunks = [motif_names[i::writer_cores] for i in range(writer_cores)]
	for TF_names_sub in TF_names_chunks:

		#Skip over any empty chunks
		if len(TF_names_sub) == 0:
			continue

		logger.debug("Creating writer queue for {0}".format(TF_names_sub))

		if args.outdir != None:
			files = [os.path.join(args.outdir, TF + ".tmp") for TF in TF_names_sub]
			temp_files.extend(files)
		elif args.outfile != None:
			files = [args.outfile + ".tmp" for TF in TF_names_sub]		#write to the same file for all
			temp_files.append(files[0])
			
		q = manager.Queue()

		TF2files = dict(zip(TF_names_sub, files))
		logger.debug("TF2files dict: {0}".format(TF2files))
		writer_pool.apply_async(file_writer, args=(q, TF2files, args)) 	#, callback = lambda x: finished.append(x) print("Writing time: {0}".format(x)))
		for TF in TF_names_sub:
			qs[TF] = q
	writer_pool.close() #no more jobs applied to writer_pool
	args.qs = qs 		#qs is a dict

	#Setup scanners pool
	input_arguments = [(chunk, args, motif_list) for chunk in region_chunks]
	task_list = [worker_pool.apply_async(motif_scanning, (chunk, args, motif_list, )) for chunk in region_chunks]
	monitor_progress(task_list, logger)
	results = [task.get() for task in task_list]	#1s

	#Wait for files to write
	for TF in qs:
		qs[TF].put((None, None))

	writer_pool.join()

	#Process each file output and write out
	logger.comment("")
	logger.info("Processing results from scanning")
	logger.debug("Running processing for files: {0}".format(temp_files))
	task_list = [worker_pool.apply_async(process_TFBS, (file, args)) for file in temp_files]
	worker_pool.close()
	monitor_progress(task_list, logger)
	worker_pool.terminate()
	results = [task.get() for task in task_list]

	logger.debug("Joining multiprocessing pools")
	worker_pool.join()	
	writer_pool.join()

	logger.end()

#----------------------------------------------------------------------------------------------------------#
if __name__ == "__main__":

	parser = argparse.ArgumentParser()
	parser = add_tfbscan_arguments(parser)
	args = parser.parse_args()

	if len(sys.argv[1:]) == 0:
		parser.print_help()
		sys.exit()

	run_tfbscan(args)

