#!/usr/bin/env python

"""
Classes for working with motifs and scanning with moods

@author: Mette Bentsen
@contact: mette.bentsen (at) mpi-bn.mpg.de
@license: MIT

"""

import numpy as np
import copy
import re
import os
import sys
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
import scipy.spatial.distance as ssd

import base64
import io

#Bio-specific packages
from Bio import motifs
import logomaker
import MOODS.scan
import MOODS.tools
import MOODS.parsers

#Internal
from tobias.utils.regions import OneRegion, RegionList
from tobias.utils.utilities import filafy, num 	#filafy for filenames

"""
def biomotif_to_gimmemotif(biomotif):

	motif_rows = list()

	for pos_id in range(bio_motif.length-1):
		row = list() # each row represents one motif index ( A C G T )
		for letter in range(4):
			row.append(bio_motif.counts[letter][pos_id])
			motif_rows.append(row)

		gimme_motif = Motif(motif_rows) 	# generate gimmemotif motif instance
		
		# Add motif name
		if format == "minimal":
			gimme_motif.id = name_list[i]
		else:
			gimme_motif.id = bio_motif.name
		gimme_motif_list.append(gimme_motif)
"""


def float_to_int(afloat):
	""" Converts integer floats (e.g. 1.0000) to integers """

	elements = str(afloat).split(".")
	
	if len(elements) == 1:	#already int, do nothing
		return(afloat)
	elif len(elements) == 2:

		if float(elements[1]) == 0:	#float is int
			return(int(afloat))
		else:
			return(afloat)
	else:
		pass #if afloat is a string with multiple "."'s

#----------------------------------------------------------------------------------------#
#List of OneMotif objects
class MotifList(list):

	# initialize vars
	#Set by setup moods scanner
	forward = {'names': [], 
			'matrices': [], # pssms
			'thresholds': []}
	reverse = {'names': [], 
			'matrices': [] , # pssms
			'thresholds': []}

	#Scanner
	moods_scanner_forward = None
	moods_scanner_reverse = None
 
	def __init__(self, lst=[]):
		""" Initialize the MotifList object with a list of OneMotif objects (lst) or empty """

		super(MotifList, self).__init__(iter(lst))

		#Update global background when joining motifs with different backgrounds
		self.set_background()

	def __str__(self):
		return("\n".join([str(onemotif) for onemotif in self]))

	def from_file(self, path):
		"""
		Read a file of motifs to MotifList format
		"""
		
		biopython_formats = ["jaspar"]

		#Establish format of motif
		content = open(path).read()
		file_format = get_motif_format(content)

		#For biopython reading
		if file_format == "pfm":
			file_format = "jaspar"

		#parse MEME file
		if file_format == "meme":
			bg_flag = False 	# read background letter frequencies
			proba_flag = False  # read letter-probabilities
			new_motif = True 	# initialize vars for new motif

			#Intialize header variables
			bases = None
			strands = None
			bg = None
			
			#Read meme input line by line
			lines = content.split("\n") + [""] #make sure file ends with empty line
			for line in lines:

				### Header content ###
				# parse alphabet
				if line.startswith("ALPHABET="): # TODO implement for custom alphabet
					bases = list(line.replace("ALPHABET=", "").lstrip().rstrip())

				# parse strands
				elif line.startswith("strands"):
					strands = line.replace("strands: ", "")	#strand string from header

				# find background freq
				elif line.startswith("Background letter frequencies"):
					bg_flag = True # next line is background frequencies

				# parse background freq
				elif bg_flag:
					# Assumes one background line. Might not be the case for custom alphabets.
					bg_flag = False

					bg_and_freq = re.split(r"(?<=\d)\s+", line.strip()) # split after every number followed by a whitespace
					bg_dict = {key: float(value) for key, value in [el.split(" ") for el in bg_and_freq]}
					if bases is None:
						bases = sorted(bg_dict.keys())
					
					bg = np.array([bg_dict[base] for base in bases])

				### Motif content ###
				# parse id, name
				elif line.startswith("MOTIF"):

					#If 'MOTIF' immediately follows previous motif, write out collected probability_matrix to self[-1]
					if proba_flag == True:

						# transpose and convert probability matrix to count using saved .n
						count_matrix = np.array(probability_matrix).T * self[-1].n
						count_matrix = np.round(count_matrix).astype(int) 	#counts are counted to integers
						count_matrix = count_matrix.tolist()

						#Set counts for current OneMotif object
						self[-1].set_counts(count_matrix)	#this also checks for format

					#Initialize new motif
					probability_matrix = []
					self.append(OneMotif(motifid=""))	#initialize dummy OneMotif for filling in

					columns = line.split()
					if len(columns) > 2: # MOTIF, ID, NAME
						motif_id, name = columns[1], columns[2]
					elif len(columns) == 2: # MOTIF, ID
						motif_id, name = columns[1], ""	# name not given
					
					self[-1].id = motif_id
					self[-1].name = name

					#Set any information collected from header (overwrites defaults from OneMotif)
					if bases is not None:
						self[-1].bases = bases
					if strands is not None:
						self[-1].strands = strands
					if bg is not None:
						self[-1].bg = bg

				# find and parse letter probability matrix header
				elif line.startswith("letter-probability matrix"):
					proba_flag = True

					# parse info
					key_value_string = re.sub("letter-probability matrix:\s*", "", line.rstrip())

					# only proceed if there is information
					if len(key_value_string) > 0:
						key_value_split = re.split(r"(?<!=)\s+", key_value_string)	#Split on any space not preceeded by =
						key_value_lists = [re.split(r"=\s*", pair) for pair in key_value_split]
						key_value_dict = {pair[0]: pair[1] for pair in key_value_lists}
						info = key_value_dict
						self[-1].info = info #Add meme-read information to info-dict

						if "nsites" in info:
							self[-1].n = int(float(info["nsites"])) #overwrites OneMotif default

				# parse probability matrix or save motif
				elif proba_flag:

					if re.match(r"^\s*(?![-\s]+)([\d\-\.\se\-]+?)$", line):
						columns = list(map(float, line.split()))

						# check for correct number of columns
						if not len(columns) == len(self[-1].bases):
							sys.exit("Error when reading probability matrix from {0}! Expected {1} columns found {2}!".format(path, len(self[-1].bases), len(columns)))

						# append a single row split into its columns
						probability_matrix.append(columns)

					# motif ended; save this motif
					else:
						proba_flag = False

						# transpose and convert probability matrix to count using saved .n
						count_matrix = np.array(probability_matrix).T * self[-1].n
						count_matrix = np.round(count_matrix).astype(int) 	#counts are counted to integers
						count_matrix = count_matrix.tolist()

						#Set counts for current OneMotif object
						self[-1].set_counts(count_matrix)	#this also checks for format


		# parse PFM/ JASPAR
		elif file_format in biopython_formats:
			with open(path) as f:
				for m in motifs.parse(f, file_format):
					self.append(OneMotif(motifid=m.matrix_id, name=m.name, counts=[m.counts[base] for base in ["A", "C", "G", "T"]]))
					self[-1].biomotifs_obj = m		#biopython motif object	

		else:
			sys.exit("Error when reading motifs from {0}! File format: {1}".format(path, file_format))

		#Check that 'w' fits the length of the read motif
		for motif in self:
			if "w" in motif.info:
				w = int(motif.info["w"])
				l = len(motif.counts[0])
				if w != l:
					sys.exit("Error reading motif '{3}' from {0}! 'w' given in 'letter-probability matrix'-line ({1}) does not match length of the motif read ({2})".format(path, w, l, motif.id))

		#Check if any motifs have length 0 (length == None)
		for motif in self:
			if motif.length == None:
				sys.exit("ERROR: No matrix could be read for motif '{0} {1}' - please check the format of the input motif file.".format(motif.id, motif.name))

		#Final changes to motifs
		for motif in self:

			#Convert float to ints
			for r in range(4):
				for c in range(motif.length):
					motif.counts[r][c] = float_to_int(motif.counts[r][c])

		return(self)

	def to_file(self, path, fmt="pfm"):
		"""
		Write MotifList to motif file

		Parameter:
		----------
		path : string
			Output path
		fmt : string
			Format of motif file (pfm/jaspar/meme)
		"""

		out_string = self.as_string(fmt)

		#Write to output file
		f = open(path, "w")
		f.write(out_string)
		f.close()

		return(self)

	def as_string(self, output_format="pfm"):
		"""
		Returns the MotifList as a string in the given output_format

		Parameter:
		-----------
		output_format : string
			Motif format ("pfm", "jaspar", "meme")
		"""

		out_string = ""

		header = True
		for motif in self:
			# Add a single header for the first motif as it is required in meme format. (For other formats ignored)
			out_string += motif.as_string(output_format=output_format, header=header)
			header = False

		return(out_string)
	
	def get_background(self):
		"""
		Combines background of all motifs to a global background.
  
		Returns:
			numpy array of frequencies or None if no motifs
		"""

		if len(self) > 0:
			global_bg = np.array([0] * len(self[0].bg), dtype=float)
		else:
			return None

		total_n = 0
		for motif in self:
			global_bg += motif.bg * motif.n
			total_n += motif.n

		return global_bg / total_n

	def set_background(self):
		"""
		Set the background using get_background. Sets self.bg to default if there are no motifs in list
		"""

		self.bg = self.get_background() if len(self) > 0 else np.array([0.25] * 4) # (A,C,G,T). Default is equal background if not overwritten by MEME (See OneMotif)

	#---------------- Functions for moods scanning ------------------------#

	def setup_moods_scanner(self, strand="."):
		""" 
		Sets self.moods_scanner_forward and self.moods_scanner_reverse objects 

		Parameters:
			strand (string): Initialize scanner for given strand. '+', '-' or '.' for both.
		"""

		# make sure everything necessary exists
		for motif in self:
			if len(motif.prefix) <= 0 or motif.threshold == None:
				raise Exception("Missing prefix and/or threshold! Please consider running 'motif.set_prefix()' and 'motif.get_threshold()' on the respective motif(s).")
		
		# forward scanner
		if strand in ["+", "."]:
			self.moods_scanner_forward, self.forward_parameters = self.__init_scanner(strand="+")

		# reverse scanner
		if strand in ["-", "."]:
			self.moods_scanner_reverse, self.reverse_parameters = self.__init_scanner(strand="-")

	def __init_scanner(self, strand="+"):
		"""
		Create a new scanner.

		Parameter:
			strand (string): Strand on which the scanner should search. Either "+" or "-".
		Returns:
			Tuple of scanner and parameter dict.
		"""
		if strand == "+":
			motifs = self
		elif strand == "-":
			motifs = self.get_reverse()
			motifs.set_background()	#updates background in case of unequal G/C

		#Calculate pssm if needed
		for motif in motifs:
			if motif.pssm is None:
				motif.get_pssm()

		#Set lists of names/thresholds used for scanning	
		parameters = {'names': [], 'matrices': [], 'thresholds': []}
		for motif in motifs:
			parameters["names"].append(motif.prefix)
			parameters["matrices"].append(motif.pssm)
			parameters["thresholds"].append(motif.threshold)

		#Setup scanner
		scanner = MOODS.scan.Scanner(7)
		scanner.set_motifs(parameters["matrices"], motifs.bg, parameters["thresholds"])

		return (scanner, parameters)

	def scan_sequence(self, seq, region, strand="."):
		"""
		Scan sequence with the motifs in self 
		
		Parameters:
		-----------
		seq : string
			DNA sequence
		region : OneRegion
			OneRegion object fitting seq
		strand : string
			One of "+", "-", "." to indiciate one which strand to search.
		"""
		sites = RegionList()	#Empty regionlist

		#Check that scanners have been initialized
		if strand in ["+", "."] and self.moods_scanner_forward == None:
			self.setup_moods_scanner("+")
			
		if strand in ["-", "."] and self.moods_scanner_reverse == None:
			self.setup_moods_scanner("-")

		#Scan sequence on either + or - strand (or both)
		if strand in ["+", "."]:
			sites += self.__stranded_scan(seq=seq, region=region, strand="+")

		if strand in ["-", "."]:
			sites += self.__stranded_scan(seq=seq, region=region, strand="-")

		return(sites)

	def __stranded_scan(self, seq, region, strand="+"):
		"""
		Scan the sequence based on the given strand. Assumes that the sequence is 5'-3' on the + strand.

		Parameter:
			seq (string): DNA sequence
			region (OneRegion): Object fitting the sequence.
			strand (string): Either "+" or "-".
		Returns:
			List of matches as RegionList object.
		"""

		sites = RegionList()

		if strand == "+":
			scanner = self.moods_scanner_forward
			parameters = self.forward_parameters
		elif strand == "-":
			scanner = self.moods_scanner_reverse
			parameters = self.reverse_parameters

		# scan sequence
		results = scanner.scan(seq)

		# combine results and parameters to RegionList
		for (matrix, name, result) in zip(parameters["matrices"], parameters["names"], results):
			motif_length = len(matrix[0])

			for match in result:
				start = region.start + match.pos # match pos is 1 based
				end = start + motif_length
				score = round(match.score, 5)

				site = OneRegion([region.chrom, start, end, name, score, strand])
				sites.append(site)

		return sites

	#---------------- Functions for motif clustering ----------------------#
	def cluster(self, threshold=0.5, metric = "pcc", clust_method = "average"):
		""" 
		Returns:
		----------
		dict
			A dictionary with keys=cluster names and values=MotifList objects
		"""

		#Needs gimmemotif
		from gimmemotifs.motif import Motif
		from gimmemotifs.comparison import MotifComparer
		sns.set_style("ticks")	#set style back to ticks, as this is set globally during gimmemotifs import

		#Fill in self.gimme_obj variable
		motif_list = [motif.get_gimmemotif().gimme_obj for motif in self]	#list of gimmemotif objects

		#Similarities between all motifs
		mc = MotifComparer()
		score_dict = mc.get_all_scores(motif_list, motif_list, match = "total", metric = metric, combine = "mean")   #metric can be: seqcor, pcc, ed, distance, wic, chisq, akl or ssd
		self.similarity_matrix = generate_similarity_matrix(score_dict)

		# Clustering
		vector = ssd.squareform(self.similarity_matrix.to_numpy())
		self.linkage_mat = linkage(vector, method=clust_method)

		# Flatten clusters
		fclust_labels = fcluster(self.linkage_mat, threshold, criterion="distance")			#cluster membership per motif
		formatted_labels = ["Cluster_{0}".format(label) for label in fclust_labels]

		# Extract motifs belonging to each cluster
		cluster_dict = {label: MotifList() for label in formatted_labels}	#initialize dictionary
		for i, cluster_label in enumerate(formatted_labels):
			cluster_dict[cluster_label].append(self[i])

		return cluster_dict

	def create_consensus(self):
		""" Create consensus motif from MotifList """

		from gimmemotifs.comparison import MotifComparer

		self = [motif.get_gimmemotif() if motif.gimme_obj is None else motif for motif in self]	#fill in gimme_obj if it is not found
		motif_list = [motif.gimme_obj for motif in self]	#list of gimmemotif objects

		if len(motif_list) > 1:
			consensus_found = False
			mc = MotifComparer()

			#Initialize score_dict
			score_dict = mc.get_all_scores(motif_list, motif_list, match = "total", metric = "pcc", combine = "mean")

			while not consensus_found:

				#Which motifs to merge?
				best_similarity_motifs = sorted(find_best_pair(motif_list, score_dict))   #indices of most similar motifs in cluster_motifs

				#Merge
				new_motif = merge_motifs(motif_list[best_similarity_motifs[0]], motif_list[best_similarity_motifs[1]]) 

				del(motif_list[best_similarity_motifs[1]])
				motif_list[best_similarity_motifs[0]] = new_motif

				if len(motif_list) == 1:    #done merging
					consensus_found = True

				else:   #Update score_dict

					#add the comparison of the new motif to the score_dict
					score_dict[new_motif.id] = score_dict.get(new_motif.id, {})

					for m in motif_list:
						score_dict[new_motif.id][m.id] = mc.compare_motifs(new_motif, m, metric= "pcc")
						score_dict[m.id][new_motif.id] = mc.compare_motifs(m, new_motif, metric = "pcc")
	
		#Round pwm values
		gimmemotif_consensus = motif_list[0]
		gimmemotif_consensus.pwm = [[round(f, 5) for f in l] for l in gimmemotif_consensus.pwm]

		#Convert back to OneMotif obj
		onemotif_consensus = gimmemotif_to_onemotif(gimmemotif_consensus)
		onemotif_consensus.gimme_obj = gimmemotif_consensus	

		#Control the naming of the new motif
		all_names = [motif.name for motif in self]
		onemotif_consensus.name = ",".join(all_names[:3])
		onemotif_consensus.name += "(...)" if len(all_names) > 3 else ""

		return(onemotif_consensus)

	def plot_motifs(self, nrow=None, ncol=None, output="motif_plot.png", figsize=None, formation = "row"):
		""" Plot list of motifs to one figure """

		n_motifs = len(self)

		# check formation or set default value
		formation, nrow, ncol = get_formation(formation, ncol, nrow, n_motifs)

		# Check if motifs fit into grid
		if nrow * ncol < n_motifs:
			sys.exit("ERROR: Insufficient space in grid. Please add more rows or columns. Number of motifs: "
					+ str(n_motifs) 
					+ ", Number of spaces: " 
					+ str(ncol*nrow))

		# Get longest motif
		longest_motif = max([len(i[0]) for i in [motif.counts for motif in self]])

		if figsize is None:
			figsize=(longest_motif*0.55*ncol, nrow*3)

		fig = plt.subplots(squeeze=False, figsize=figsize)

		for x, motif in enumerate(self):

			# create axes object for specified position
			ax = plt.subplot2grid((nrow, ncol), formation[x])
			#plot logo to axes object
			motif.create_logo(ax, longest_motif)

		plt.savefig(output)

		return fig

	def make_unique(self):
		""" Make motif ids unique for MotifList """

		seen = {}

		for motif in self:
			m_id = motif.id
			if m_id not in seen:
				seen[m_id] = 1
			else:
				new_id = motif.id + "_" + str(seen[m_id])
				motif.id = new_id
				seen[m_id] += 1

		return(self)

	def get_reverse(self):
		""" Returns a motiflist of reversed motifs."""

		return MotifList([motif.get_reverse() for motif in self])

	def __add__(self, other):
		"""
		Enable + operator to return MotifList.
		"""
		return MotifList(list(self) + list(other))

#--------------------------------------------------------------------------------------------------------#
def gimmemotif_to_onemotif(gimmemotif_obj):
	""" Convert gimmemotif object to OneMotif object """

	# get count matrix
	counts = np.array(gimmemotif_obj.pfm).T.tolist()

	onemotif_obj = OneMotif(motifid=gimmemotif_obj.id, counts=counts)

	return(onemotif_obj)

#--------------------------------------------------------------------------------------------------------#
def generate_similarity_matrix(score_dict):
	"""Generate a similarity matrix from the output of get_all_scores()

	Parameter:
	----------
	score_dict : dict
		a dictionary of dictionarys containing a list of similarity scores

	Returns:
	--------
	DataFrame
		a DataFrame (Pandas) with motif 1 a columns and motif 2 as rows
	"""

	m1_keys = list(score_dict.keys())
	m2_keys = list(score_dict.values())[0].keys()   #should be similar to m1_keys

	m1_labels = [s.replace('\t', ' ') for s in m1_keys] # replacing tabs with whitespace
	m2_labels = [s.replace('\t', ' ') for s in m2_keys]
	
	#Make sure similarity dict is symmetrical:
	similarity_dict = {m:{} for m in m1_labels}  #initialize dict
	for i, m1 in enumerate(m1_keys):
		for j, m2 in enumerate(m2_keys):    
			score = round(1 - np.mean([score_dict[m1][m2][0], score_dict[m2][m1][0]]), 3)
			score = min(score, 1) #fixes that score can go above 1 if the similarity was negative
			
			similarity_dict[m1_labels[i]][m2_labels[j]] = score
			similarity_dict[m2_labels[j]][m1_labels[i]] = score

	#Format similarity dict to dataframe
	similarity_dict_format = {m1: [similarity_dict[m1][m2] for m2 in m2_labels] for m1 in m1_labels}
	dataframe = pd.DataFrame(similarity_dict_format, index = m2_labels).replace(-0, 0)

	return dataframe

#--------------------------------------------------------------------------------------------------------#
def merge_motifs(motif_1, motif_2):
	"""Creates the consensus motif from two provided motifs, using the pos and orientation calculated by gimmemotifs get_all_scores()

	Parameter:
	----------
	motif_1 : Object of class Motif
		First gimmemotif object to create the consensus.
	motif_2 : Object of class Motif
		Second gimmemotif object to create consensus.
	Returns:
	--------
	consensus : Object of class Motif
		Consensus of both motifs with id composed of ids of motifs it was created.
	"""
	from gimmemotifs.comparison import MotifComparer

	mc = MotifComparer()
	_, pos, orientation = mc.compare_motifs(motif_1, motif_2, metric= "pcc")
	consensus = motif_1.average_motifs(motif_2, pos = pos, orientation = orientation)
	consensus.id = motif_1.id + "+" + motif_2.id

	return consensus

#--------------------------------------------------------------------------------------------------------#
def find_best_pair(cluster_motifs, score_dict):
	"""Finds the best pair of motifs based on the best similarity between them im comparison to other motifs in the list.
	Parameter:
	----------
	clusters_motifs : list
		List of motifs assigned to the current cluster.
	score_dict : dict
		Dictionary conatining list of [similarity_score, pos, strand] as values and motif names as keys.
	Returns:
	--------
	best_similarity_motifs : list of two elements
		List of the best pair of motifs found based on the similarity.
	"""

	best_similarity = 0
	for i, m in enumerate(cluster_motifs):
		for j, n in enumerate(cluster_motifs):
			if m.id is not n.id: 
				this_similarity = score_dict[m.id][n.id][0]
				if this_similarity > best_similarity:
					best_similarity = this_similarity
					best_similarity_motifs = [i, j] #index of the most similar motifs in cluster_motifs

	return best_similarity_motifs

#--------------------------------------------------------------------------------------------------------#
def get_formation(formation, ncol, nrow, nmotifs):
	""" check formation or set formation to one of the existing options """

	# if ncol and/or nrow is missing automatically set fitting parameters 
	if formation != "alltoone":
		if ncol is None and nrow is None:
			half_nmotifs = math.ceil(math.sqrt(nmotifs))
			ncol, nrow = half_nmotifs, half_nmotifs
		else:
			if ncol is None:
				ncol = math.ceil(nmotifs/nrow)
			if nrow is None:
				nrow = math.ceil(nmotifs/ncol)

	if isinstance(formation, str):
		
		if formation == "row":

			# fill plot left to right

			formation = list()
			rows = list(range(nrow))
			for row in rows:
				for col in range(ncol):
					formation.append((row,col))

		elif formation == "col":

			# fill plot top to bottom

			formation = list()
			rows = list(range(nrow))
			for col in range(ncol):
				for row in rows:
					formation.append((row,col))

		elif formation == "alltoone":

			# fill first column execpt for one motif
			# ignores parameter ncol and nrow

			formation = list()
			rows = list(range(nmotifs-1))
			for row in rows:
				formation.append((row,0))
			formation.append((math.ceil(len(rows)/2)-1, 1))

			ncol = 2
			nrow = len(rows)

		else:
			sys.exit("ERROR: Unknown formation setting.")
	else:

		# Check if formation fits to grid
		formation_max_row = max([i[0] for i in formation])
		formation_max_col = max([i[1] for i in formation])
		if nrow < formation_max_row or ncol < formation_max_col:
			sys.exit("ERROR: Grid is to small for specified formation")

	return formation, nrow, ncol


#----------------------------------------------------------------------------------------#
#Contains info on one motif formatted for use in moods
class OneMotif:

	id = "" # motif id (should be unique)
	name = "" # motif name (does not have to be unique)
	bases = ["A", "C", "G", "T"] # alphabet order must correspont with counts matrix!
	bg = np.array([0.25,0.25,0.25,0.25]) # background set to equal by default
	strands = "+ -"		# default meme strands
	n = 20				# default number of sites used for creating motif
	length = None 		# length of the motif
 
	threshold = None # moods hit threshold calculated in get_threshold()
	gimme_obj = None # gimmemotif obj
	info = {} # additional motif information used in meme format
	prefix = "" # output prefix set in set_prefix()
 
	# matrices
	counts = None # base counts per position; list of 4 lists (A,C,G,T) (each as long as motif)
	pfm = None # position frequency matrix, i.e. counts / sum of counts per position
	pssm = None # The log-odds scoring matrix (pssm) calculated from get_pssm.

	def __init__(self, motifid, counts=None, name=None):

		self.id = motifid if motifid != None else ""		#should be unique
		self.name = name if name != None else "" 			#does not have to be unique

		# sets counts, length and n
		if not counts is None:
			self.set_counts(counts)

	def __str__(self):
		""" Format used for printing """
		return("{0}".format(self.__dict__))

	def set_prefix(self, naming="name_id"):
		""" Set name to be used in 4th column and as output prefix """

		if naming == "name":
			prefix = self.name
		elif naming == "id":
			prefix = self.id
		elif naming == "name_id":
			prefix = self.name + "_" + self.id
		elif naming	== "id_name":
			prefix = self.id + "_" + self.name
		else:
			prefix = "None"

		self.prefix = filafy(prefix)
		return(self)

	def get_pfm(self):
		""" Set the position frequency matrix (self.pfm) from self.count. The frecquency matrix is set as the relative frequency.
			This matrix is also sometimes called the PPM or probability matrix.
		"""

		self.pfm = self.counts / np.sum(self.counts, axis=0)
		return(self)

	def get_gimmemotif(self):
		""" Get gimmemotif object for motif 
			Reads counts from self.counts 
		"""

		from gimmemotifs.motif import Motif

		self.length = len(self.counts[0])

		motif_rows = []
		for pos_id in range(self.length):
			row = [self.counts[letter][pos_id] for letter in range(len(self.bases))] 	# each row represents one position in motif ( A C G T )
			motif_rows.append(row)

		# generate gimmemotif motif instance
		self.gimme_obj = Motif()
  
		# populate empty object
		self.gimme_obj.id = self.id + " " + self.name
		self.gimme_obj.pfm = motif_rows
		self.gimme_obj.pwm = self.gimme_obj.pfm_to_pwm(motif_rows)

		return(self)
		
	def get_biomotif(self):
		""" Get biomotif object for motif """

		# TODO implement
		self.biomotif_obj = ""

	def get_reverse(self):
		""" Reverse complement motif """

		# reverse counts
		rev_counts = [[],[],[],[]]
		rev_counts[0] = self.counts[3][::-1] # rev T => A
		rev_counts[1] = self.counts[2][::-1] # rev G => C
		rev_counts[2] = self.counts[1][::-1] # rev C => G
		rev_counts[3] = self.counts[0][::-1] # rev A => T
		rev_bg = self.bg[[3, 2, 1, 0]]	# reverse background; will be the same if A/T, C/G ratios are balanced

		# Create reverse motif obj and fill in 
		reverse_motif = OneMotif(motifid=self.id, counts=rev_counts, name=self.name)
		reverse_motif.info = self.info 	# add info from original motif
		reverse_motif.bg = rev_bg		# add background
		reverse_motif.prefix = self.prefix
		reverse_motif.threshold = self.threshold
		# TODO reverse strand; only applicable for meme files

		return(reverse_motif)	#OneMotif object

	def get_pssm(self, ps=0.01):
		""" Calculate pssm from pfm. ps is a pseudocount. Note: PSSM is sometimes called PWM in literature"""

		if self.pfm is None:
			self.get_pfm()

		bg_col = self.bg.reshape((-1,1))
		pseudo_vector = ps * bg_col

		pfm_pseudocount = np.true_divide(self.pfm + pseudo_vector, np.sum(self.pfm + pseudo_vector, axis=0)) #pfm with added pseudocounts
		self.pssm = np.log(pfm_pseudocount) - np.log(bg_col) #pfm/bg

		return(self)

	def get_threshold(self, pvalue=0.0001):
		""" Get threshold for moods scanning """

		if self.pssm is None:
			self.get_pssm()

		pssm_tuple = tuple([tuple(row) for row in self.pssm])
		self.threshold = MOODS.tools.threshold_from_p(pssm_tuple, self.bg, pvalue, 4)

		return(self)

	def information_content(self, ps=0.01):
		'''
		Calculate the information content of the given motif.
		The information content is calculated in bits and is saved in the .ic-variable
		
		Uses:
		:ps: a pseudocount to add to the pfm
		:self.pfm: Position probability matrix as dataframe where rows are positions and columns bases.
		:self.bg: Background frequencies
		'''

		if self.pfm is None:
			self.get_pfm()
	
		#Calculate information content
		bg_col = self.bg.reshape((-1,1))
		self.information = self.pfm * (np.log2(self.pfm + ps) - np.log2(bg_col + ps))
		self.ic = np.sum(self.information)	#sum of all information content
		
		self.bits = self.pfm * np.sum(self.information, axis=0)	#bits is the pfm-fraction of information content per column

		return(self)

	def gc_content(self):
		'''
		Calculate the GC content of the motif.
		The GC content is calculated in percent (0-1) and is stored in the .gc-variable.

		Uses:
		:self.pfm: Position probability matrix as dataframe where rows are positions and columns bases.
		'''

		if self.pfm is None:
			self.get_pfm()

		# Calculate GC content per position
		self.gc_positions = self.pfm[self.bases.index('G')] + self.pfm[self.bases.index('C')]

		self.gc = np.mean(self.gc_positions)

		return(self)

	def logo_to_file(self, filename, ylim=(0, 2)):
		"""
		Plots the motif to pdf/png/jpg file

		Parameters:
			filename (string): Name of the output file.
			ylim (tuple, list or string): Forwarded to create_logo().
		"""

		ext = os.path.splitext(filename)[-1]

		#Currently only working with pdf
		#filename = filename.replace(ext, ".pdf")	#hack

		if ext == "jpg" :
			filename[-3:] = "png"
			warnings.warn("The 'jpg' format is not supported for motif image. Type is set to 'png'")

		#self.gimme_obj.to_img(filename)	
		logo = self.create_logo(ylim=ylim)	#returns a logo object
		logo.fig.savefig(filename)
		plt.close(logo.fig)

	def get_base(self):
		""" Get base64 string for plotting in HTML """

		image = io.BytesIO()

		logo = self.create_logo()
		logo.fig.savefig(image)
		self.base = base64.encodestring(image.getvalue()).decode("utf-8") 
		self.base = self.base.replace("\n", "")	#replace automatic \n

		return(self)

	def create_logo(self, ax = None, motif_len=None, ylim=(0, 2)):
		"""
		Creates motif logo in axes object

		Parameters:
			ax (matplitlib.axes): Axes object, where to logo should live. Default = None
			motif_len (int): Number of bases displayed in plot. Default length of motif (all).
			ylim (tuple, list or string): Either 'auto' to autoscale yaxis to max or tuple/ list of parameters (*args) to matplotlib.axes.Axes.set_ylim. Default 0 to 2.

		Returns:
			logomaker.Logo object
		"""

		# convert to pandas dataframe
		df = pd.DataFrame(self.counts).transpose()
		df.columns = self.bases

		if not motif_len:
			motif_len = df.shape[0]

		# transform matrix to information based values
		info_df = logomaker.transform_matrix(df, from_type="counts", to_type="information")
		self.info_df = info_df

		# compute y-axis limits and ticks
		ylim = (0, info_df.sum(axis=1).max()) if ylim == "auto" else ylim
		# validate ylim
		if not isinstance(ylim, (list, tuple)):
			raise ValueError("Parameter ylim should be either 'auto' or a list/ tuple of arguments to matplotlib.axes.Axes.set_ylim.")
		tick_num = 4
		step = (ylim[1] - ylim[0]) / tick_num
		yticks = [ylim[0]+step*i for i in range(tick_num + 1)]

		# create Logo object
		logo = logomaker.Logo(info_df, ax = ax)

		# style
		logo.style_xticks(rotation=0, fmt='%d', anchor=0)
		logo.ax.set_ylim(*ylim)
		logo.ax.set_xlim(-0.5, motif_len-0.5)
		logo.ax.set_yticks(yticks, minor=False)
		logo.ax.set_xticklabels(range(1,motif_len+1)) #0 -> 1-based coordinates
		logo.style_spines(visible=False)
		logo.style_spines(spines=['left', 'bottom'], visible=True)
		logo.ax.xaxis.set_ticks_position('none')
		logo.ax.xaxis.set_tick_params(pad=-1)

		return logo

	def set_counts(self, counts):
		""" 
		Set the count matrix. Also updates number of sites (n) and motif length (length).

		Parameters:
			counts (list): List of length 4. Contains lists of base counts for each position.

		Returns:
			self (OneMotif): OneMotif object with modified counts.
  		"""
    	# check counts input
		if len(counts) != 4:
			raise ValueError("Input counts must be of length 4. Received length {0} for motif (id='{1}', name='{2}'). Input counts are: {3}".format(len(counts), self.id, self.name, counts))
		lengths = [len(base) for base in counts]
		if len(set(lengths)) != 1:
			raise ValueError("All lists in counts must be of same length.")

		# add counts and associated length/n to OneMotif object
		self.counts = counts
		self.length = lengths[0]	#update motif length
		self.n = np.sum([row[0] for row in counts])

		return self

	def as_string(self, output_format="pfm", header=True):
		"""
		Returns the OneMotif object as a string in the given output_format.
  
		Parameters:
			output_format (string): Set the output format one of ["pfm", "jaspar", "meme"].

			header (bool): Whether a header should be generated (meme only).
   
		Returns:
			(string): String in the chosen file format.
		"""
		out_string = ""
  
		if output_format in ["pfm", "jaspar"]:
			# create id line
			out_string += f">{self.id}\t{self.name}\n"

			# add counts
			for index, base in enumerate(self.bases):
				row = " ".join(map(str, self.counts[index]))

				if output_format == "jaspar":
					row = f"{base} [{row} ]"

				out_string += row + "\n"
			
		elif output_format == "meme":

			# create meme header
			if header:
				# TODO read meme version from original file (or default to version 4)
				meme_header = "MEME version 4\n\n"
				meme_header += "ALPHABET= {0}\n\n".format("".join(self.bases))
				meme_header += "strands: {0}\n\n".format(self.strands)
				meme_header += "Background letter frequencies\n"
				meme_header += " ".join([f"{self.bases[i]} {self.bg[i]}" for i in range(4)]) + "\n"

				out_string += meme_header

			# add motif information to string
			out_string += f"\nMOTIF {self.id} {self.name}\n"
			out_string += f"letter-probability matrix: alength= {len(self.bases)} w= {self.length} nsites= {int(round(self.n))} E= 0\n" # omit e_value as it could be out of context
			
			# add pfm
			if self.pfm is None:
				self.get_pfm()

			precision = 6
			for row in self.pfm.T: #add frequency per position
				out_string += " {0}\n".format("  ".join(map(lambda f: format(round(f, precision), f".{precision}f"), row)))

		# TODO also implementation in from_file needed
		#elif output_format == "transfac":
		#	out_string += self.get_gimmemotif().gimme_obj.to_transfac()

		else:
			raise ValueError("Format " + output_format + " is not supported")

		return out_string

	def to_file(self, output_file, fmt="pfm"):
		"""
		Save motif as a file.

		Parameters:
			output_file (string): Name of the output file.

			fmt (string): Format of the output file. One of ["pfm", "jaspar", "meme"].

		Returns:
			self (OneMotif)
		"""

		with open(output_file, "w") as f:
			f.write(self.as_string(output_format=fmt))

		return self

	@staticmethod
	def from_fasta(fasta, motifid, name=None):
		"""
		Create motif from fasta.
		Will use captital letters as motif sites (see JASPAR sites format).

		Parameters:
			fasta (string): Path to fasta file.

			motifid (string): Unique id of the motif.

			name (string): Name of the motif. Defaults to 'None'.

		Returns:
			OneMotif object
		"""
		with open(fasta) as handle:
			motif = motifs.read(handle, "sites")

		return OneMotif(motifid=motifid, counts=[motif.counts[base] for base in ["A", "C", "G", "T"]], name=name)

	def __repr__(self):
		"""
		OneMotif representation.
		"""
		return f"<OneMotif: {self.id}{' ' + self.name if len(self.name) > 0 else ''}>"


###########################################################

def get_motif_format(content):
	""" Get motif format from string of content """
	
	#Estimate input format
	if re.match(r".*MEME version.+", content, re.DOTALL) is not None: # MOTIF\s.+letter-probability matrix.+[\d\.\s]+", content, re.MULTILINE) is not None:
		motif_format = "meme"

	elif re.match(r">.+A.+\[", content, re.DOTALL) is not None:
		motif_format = "jaspar"

	elif re.match(r">.+", content, re.DOTALL) is not None:
		motif_format = "pfm"

	elif re.match(r"AC\s.+", content, re.DOTALL) is not None:
		motif_format = "transfac"
	
	else:
		motif_format = "unknown"

	return(motif_format)
