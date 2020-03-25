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
#from matplotlib.text import TextPath
#from matplotlib.patches import PathPatch
#from matplotlib.font_manager import FontProperties
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

	def __init__(self, lst=[]):
		""" Initialize the MotifList object with a list of OneMotif objects (lst) or empty """

		super(MotifList, self).__init__(iter(lst))

		self.bg = np.array([0.25,0.25,0.25,0.25])

		#Set by setup moods scanner
		self.names = []
		self.matrices = [] 	#pssms
		self.strands = []
		self.thresholds = []

		#Scanner
		self.moods_scanner = None

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

		#Read motifs
		if file_format == "meme":
			
			lines = content.split("\n")
			for line in lines:
				columns = line.strip().split()

				if line.startswith("MOTIF"):
					self.append(OneMotif()) 				#create new motif
					self[-1].input_format = file_format
					self[-1].counts = [[] for _ in range(4)]
					self[-1].n = 1 #preliminarily 1

					#Get id/name of motif
					if len(columns) > 2: #MOTIF, ID, NAME
						motif_id, name = columns[1], columns[2]
					elif len(columns) == 2: # MOTIF, ID
						motif_id, name = columns[1], ""	#name not given

					self[-1].id = motif_id
					self[-1].name = name
				
				else:
					if len(self) > 0: #if there was already one motif header found
					
						#If line contains counts
						if re.match(r"^[\s]*([\d\.\s]+)$", line):	#starts with any number of spaces (or none) followed by numbers
							for i, col in enumerate(columns):
								self[-1].counts[i].append(num(col))
						elif re.match("^letter-probability", line):
							#example: "letter-probability matrix: alength= 4 w= 6 nsites= 24 E= 0"

							m = re.search(r"nsites= ([0-9]+)", line)
							if m is not None:
								self[-1].n = int(m.group(1))

			#Multiply all counts with number of sequences
			for motif in self:
				for i, nuc_counts in enumerate(motif.counts):
					motif.counts[i] = [pos * motif.n for pos in nuc_counts]

		elif file_format in biopython_formats:
			with open(path) as f:
				for m in motifs.parse(f, file_format):
					self.append(OneMotif(motifid=m.matrix_id, name=m.name, counts=[m.counts[base] for base in ["A", "C", "G", "T"]]))
					self[-1].biomotifs_obj = m		#biopython motif object	

		else:
			sys.exit("Error when reading motifs from {0}! File format: {1}".format(path, file_format))

		#Gimmemotifs functionality removed in 0.11.0
		"""
		elif file_format in gimmemotif_formats:
			gimme_motif_list = read_motifs(infile = path, fmt = file_format)
			for gimmemotif in gimme_motif_list:
				onemotif_obj = gimmemotif_to_onemotif(gimmemotif)
				self.append(onemotif_obj)	#add OneMotif object to list
		"""

		#Check correct format of pfms
		for motif in self:
			nuc, length = np.array(motif.counts).shape
			motif.length = length
			if nuc != 4:
				sys.exit("ERROR: Motif {0} has an unexpected format and could not be read".format(motif))

		#Convert float to ints
		for motif in self:
			for r in range(4):
				for c in range(motif.length):
					motif.counts[r][c] = float_to_int(motif.counts[r][c])
	
		#Fill in motifs with additional parameters; Estimate widths and n_sites
		for i, motif in enumerate(self):
			self[i].n = int(round(sum([base_counts[0] for base_counts in motif.counts])))
			self[i].length = len(motif.counts[0])

		return(self)

	def to_file(self, path, fmt="pfm"):
		"""
		Write MotifList to motif file

		Parameter:
		----------
		path : string
			Output path
		fmt : string
			Format of motif file (pfm/jaspar/meme/transfac)
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

		bases = ["A", "C", "G", "T"]
		out_string = ""

		#Establish which output format
		if output_format in ["pfm", "jaspar"]:
			for motif in self:
				out_string += ">{0}\t{1}\n".format(motif.id, motif.name)

				#Convert counts to string
				counts_type = float if float in [type(element) for element in motif.counts[0]] else int	#Estimate from first base counts
				if counts_type == float:
					counts_string = [["{0:.5f}".format(element) for element in row_counts] for row_counts in motif.counts]
				else:
					counts_string = [[str(element) for element in row_counts] for row_counts in motif.counts]
					element_width = max(sum([[len(element) for element in row_counts] for row_counts in counts_string], []))	#Length of longest element to align to
					counts_string = [[" " * (element_width - len(element)) + element for element in row_counts] for row_counts in counts_string]

				#Write out to format
				for i, base_counts in enumerate(counts_string):
					out_string += "{0} [ {1} ]\n".format(bases[i], "  ".join(base_counts)) if output_format == "jaspar" else "  ".join(base_counts) + "\n"
				out_string += "\n"

		elif output_format == "meme":
			
			meme_header = "MEME version 4\n\n"
			meme_header += "ALPHABET=ACGT\n\n"
			meme_header += "strands: + -\n\n"
			meme_header += "Background letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n\n"
			out_string += meme_header

			for motif in self:
				out_string += "MOTIF\t{0}\t{1}\n".format(motif.id, motif.name)
				out_string += "letter-probability matrix: alength=4 w={0} nsites={1} E=0\n".format(motif.length, motif.n)

				for i in range(motif.length):
					row = [float(motif.counts[j][i]) for j in range(4)] 	#row contains original row from content
					n_sites = round(sum(row), 0)
					row_freq = ["{0:.5f}".format(num/n_sites) for num in row] 
					out_string += "  ".join(row_freq) + "\n"
				
				out_string += "\n"

		#elif output_format == "transfac":
		#	for motif in self:
		#		out_string += motif.to_transfac()

		else:
			raise ValueError("Format " + output_format + " is not supported")	

		return(out_string)			


	#---------------- Functions for moods scanning ------------------------#

	def setup_moods_scanner(self):
		""" Sets self.moods_scanner object """

		tups = [(motif.prefix, motif.strand, motif.pssm, motif.threshold) for motif in self] 		#list of tups
		if len(tups) > 0:
			self.names, self.strands, self.matrices, self.thresholds = list(map(list, zip(*tups))) 	#get "columns"
		else:
			self.names, self.strands, self.marices, self.thresholds = ([], [], [], [])

		scanner = MOODS.scan.Scanner(7)
		scanner.set_motifs(self.matrices, self.bg, self.thresholds)

		self.moods_scanner = scanner

	def scan_sequence(self, seq, region):
		"""
		Scan sequence with the motifs in self 
		
		Parameters:
		-----------
		seq : string
			DNA sequence
		region : OneRegion
			OneRegion object fitting seq
		"""

		if self.moods_scanner == None:
			self.setup_moods_scanner()

		#Scan sequence
		results = self.moods_scanner.scan(seq)

		#Convert results to RegionList
		sites = RegionList()	#Empty regionlist
		for (matrix, name, strand, result) in zip(self.matrices, self.names, self.strands, results):
			motif_length = len(matrix[0])
			for match in result:
				start = region.start + match.pos 	#match pos is 1 based
				end = start + motif_length		
				score = round(match.score, 5)

				site = OneRegion([region.chrom, start, end, name, score, strand])	#Create OneRegion obj
				sites.append(site)

		return(sites)

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

		self = [motif.get_gimmemotifs() if motif.gimme_obj is None else motif for motif in self]	#fill in gimme_obj if it is not found
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

#--------------------------------------------------------------------------------------------------------#
def gimmemotif_to_onemotif(gimmemotif_obj):
	""" Convert gimmemotif object to OneMotif object """

	length = len(gimmemotif_obj.pwm)

	onemotif_obj = OneMotif(motifid=gimmemotif_obj.id)
	for pos in range(length):
		for base in range(4):
			onemotif_obj.counts[base].append(gimmemotif_obj.pfm[pos][base])

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

	bases = ["A", "C", "G", "T"]

	def __init__(self, motifid=None, name=None, counts=None):
		
		self.id = motifid if motifid != None else ""		#should be unique
		self.name = name if name != None else "" 			#does not have to be unique

		self.prefix = "" 		#output prefix set in set_prefix
		self.counts = counts if counts != None else [[] for _ in range(4)] 	#counts, list of 4 lists (A,C,G,T) (each as long as motif)
		self.strand = "+"		#default strand is +
		self.length = len(counts[0]) if counts != None else None	#length of motif

		#Set later
		self.pfm = None
		self.bg = np.array([0.25,0.25,0.25,0.25]) 	#background set to equal by default
		self.pssm = None 							#pssm calculated from get_pssm
		self.threshold = None 						#threshold calculated from get_threshold
		self.gimme_obj = None						#gimmemotif obj

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
		""" Set self.pfm from self.counts """

		self.pfm = self.counts / np.sum(self.counts, axis=0)
		return(self)

	def get_gimmemotif(self):
		""" Get gimmemotif object for motif 
			Reads counts from self.counts """

		from gimmemotifs.motif import Motif

		self.length = len(self.counts[0])

		motif_rows = []
		for pos_id in range(self.length):
			row = [self.counts[letter][pos_id] for letter in range(4)] 	# each row represents one position in motif ( A C G T )
			motif_rows.append(row)

		self.gimme_obj = Motif(motif_rows) 	# generate gimmemotif motif instance
		self.gimme_obj.id = self.id + " " + self.name

		return(self)
		
	def get_biomotif(self):
		""" Get biomotif object for motif """

		self.biomotif_obj = ""

	def get_reverse(self):
		""" Reverse complement motif """
		if self.pfm is None:
			self.get_pfm()

		#Create reverse motif obj
		reverse_motif = OneMotif()	#empty
		for att in ["id", "name", "prefix", "strand", "length"]:
			setattr(reverse_motif, att, getattr(self, att))

		reverse_motif.strand = "-" if self.strand == "+" else "+"
		reverse_motif.pfm = MOODS.tools.reverse_complement(self.pfm, 4)
		return(reverse_motif)	#OneMotif object

	def get_pssm(self, ps=0.01):
		""" Calculate pssm from pfm """

		if self.pfm is None:
			self.get_pfm()

		bg_col = self.bg.reshape((-1,1))
		pseudo_vector = ps * bg_col

		pssm = np.log(np.true_divide(self.pfm + pseudo_vector, np.sum(self.pfm + pseudo_vector, axis=0))) - np.log(bg_col)
		pssm = tuple([tuple(row) for row in pssm])
		self.pssm = pssm
		return(self)

	def get_threshold(self, pvalue):
		""" Get threshold for moods scanning """
		if self.pssm is None:
			self.get_pssm()

		self.threshold = MOODS.tools.threshold_from_p(self.pssm, self.bg, pvalue, 4)
		return(self)

	def calc_bit_score(self):
		""" Bits for logo plots (?) """
		if self.pfm is None:
			self.get_pfm()

		pfm_arr = np.copy(self.pfm)
		pfm_arr[pfm_arr == 0] = np.nan

		#Info content per pos
		entro = pfm_arr * np.log2(pfm_arr)
		entro[np.isnan(entro)] = 0
		info_content = 2 - (- np.sum(entro, axis=0))		#information content per position in motif
		self.ic = info_content
		self.bits = self.pfm * info_content	

	def logo_to_file(self, filename):
		""" Plots the motif to pdf/png/jpg file """

		ext = os.path.splitext(filename)[-1]

		#Currently only working with pdf
		#filename = filename.replace(ext, ".pdf")	#hack

		if ext == "jpg" :
			filename[-3:] = "png"
			warnings.warn("The 'jpg' format is not supported for motif image. Type is set tp 'png'")

		#self.gimme_obj.to_img(filename)	
		logo = self.create_logo()	#returns a logo object
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

	def create_logo(self, ax = None, motif_len=None):
		""" Creates motif logo in axes object """

		# convert to pandas dataframe
		df = pd.DataFrame(self.counts).transpose()
		df.columns = ["A", "C", "G", "T"]

		if not motif_len:
			motif_len = df.shape[0]

		# transform matrix to information based values
		info_df = logomaker.transform_matrix(df, from_type="counts", to_type="information")

		# create Logo object
		logo = logomaker.Logo(info_df, ax = ax)

		# style
		logo.style_xticks(rotation=0, fmt='%d', anchor=0)
		logo.ax.set_ylim(0, 2)
		logo.ax.set_xlim(-0.5, motif_len-0.5)
		logo.ax.set_yticks([0, 0.5, 1, 1.5, 2], minor=False)
		logo.style_spines(visible=False)
		logo.style_spines(spines=['left', 'bottom'], visible=True)
		logo.ax.xaxis.set_ticks_position('none')
		logo.ax.xaxis.set_tick_params(pad=-1)

		return logo

###########################################################

def get_motif_format(content):
	""" Get motif format from string of content """
	
	#Estimate input format
	if re.match(r"MEME version.+", content, re.DOTALL) is not None: # MOTIF\s.+letter-probability matrix.+[\d\.\s]+", content, re.MULTILINE) is not None:
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
