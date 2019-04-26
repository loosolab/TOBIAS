#!/usr/bin/env python

"""
Classes for working with motifs and scanning with moods

@author: Mette Bentsen
@contact: mette.bentsen (at) mpi-bn.mpg.de
@license: MIT

"""

import numpy as np
import copy

import MOODS.scan
import MOODS.tools
import MOODS.parsers

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.text import TextPath
from matplotlib.patches import PathPatch
from matplotlib.font_manager import FontProperties

#Internal 
from tobias.utils.regions import * 
from tobias.utils.utilities import filafy 	#filafy for filenames

#----------------------------------------------------------------------------------------#
#List of OneMotif objects
class MotifList(list):

	def __init__(self, lst=[]):

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

	def setup_moods_scanner(self):

		tups = [(motif.prefix, motif.strand, motif.pssm, motif.threshold) for motif in self] 		#list of tups
		if len(tups) > 0:
			self.names, self.strands, self.matrices, self.thresholds = list(map(list, zip(*tups))) 	#get "columns"
		else:
			self.names, self.strands, self.marices, self.thresholds = ([], [], [], [])

		scanner = MOODS.scan.Scanner(7)
		scanner.set_motifs(self.matrices, self.bg, self.thresholds)

		self.moods_scanner = scanner


	def scan_sequence(self, seq, region):

		if self.moods_scanner == None:
			self.setup_moods_scanner()

		#Scan sequence
		results = self.moods_scanner.scan(seq)

		#Convert results to RegionList
		sites = RegionList()
		for (matrix, name, strand, result) in zip(self.matrices, self.names, self.strands, results):
			motif_length = len(matrix[0])
			for match in result:
				start = region.start + match.pos 	#match pos is 1 based
				end = start + motif_length		
				score = round(match.score, 5)

				site = OneRegion([region.chrom, start, end, name, score, strand])
				sites.append(site)

		return(sites)

#----------------------------------------------------------------------------------------#
#Contains info on one motif formatted for use in moods
class OneMotif:

	bases = ["A", "C", "G", "T"]

	def __init__(self, motifid, name, counts):
		
		self.id = motifid			#must be unique
		self.name = name	#does not have to be unique

		self.prefix = "" 	#output prefix set in set_prefix
		self.counts = counts  	#counts
		self.strand = "+"		#default strand is +

		#Set later
		self.pfm = None
		self.bg = np.array([0.25,0.25,0.25,0.25]) 	#background set to equal by default
		self.pssm = None 							#pssm calculated from get_pssm
		self.threshold = None 						#threshold calculated from get_threshold
		self.bits = None

	def __str__(self):

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
		self.pfm = self.counts / np.sum(self.counts, axis=0)


	def get_reverse(self):
		if self.pfm is None:
			self.get_pfm()

		reverse_motif = copy.deepcopy(self)
		reverse_motif.strand = "-"
		reverse_motif.pfm = MOODS.tools.reverse_complement(self.pfm,4)
		return(reverse_motif)	#OneMotif object


	def get_pssm(self, ps=0.01):
		""" """
		if self.pfm is None:
			self.get_pfm()

		bg_col = self.bg.reshape((-1,1))
		pseudo_vector = ps * bg_col

		pssm = np.log(np.true_divide(self.pfm + pseudo_vector, np.sum(self.pfm + pseudo_vector, axis=0))) - np.log(bg_col)
		pssm = tuple([tuple(row) for row in pssm])
		self.pssm = pssm


	def get_threshold(self, pvalue):
		if self.pssm is None:
			self.get_pssmm()

		self.threshold = MOODS.tools.threshold_from_p(self.pssm, self.bg, pvalue, 4)
		return(self)


	def calc_bit_score(self):

		if self.pfm is None:
			self.get_pfm()

		pfm_arr = np.copy(self.pfm)
		pfm_arr[pfm_arr == 0] = np.nan

		#Info content per pos
		entro = pfm_arr * np.log2(pfm_arr)
		entro[np.isnan(entro)] = 0
		info_content = 2 - (- np.sum(entro, axis=0))		#information content per position in motif
		self.bits = self.pfm * info_content	


	def plot_logo(self):

		LETTERS = { "T" : TextPath((-0.305, 0), "T", size=1, prop=fp),
		            "G" : TextPath((-0.384, 0), "G", size=1, prop=fp),
		            "A" : TextPath((-0.35, 0), "A", size=1, prop=fp),
		            "C" : TextPath((-0.366, 0), "C", size=1, prop=fp) }
		COLOR_SCHEME = {'G': 'orange', 
		                'A': "#CC0000", 
		                'C': 'mediumblue', 
		                'T': 'darkgreen'}

		def add_letter(base, x, y, scale, ax):
			""" Add letter to axis at positions x/y"""
			
			text = LETTERS[base]
			t = mpl.transforms.Affine2D().scale(1*globscale, scale*globscale) + \
				mpl.transforms.Affine2D().translate(x,y) + ax.transData
			p = PathPatch(text, lw=0, fc=COLOR_SCHEME[base], transform=t)
			if ax != None:
				ax.add_artist(p)
			return p

		self.calc_bit_score()
		self.length = self.bits.shape[1]

		fp = FontProperties(family='sans-serif', weight="bold") 
		globscale = 1.35
	
		#Plot logo
		fig, ax = plt.subplots(figsize=(10,3))	
		max_y = 0
		for position in range(self.length):	#0-based positions

			base_bit_tups = zip(OneMotif.bases, self.bits[:,position])
			
			#Go through bases sorted from lowest to highest bit score
			y = 0	#position to place letter
			for (base, score) in sorted(base_bit_tups, key=lambda tup: tup[1]):
				add_letter(base, position+1, y, score, ax)
				y += score
			max_y = max(max_y, y)

		plt.xticks(range(1,self.length+1))
		plt.xlim((0.2, self.length+0.8)) 
		plt.ylim((0, max_y)) 
		plt.tight_layout()  	
		
		return(fig, ax)	


###########################################################

def get_motif_format(content):

	#Estimate input format
	motif_format = "unknown"
	
	if re.match("MEME version.+", content, re.DOTALL) is not None: # MOTIF\s.+letter-probability matrix.+[\d\.\s]+", content, re.MULTILINE) is not None:
		motif_format = "meme"

	elif re.match(">.+A.+\[", content, re.DOTALL) is not None:
		motif_format = "jaspar"

	elif re.match(">.+", content, re.DOTALL) is not None:
		motif_format = "pfm"

	return(motif_format)


def convert_motif(content, output_format):
	""" Output formats are "pfm", "jaspar" or "meme" """

	bases = ["A", "C", "G", "T"]
	input_format = get_motif_format(content)
	converted_content = ""

	if input_format == output_format:

		#remove any meme headers 
		m = re.match("^(MEME.*?)(MOTIF.*)", content, re.DOTALL)
		if m:
			converted_content = m.group(2) + "\n"
		else:	
			converted_content = content + "\n"

	################ pfm <-> jaspar ################
	elif (input_format == "pfm" or input_format == "jaspar") and (output_format == "pfm" or output_format == "jaspar"):
		
		for line in content.split("\n"):
			if line.startswith(">"):
				converted_content += line + "\n"	#header line + \n as this was removed in split
				i = -1
			
			else:
				m = re.match(".*?([\d]+[\d\.\s]+).*?", line)

				if m:
					i += 1	#i is 0 for first pfm line
					pfm_line = m.group(1)
					fields = [field for field in pfm_line.rstrip().split()]
					
					converted_line =  "{0} [ {1} ] \n".format(bases[i], "\t".join(fields)) if output_format == "jaspar" else "\t".join(fields) + "\n"
					converted_content += converted_line

					if i == 3: # last line
						converted_content += "\n"

				else:
					continue


	################ meme -> jaspar/pfm ################
	elif input_format == "meme" and (output_format == "jaspar" or output_format == "pfm"): 
				
		motif_content = []
		header = ""
		lines = content.split("\n") + ["MOTIF"]		#add motif to end to write out motif
		for idx, line in enumerate(lines):
			if line.startswith("MOTIF"):
				
				#Write any previous motif saved
				if len(motif_content) > 0:
					for i, column in enumerate(motif_content):	#column = list of values
						converted_line = "{0} [ {1} ] \n".format(bases[i], "\t".join(column)) if output_format == "jaspar" else "\t".join(column) + "\n"
						converted_content += converted_line 	#contains \n

				#Get ready for this motif 
				if idx < len(lines) - 1:		#Up until the last line, it is possible to save for next	
					columns = line.strip().split()
					if len(columns) > 2: #MOTIF, ID, NAME
						motif_id, name = columns[1], columns[2]
					elif len(columns) == 2: # MOTIF, ID
						motif_id, name = columns[1], columns[1]

					header = ">{0}\t{1}\n".format(motif_id, name)

					converted_content += header
					motif_content = [[] for _ in range(4)] 	#ACGT

			elif re.match("^[\s]*([\d\.\s]+)$", line):	#starts with any number of spaces (or none) followed by numbers
				columns = line.rstrip().split()
				for i, col in enumerate(columns):
					motif_content[i].append(col)

		

	################ jaspar/pfm -> meme ################
	elif (input_format == "jaspar" or input_format == "pfm") and output_format == "meme":
			
		motif_content = [] 	#no motifs found yet, this is empty

		lines = content.split("\n") + [">"] 	#add ">" at the end to make sure that the last motif is saved
		for idx, line in enumerate(lines):

			m = re.match(".*?([\d]+[\d\.\s]+).*?", line)

			if line.startswith(">"):
				
				#Write any previous motif saved
				if len(motif_content) > 0:
					motif_w = len(motif_content[0])		
					n_sites = int(round(sum(float(motif_content[i][0]) for i in range(4)), 0)) 	#sum of first site freqs 
					
					converted_content += "letter-probability matrix: alength=4 w={0} nsites={1} E=0\n".format(motif_w, n_sites)
					for i in range(motif_w):
						row = [float(motif_content[j][i]) for j in range(4)] 	#row contains original row from content
						n_sites = round(sum(row), 0)
						row_freq = ["{0:.5f}".format(num/n_sites) for num in row] 
						converted_content += "  ".join(row_freq) + "\n"
					converted_content += "\n"

				if idx < len(lines) - 1:		#Up until the last line, it is possible to save for next
					columns = line[1:].strip().split()		#[1:] to remove > from header
					try:
						motif_id, name = columns[0], columns[1]
					except:
						motif_id, name = ".", "."
						print(line)

					converted_content += "MOTIF {0} {1}\n".format(motif_id, name)
					motif_content = [] 	#list of rows from jaspar format motif

			elif m:
				columns = [field for field in m.group(1).rstrip().split()]
				motif_content.append(columns)

	return(converted_content)


def pfm_to_motifs(content):
	""" Content of a pfm motif file to MotifList format """

	#Read motifs to moods
	pfm_names = []
	pfms = []
	idx = -1

	motiflist = MotifList([])
	for line in content.split("\n"):
		if line.startswith(">"):

			#Read name for this motif
			columns = line.replace(">", "").rstrip().split()
			motifid, alt_name = columns[0], columns[1]
			
			motif_obj = OneMotif(motifid, alt_name, [])		#pfm is set to empty list

			motiflist.append(motif_obj)

		elif len(motiflist) > 0:  #if at least one header line was found
			m = re.match(".*?([\d]+[\d\.\s]+).*?", line)

			if m:
				pfm_line = m.group(1)
				pfm_fields = [float(field) for field in pfm_line.rstrip().split()]
				motif_obj.counts.append(pfm_fields)
			else:
				continue

	#check correct format of pfms
	for motif in motiflist:
		rows, cols = np.array(motif.counts).shape	
		if rows != 4:
			sys.exit("ERROR: Motif {0} has an unexpected format and could not be read")

	return(motiflist)

