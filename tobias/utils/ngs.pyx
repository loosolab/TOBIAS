# cython: language_level=3

"""
Classes to work with NGS data (reads, readlists, etc.)

@author: Mette Bentsen
@contact: mette.bentsen (at) mpi-bn.mpg.de
@license: MIT	
"""

import pysam
from collections import Counter

from tobias.utils.regions import *
import time
import sys

import numpy as np
cimport numpy as np
import cython


#--------------------------------------------------------------------------------------------------#
class OneRead:
	""" OneRead class - takes many attributes from pysam read class"""

	def __init__(self):

		#self.read = 

		#New info on read
		self.cutsite = None		#no cutsite yet
		self.kmer = []			#kmer in integer nucleotide code around cutsite
		self.bias = 1  			#value of PSSM scoring for this read (changed during bias assignment)
		self.weight = 1			#the count-weight of this read (changed during correction)

	def __str__(self):
		return("OneRead object:\n{0}".format(vars(self)))    #)name: {0}\nchrom: {1}\nstart: {2}\nend: {3}\n".format(self.query_name, self.chrom, self.reference_start, self.reference_end))

	def from_read(self, read):
		""" Fill with attributes from read """

		#Taken from AlignedSegment attributes
		self.chrom = read.reference_name
		self.query_name = read.query_name
		self.reference_start = read.reference_start
		self.reference_end = read.reference_end
		self.query_alignment_start = read.query_alignment_start
		self.query_alignment_end = read.query_alignment_end
		self.template_length = read.template_length
		self.is_reverse = read.is_reverse 
		self.query_length = read.query_length
		self.flag = read.flag
		self.cigartuples = read.cigartuples

		return(self)

	def get_cutsite(self, read_shift):
		""" Finds the cutsite position from read taking direction & read shift into account """

		#Full length of read including soft clipped bases
		pos_shift, neg_shift = read_shift
		positions = [self.reference_start - self.query_alignment_start, self.reference_end + self.query_length - self.query_alignment_end]
		
		#- strand
		if self.is_reverse == True:
			self.cutsite = positions[1] + 1 + neg_shift  	#Position to the right of read start denotes cutsite
			
		#+ strand
		else:				
			self.cutsite = positions[0] + 1 + pos_shift		#.bam start is 0 based, so +1 converts to true genomic coordinates. 



	def get_kmer(self, genomic_sequence, k_flank):
		""" Extract kmer from genomic_sequence around cutsite """

		seq_start, seq_end = genomic_sequence.region.start, genomic_sequence.region.end		#excluding start, including end
		cutsite = self.cutsite

		if cutsite > seq_start + k_flank + 1 and cutsite <= seq_end - k_flank:

			#reverse read
			if self.is_reverse == True:
				i = seq_end - cutsite  		#position in revcomp_num
				self.kmer = genomic_sequence.revcomp[i-k_flank:i+k_flank+1] 

			#forward read
			else:
				i = cutsite - 1 - seq_start  	#-1 for 0-based python coordinates	
				self.kmer = genomic_sequence.sequence[i-k_flank:i+k_flank+1]

		else:
			self.kmer = np.array([4] * (k_flank*2 + 1))

	def shift_cutsite(self, bp):

		if self.is_reverse == True:
			self.cutsite -= bp
		else:
			self.cutsite += bp


#-------------------------------------------------------------#
class ReadList(list):
	""" List of OneRead objects """

	def __init__(self, data_array=[]):
		list.__init__(self, data_array)

	def __getitem__(self, sliced):
		return ReadList(self[sliced])

	def from_bam(self, bam_obj, region):
		try:
			for read in bam_obj.fetch(region.chrom, region.start, region.end):
				if read.is_unmapped == False and read.is_duplicate == False: 
					self.append(OneRead().from_read(read))

		except:
			sys.exit("Error reading {0} from bam object".format(region))

		return(self)
		

	def split_strands(self):

		strand_reads = [ReadList(), ReadList()] 	#forward and reverse
		for read_obj in self:
			if read_obj.is_reverse == True:
				strand_reads[1].append(read_obj)
			else:
				strand_reads[0].append(read_obj)
		return(strand_reads)


	def get_insert_sizes(self):
		""" Extract insert sizes from reads in dict """

		insert_sizes = Counter({})
		for read in self:

			size = abs(read.template_length) - 9
			insert_sizes[size] = insert_sizes.get(size, 0) + 1

		return(insert_sizes)

	@cython.boundscheck(False)
	@cython.wraparound(False) 	#dont deal with negative indices
	def signal(self, region):

		cdef int reg_start = region.start
		cdef int reg_end = region.end
		cdef int cut
		cdef int reg_len = reg_end - reg_start
		cdef np.ndarray[np.float64_t, ndim=1] values = np.zeros(reg_len)		#array to fill in

		for read_obj in self:
			cut = read_obj.cutsite
			if cut > reg_start and cut <= reg_end:
				values[cut - reg_start - 1] += read_obj.weight

		return(values)