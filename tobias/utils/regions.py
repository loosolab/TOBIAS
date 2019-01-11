#!/usr/bin/env python

"""
Classes for working with genomic regions

@author: Mette Bentsen
@contact: mette.bentsen (at) mpi-bn.mpg.de
@license: MIT

"""

import numpy as np
import sys
import re
from copy import deepcopy
import pyBigWig
from collections import Counter

#--------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------#

class OneRegion(list):

	nuc_to_pos = {"A":0, "T":1, "C":2, "G":3}

	def __init__(self, lst=["",0,0]):

		super(OneRegion, self).__init__(iter(lst))
		no_fields = len(lst)

		#Required
		self.chrom = lst[0]
		self.start = int(lst[1])		#exclude start
		self.end = int(lst[2])			#include end

		#Optional
		self.name = lst[3] if no_fields > 3 else ""
		self.score = lst[4] if no_fields > 4 else ""
		self.strand = lst[5] if no_fields > 5 else "."

	def __str__(self):
		return("{0}".format("\t".join(str(x) for x in self)))
	
	def tup(self):
		return((self.chrom, self.start, self.end, self.strand))
	
	def update(self):

		self[0] = self.chrom
		self[1] = self.start
		self[2] = self.end


	def get_length(self):
		return(self.end - self.start)

	def get_width(self):
		return(self.end - self.start)

	def extend_reg(self, bp):
		""" Extend region with bp to either side """

		self.start -= bp
		self[1] = self.start

		self.end += bp
		self[2] = self.end

		return(self)

	def set_width(self, bp):
		""" Set with of region centered on original region """

		flank_5, flank_3 = int(np.floor(bp/2.0)), int(np.ceil(bp/2.0))  #flank 5', flank 3'

		if self.strand == "-":
			mid = int(np.ceil((self.start + self.end) / 2.0))
			self.start = mid - flank_5
			self.end = mid + flank_3
		else:
			mid = int(np.floor((self.start + self.end) / 2.0))
			self.start = mid - flank_3
			self.end = mid + flank_5

		self[1] = self.start
		self[2] = self.end

		return(self)


	def split_region(self, max_length):
		""" Split genomic region into smaller subsets. Returns a RegionList of OneRegion objects """

		regions = RegionList() 	#Empty regions object

		starts = list(range(self.start, self.end, max_length))
		ends = list(range(self.start + max_length, self.end, max_length)) + [self.end]

		for (start, end) in zip(starts, ends):
			regions.append(OneRegion([self.chrom, start, end]))

		return(regions)


	def check_boundary(self, boundaries_dict, action="cut"):
		""" Check if regions are within chromosome boundaries. Actions:
				- "cut": cut region to bounds
				- "remove": remove region outside bounds (returns None)
		"""

		if (action == "cut" or action == "remove"):

			chrom = self.chrom

			#Check start
			if self.start < 0:
				if action == "cut":
					self.start = 0
				elif action == "remove":
					self = None

			#Check end
			if self.end > int(boundaries_dict[chrom]):
				if action == "cut":
					self.end = int(boundaries_dict[chrom])

				elif action == "remove":
					self = None

			if self.get_length() < 0:
				self = None

		else:
			exit("Error in regions.check_boundary: unknown action")

		return(self)


	def get_signal(self, pybw, numpy_bool = True):
		""" Get signal from bigwig in region """

		try:
			values = pybw.values(self.chrom, self.start, self.end, numpy=numpy_bool)
			values = np.nan_to_num(values)	#nan to 0
			
			if self.strand == "-":
				signal = values[::-1]
			else:
				signal = values
				
		except:
			print("Error reading region: {0} from pybigwig object".format(self.tup()))
			signal = np.zeros(self.end - self.start)
			signal[:] = np.nan

		return(signal)	


#--------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------#

class RegionList(list):
	""" Class Regions is a list of OneRegion objects """

	def __str__(self):
		return("\n".join([str(oneregion) for oneregion in self]))
	
	def __init__(self, lst=[]):
		super(RegionList, self).__init__(iter(lst))

	# ---- I/O ---- #
	def from_list(self, lst):
		""" Initialize object from list of OneRegion objects """

		for obj in lst:
			self.append(obj)
		return(self)


	def from_bed(self, bedfile_f):
		""" Initialize Object from bedfile """

		#Read all lines
		bedlines = open(bedfile_f).readlines()
		self = RegionList([None]*len(bedlines))
		for i, line in enumerate(bedlines):
		
			if line.startswith("#"): #comment lines are excluded
				continue

			#Test line format
			if re.match("[^\s]+\t\d+\t\d+.", line) == None:
				print("ERROR: Line {0} in {1} is not proper bed format:\n{2}".format(i+1, bedfile_f, line))
				sys.exit()

			columns = line.rstrip().split("\t")
			columns[1] = int(columns[1]) #start
			columns[2] = int(columns[2]) #end
			
			if columns[1] >= columns[2]:
				print("ERROR: Line {0} in {1} is not proper bed format:\n{2}".format(i+1, bedfile_f, line))
				sys.exit()

			region = OneRegion(columns)
			self[i] = region

		return(self)


	def as_bed(self, additional=True):

		bed = ""
		for region in self:

			#Update info
						
			line = "{0}\n".format("\t".join([str(col) for col in region]))

			#columns = [region.chrom, region.start, region.end, region.name, region.score, region.strand]

			#if additional == True:
			#line = "{0}\n".format("\t".join([str(col) for col in columns]))
			#else:
			#line = "{0}\n".format("\t".join([str(col) for col in columns[:3]]))
			bed += line

		return(bed)


	def write_bed(self, bed_f):
		""" Write regions to bedfile """

		bed_content = self.as_bed()
		out = open(bed_f, "w")
		out.write(bed_content)
		out.close()


	# ---- Working with regions ---- #
	def count(self):
		""" Number of regions """
		return(len(self))


	def loc_sort(self, contig_list=[]):
		""" Sorts list of region objects based on genomic location """
		if len(contig_list) > 0:
			order_dict = dict(zip(contig_list, range(len(contig_list))))
			self.sort(key=lambda region: (order_dict.get(region.chrom, region.chrom), region.start, region.end, region.name))
		else:
			self.sort(key=lambda region: (region.chrom, region.start, region.end, region.name))

	def score_sort(self):
		""" Sorts list of region objects based on score in last column """
		self.sort(key=lambda region: float(region[-1]), reverse=True)

	def get_chroms(self):

		chroms = {}
		for region in self:
			chroms[region.chrom] = ""
			
		return(list(chroms.keys()))


	def split_chrom(self):
		""" Split list of regions into separate chromosomes (output is RegionLists in dict) """

		output_dict = {}
		for region in self:
			if region.chrom not in output_dict:
				output_dict[region.chrom] = RegionList()

			output_dict[region.chrom].append(region)

		return(output_dict)


	def chunks(self, n):
		""" Split list of regions into n chunks (e.g. to use in multiprocessing) """

		no_reg = self.count()

		if no_reg > 0:
			per_chunk = int(np.ceil(no_reg/float(n)))
			chunk_lst = [RegionList(self[i:i+per_chunk]) for i in range(0, no_reg, per_chunk)]
		else:
			chunk_lst = [self]
		return(chunk_lst)


	def merge(self, name=False): 	
		""" Merge overlapping genomic regions. If name == True, regions will only be merged if the name is the same (used in count_overlaps) """

		self.loc_sort()		#sort before overlapping
		no = self.count()	#number of regions

		prev_chrom, prev_start, prev_end = "", 0, 0
		i = 1
		while i < no:

			prev_chrom, prev_start, prev_end = self[i-1].chrom, self[i-1].start, self[i-1].end
			curr_chrom, curr_start, curr_end = self[i].chrom, self[i].start, self[i].end

			#Check naming
			merge_flag = True
			if name == True:
				merge_flag = True if self[i-1].name == self[i].name else False

			if curr_chrom == prev_chrom and curr_start < prev_end and merge_flag == True:		#if overlapping regions
					self[i].start = prev_start
					self[i][1] = prev_start

					self[i].end = curr_end
					self[i][2] = curr_end

					del self[i-1]
					no -= 1
			else:
				i += 1  #compare next

		return(self)

	def remove_chroms(self, chromlist):
		""" Removes regions within chromlist chromosomes """

		no = self.count()
		i = 0
		while i < no:
			if self[i].chrom in chromlist:
				del self[i]
				no -= 1
			else:
				i += 1

	def keep_chroms(self, chromlist):
		no = self.count()
		i = 0
		while i < no:
			if self[i].chrom not in chromlist:
				del self[i]
				no -= 1
			else:
				i += 1

	def subset(self, no):
		""" Take no-size subset of regions from RegionList """

		if self.count() > no:
			del self[no:]
		return(self)


	def remove_duplicates(self):
		""" Returns a unique list of sites """

		self.loc_sort()
		prev_chrom, prev_start, prev_end, prev_strand = "", 0, 0, ""
		unique = RegionList()

		for i, region in enumerate(self):

			curr_chrom, curr_start, curr_end, curr_strand = region.chrom, region.start, region.end, region.strand

			if curr_chrom == prev_chrom and curr_start == prev_start and curr_end == prev_end and curr_strand == prev_strand:  #not unique
				pass
			else:
				unique.append(region)

			#Save as previous for next comparison
			prev_chrom, prev_start, prev_end, prev_strand = region.chrom, region.start, region.end, region.strand

		return(unique)

	"""
	def intersect(self, b_regions):
		# Overlap regions with other regions - only return those regions in a with an overlap in b

		#Sort before starting
		self.loc_sort()
		b_regions.loc_sort()

		a_len = self.count()
		b_len = b_regions.count()

		chromlist = sorted(list(set([region.chrom for region in self] + [region.chrom for region in b_regions])))
		chrom_pos = {chrom:chromlist.index(chrom) for chrom in chromlist}

		a = self
		b = b_regions

		a_i = 0
		b_i = 0

		while a_i < a_len and b_i < b_len:

			a_chrom, a_start, a_end = a[a_i].chrom, a[a_i].start, a[a_i].end
			b_chrom, b_start, b_end = b_regions[b_i].chrom, b_regions[b_i].start, b_regions[b_i].end

			#Do comparison
			if a_chrom == b_chrom:



			#a is ahead of b -> check next b
			elif chrom_pos[a_chrom] > chrom_pos[b_chrom]: 	#if a_chrom is after current b_chrom
				b_i += 1

			#b is ahead of a -> check next a
			elif chrom_pos[b_chrom] > chrom_pos[a_chrom]:	#if b_chrom is after current a_chrom
				a_i += 1
	"""


	def subtract(self, b_regions):
		""" Subtract b_regions from self regions """

		#Sort before starting
		self.loc_sort()
		b_regions.loc_sort()

		a_len = self.count()
		b_len = b_regions.count()

		chromlist = sorted(list(set([region.chrom for region in self] + [region.chrom for region in b_regions])))
		chrom_pos = {chrom:chromlist.index(chrom) for chrom in chromlist}

		a = self
		b = b_regions

		a_i = 0
		b_i = 0
		while a_i < self.count() and b_i < b_regions.count():

			a_chrom, a_start, a_end = a[a_i].chrom, a[a_i].start, a[a_i].end
			b_chrom, b_start, b_end = b_regions[b_i].chrom, b_regions[b_i].start, b_regions[b_i].end

			if a_chrom == b_chrom:

				if a_end <= b_start:	#current a is placed before current b
					a_i += 1	

				elif a_start >= b_end:	#current a is placed after current b 
					b_i += 1

				elif a_start >= b_start and a_end <= b_end: 	#a is completely contained within b (a is removed)
					del a[a_i]

				elif a_start < b_start and a_end > b_end: 		#b is completely contained within a (a is split into two)
					a[a_i] = OneRegion([a_chrom, a_start, b_start])				#first new region
					a.insert(a_i + 1, OneRegion([a_chrom, b_end, a_end]))			#second new region

				elif a_end > b_end: 		#delete left side
					a[a_i] = OneRegion([a_chrom, b_end, a_end])

				elif a_start < b_start: 	#delete right side
					a[a_i] = OneRegion([a_chrom, a_start, b_start])
			
			elif chrom_pos[a_chrom] > chrom_pos[b_chrom]: 	#if a_chrom is after current b_chrom
				b_i += 1

			elif chrom_pos[b_chrom] > chrom_pos[a_chrom]:	#if b_chrom is after current a_chrom
				a_i += 1

		return(self)

	def apply_method(self, method, *args):
		""" Applies a method to every OneRegion object in regions list """ 

		reglist = RegionList()
		for i in range(self.count()):
			out = method(self[i], *args) 	#apply method to region

			if out == None:			#no return, change OneRegion object itself
				pass
			else:
				if type(out) is OneRegion:
					reglist.append(out)
				elif type(out) is RegionList:
					reglist.extend(out)

		return(reglist)


	def resolve_overlaps(self, priority="higher"):
		""" Priority "higher" or "lower" """

		self.loc_sort()
		no_regions = len(self)

		i = 0
		j = 1

		while i + j < no_regions:

			reg_a = self[i]
			reg_b = self[i+j]

			if reg_a == None:
				i += 1
				j = 1
			elif reg_b == None:
				j += 1
			else:
				if reg_a.chrom == reg_b.chrom:
					if reg_b.start < reg_a.end:
						scores = [reg_a.score, reg_b.score]

						if priority == "higher":
							worst = scores.index(min(scores))
						else:
							worst = scores.index(max(scores))

						if worst == 0:
							self[i] = None
						else:
							self[i+j] = None
					else:
						i += 1
						j = 1
				else:
					i += 1
					j = 1

		#Remove all None
		non_overlapping = RegionList()
		for reg in self:
			if reg != None:
				non_overlapping.append(reg)

		return(non_overlapping)


	def get_signal_dict(self, bigwigs):
		""" Get dict of signal[region.tup][bigwig] = signal """

		signal_dict = {region.tup():{bigwig:[] for bigwig in bigwigs} for region in self}
		for bigwig in bigwigs:
			pybw = pyBigWig.open(bigwig, "rb")
			for region in self:
				signal_dict[region.tup()][bigwig] = region.get_signal(pybw)
			pybw.close()

		return(signal_dict)


	def get_width_distri(self):

		c = Counter()
		for region in self:
			reglen = region.get_width()
			c[reglen] += 1
		return(c)


	
	def count_overlaps(self):
		""" Returns a dictionary of strings and tuples - string-keys represent total bp per TF, and tuples represent the overlap in bp between two TFs """

		#Join all with similar name first
		self.loc_sort()
		self = self.merge(name=True) 

		#Count overlap
		i = 0  #index of current site
		j = 1  #index of next/compared site

		no_sites = len(self)
		
		overlap = {}

		while i < no_sites:

			s1_chrom, s1_start, s1_end, s1_id = self[i].chrom, self[i].start, self[i].end, self[i].name

			reglen = s1_end - s1_start
			overlap[s1_id] = overlap.get(s1_id,0) + reglen

			#Compare to next sites
			flag = 1	#overlapping
			j = 1
			while flag == 1 and i+j < no_sites:

				s2_chrom, s2_start, s2_end, s2_id = self[i+j][:4]

				#Do the regions overlap?
				if s1_chrom == s2_chrom:
					if s1_start < s2_end and s1_end > s2_start+1:	#if overlap
						bp_overlap = min([s1_end, s2_end]) - max([s1_start, s2_start]) 
						overlap[(s1_id, s2_id)] = overlap.get((s1_id, s2_id), 0) + bp_overlap
						overlap[(s2_id, s1_id)] = overlap.get((s2_id, s1_id), 0) + bp_overlap
						j += 1
					else:
						flag = 0
				else:
					flag = 0
			i += 1

		return(overlap)



#--------------------------------------------------------------------------------------------------#
#------------------------------ Additional functions related to regions ---------------------------#
#--------------------------------------------------------------------------------------------------#

def overlap_to_distance(overlap_dict):
	""" """

	#Find all TF names
	names = [key for key in overlap_dict.keys() if isinstance(key, str)]
	names = sorted(names)
	no_ids = len(names)

	similarity_mat = np.zeros((no_ids, no_ids))

	for i, id1 in enumerate(names): #rows
		for j, id2 in enumerate(names): #columns
			if i != j:
				tot_overlap = overlap_dict.get((id1,id2), 0)	#Find key-pair otherwise no overlap

				tot_id1 = overlap_dict[id1]
				tot_id2 = overlap_dict[id2]

				id1_dist = 1 - tot_overlap / float(tot_id1)
				id2_dist = 1 - tot_overlap / float(tot_id2)

				dist = min([id1_dist, id2_dist])
				similarity_mat[i,j] = dist

	return((similarity_mat, names))