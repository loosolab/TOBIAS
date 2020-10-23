#!/usr/bin/env python

"""
cluster_sites_by_overlap.py: Utility script to cluster .bed-file input and create a distance-matrix and dendrogram similar to the BINDetect output.

@author: Mette Bentsen
@contact: mette.bentsen (at) mpi-bn.mpg.de
@license: MIT
"""

import os
import argparse
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram
from tobias.utils.regions import RegionList, RegionCluster
from tobias.utils.logger import TobiasLogger
from tobias.utils.utilities import check_required, make_directory
from tobias.utils.logger import TobiasLogger

#----------------------------- Input arguments -----------------------------#

parser = argparse.ArgumentParser()
parser.add_argument("--bedfiles", nargs="*", help="Bedfiles with ID in the 4th column")
parser.add_argument("--outdir", help="Output directory (default: bedfile_clustering_output)", default="bedfile_clustering_output")
args = parser.parse_args()
check_required(args, ["bedfiles"])

logger = TobiasLogger()
make_directory(args.outdir) #Create dir if not existing

#------------------------- Read bedfiles to cluster ------------------------#

all_sites = RegionList()
for bedfile in args.bedfiles:
	with open(bedfile) as f:
		sites = RegionList().from_bed(bedfile)
		all_sites.extend(sites)

logger.info("Read {0} regions from {1} files".format(len(all_sites), len(args.bedfiles)))


#---------- Calculate distances between sites in RegionList object ---------#

logger.info("Calculating distances between regions")
overlaps = all_sites.count_overlaps()
clustering = RegionCluster(overlaps)
clustering.overlap_to_distance()

f_out = os.path.join(args.outdir, "distance_matrix.txt")
clustering.write_distance_mat(f_out)
logger.info("Wrote distance matrix to: {0}".format(f_out))


#------------------------------ Plot dendrogram ----------------------------#

logger.info("Plotting dendrogram")
clustering.cluster()	#provides linkage matrix

#Plot figure
fig = plt.figure(figsize = (5, len(clustering.names)/7))
dendro_dat = dendrogram(clustering.linkage_mat, 
						labels=clustering.names, 
						orientation="right", 
						above_threshold_color="black", 
						)

plt.title("Clustering of input regions")
plt.xlabel("Distance\n(based on overlap of .bed-file regions)")
plt.ylabel("Region names")

f_out = os.path.join(args.outdir, "dendrogram.pdf")
plt.savefig(f_out, bbox_inches="tight")
logger.info("Plotted dendrogram to: {0}".format(f_out))
logger.info("Done!")