#!/usr/bin/env python

"""
Skript to compare two sets of motifs with each other. Generating similarity matrix and a clustered heatmap.
If only one motif file is given, it will be compared with itself.

@author: René Wiegandt
@contact: rene.wiegandt (at) mpi-bn.mpg.de
"""

import argparse
import numpy as np
import seaborn as sns
import pandas as pd
import datetime
import os
import yaml
import re
import math
import itertools
import warnings
from Bio import motifs
from matplotlib import pyplot as plt
from gimmemotifs.motif import Motif,read_motifs
from gimmemotifs.comparison import MotifComparer
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
import scipy.spatial.distance as ssd
from collections import defaultdict

from tobias.utils.utilities import *
from tobias.utils.logger import *


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
    description += "Usage:\nTOBIAS MotifClust --motifs1 <motifs1.jaspar>\n\n"
    parser.description = format_help_description("MotifClust", description) 

    parser._action_groups.pop() #pop -h

    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    visualisation = parser.add_argument_group('visualisation arguments')

    required.add_argument("-m1", "--motifs1", dest="motifs1", required=True, help="A motif file containing a set of motifs", metavar="")
    
    optional.add_argument("-m2", "--motifs2", dest="motifs2", help="A motif file containing a set of motifs (common: Motif Database like jaspar or hocomoco)", metavar="")
    optional.add_argument("-f", "--format", dest="Format", choices= ['pwm', 'transfac', 'xxmotif', 'jaspar', 'minimal', 'meme', 'align'], help="Format of the first motif file [‘pwm’, ‘transfac’, ‘xxmotif’, ‘jaspar’, ‘minimal’, ‘meme’ or ‘align’] (Default: ‘jaspar’)", default="jaspar")
    optional.add_argument("-f2", "--format_2", dest="Format_2", choices= ['pwm', 'transfac', 'xxmotif', 'jaspar', 'minimal', 'meme', 'align'], help="Format of the second motif file [‘pwm’, ‘transfac’, ‘xxmotif’, ‘jaspar’, ‘minimal’, ‘meme’ or ‘align’] (Default: ‘jaspar’)", default="jaspar")
    optional.add_argument("-t", "--threshold", dest="threshold", help="Cluster threshold (Default = 0.5)", type=float, default=0.5)
    optional.add_argument("-s", "--ds_threshold", dest="d_thresh", help="Threshold for dissimilar cutoff", type=float, default=0.5)
    optional.add_argument("-l", "--method", dest="method", help="Method for clustering (See: https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html)", default="average")
    optional.add_argument("-a", "--cons_format", dest="cons_format", choices= ['transfac', 'meme', 'pwm'], help="Format of consensus motif file [‘transfac’, ‘meme’, ‘pwm’] (Default: pwm)", default="pwm")
    optional.add_argument("-p", "--prefix", dest="name", help="Output prefix (Default: ‘motif_comparison’)", default="motif_comparison")
    optional.add_argument("-o", "--outdir", dest="out", help="Output directory (Default: ‘./‘)", default="./")
    optional.add_argument("-m", "--merge", dest="merge", help="Merge both motif files and compare 'all to all'", action="store_true")
    optional.add_argument("-cc", "--no_col_clust", dest="ncc", help="No column clustering", action="store_false")
    optional.add_argument("-rc", "--no_row_clust", dest="nrc", help="No row clustering", action="store_false")
    optional.add_argument("-z", "--z_score", dest="zscore", choices= ['row', 'col', 'None'], help="Calculate the z-score for row or column [‘col’, ‘row’, ‘None’] (Default: None)", default="None")
    
    visualisation.add_argument("-nh", "--no_heatmap", dest="no_heatmap", help="Disable heatmap", action="store_false")
    visualisation.add_argument("-e", "--type", dest="type", choices= ['png', 'pdf', 'jpg'], help="Plot file type [png, pdf, jpg] (Default: pdf)", default="pdf")
    visualisation.add_argument("-x", "--width", dest="width", help="Width of Heatmap (Default: autoscaling)", type=int)
    visualisation.add_argument("-y", "--height", dest="height", help="Height of Heatmap (Default: autoscaling)", type=int)
    visualisation.add_argument("-d", "--dpi", dest="dpi", help="Dpi for plots (Default: 100)", type=int, default=100)
    visualisation.add_argument("-c", "--color_palette", dest="color", help="Color palette (All possible paletts: https://python-graph-gallery.com/197-available-color-palettes-with-matplotlib/. Add '_r' to reverse palette.)", default="YlOrRd_r")
    args = parser.parse_args()
    return args


#--------------------------------------------------------------------------------------------------------#
def get_motifs(path, format):
    """Read motif file via gimmemotifs or biopython. 
    If motif is read via biopython the bio.motif instace gets converted to an gimmemotif instance.

    Parameter:
    ----------
    path : String
        a path to a motif file
    format : String
        format of motif file

    Returns:
    --------
    list
        a list of motif instances
    """


    gimme_formats = ["pwm", "transfac", "xxmotif", "align"]

    # Read motif via gimmemotifs if format is supported
    if(format in gimme_formats):
        return read_motifs(infile = path, fmt = format, as_dict = False)
    else: 
        # Otherwise Bio.python is used
        with open(path) as f:
            bio_motif_list = motifs.parse(f, format)

        # Get full motif name if format is minimal
        if format == "minimal":
            name_list = get_minimal_names(path)

        # Converting Bio.motif instance to gimmemotif motif instance
        gimme_motif_list = list() # list of all motifs
        current_motif = 0
        for bio_motif in bio_motif_list:
            motif_rows = list()
            for pos_id in range(0, bio_motif.length-1):
                row = list() # each row represents one motif index ( A C G T )
                for letter in range(0,4):
                    row.append(bio_motif.counts[letter][pos_id])
                motif_rows.append(row)
            gimme_motif = Motif(motif_rows) # generate gimmemotif motif instance
            # Add motif name
            if format == "minimal":
                gimme_motif.id = name_list[current_motif]
            else:
                gimme_motif.id = bio_motif.name
            gimme_motif_list.append(gimme_motif)
            current_motif+=1
    return gimme_motif_list


#--------------------------------------------------------------------------------------------------------#
def get_minimal_names(path):
    """Get list of full motif names from minimal meme files

    Parameter:
    ----------
    path : String
        a path to a motif file

    Returns:
    --------
    List
        a list containing all motif names
    """

    motif_names = list()
    with open(path) as f:
        for line in f.readlines():
            if line.startswith("MOTIF"): # get motif header
                motif_names.append(re.sub("MOTIF\s+", "", line).rstrip())
    return motif_names


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

    similarity_dict = dict()
    
    m2_keys = list(score_dict.values())[0].keys()
    m2_label = [s.replace('\t', ' ') for s in m2_keys] # replacing tabs with whitespace

    for m1, m2_dict in score_dict.items():
        m1_space = m1.replace('\t', ' ') # replacing tabs with whitespace
        sim_list = list()
        for score_list in m2_dict.values():
            score = float('%.3f'%(1 - score_list[0]))
            sim_list.append(score)
        similarity_dict[m1_space] = sim_list
    
    return pd.DataFrame(similarity_dict, index = m2_label).replace(-0, 0)


#--------------------------------------------------------------------------------------------------------#
def subset_matrix(matrix, name_list_1, name_list_2):
    """Subsetting matrix into two symetric matrices given two label lists

    Parameter:
    ----------
    matrix :  Pandas DataFrame
        Full similarity matrix
    name_list_1 : list
        list of names
    name_list_2 : list
        list of names

    Returns:
    --------
    list:
        list with subsets of matrix
    """

    sim_matrix_1 = matrix[name_list_1].loc[name_list_1]
    sim_matrix_2 = matrix[name_list_2].loc[name_list_2]
    sim_matrix = matrix[name_list_1].loc[name_list_2]
    return sim_matrix_1, sim_matrix_2, sim_matrix


#--------------------------------------------------------------------------------------------------------#
def get_dissimilar_motifs(matrix, threshold):
    """Get motifs that are not similar to any of the motifs compared to them.

    Parameter:
    ----------
    matrix: DataTable
        a DataTable contaning the similarity scores of the motifs.
    threshold : float
        threshold for similarity score
    
    Return:
    ---------
    dissimilar_motifs : list
        list of motifs that are not similar to to any of the motifs compared to them.
    """

    columns = matrix.shape[1]
    dissimilar_motifs = list()
    for col in range(columns-1):
        col_vec = matrix.iloc[:,col] > threshold
        if(col_vec.eq(True).all()):
            dissimilar_motifs.append(matrix.columns[col])
    return dissimilar_motifs


#--------------------------------------------------------------------------------------------------------#
def get_gc_content(motifs):
    """Get the GC content and the length for a list of motifs

    Parameter:
    ----------
    motifs : list
        A list of gimmemotifs motif objects.
    
    Return:
    -------
    motif_stats : dict
        A dictionary containing the motif name as key and
        a list with GC content and motif length as value.
    """

    motif_stats = dict()
    for m in motifs:
        # Remove first row
        m_woh = re.sub('^>.*\n','',m.to_pwm())
        m_pwm = m_woh.split()
        motif_length = len(m)
        gc_sum = 0
        # Get sum of G and C ratios
        for x in range(0, len(m_pwm)-1, 4):
            gc_sum += float(m_pwm[x+1]) + float(m_pwm[x+2])
        # Calc GC ratio of the motif
        gc_content = gc_sum/motif_length
        motif_name = m.id.replace('\t', ' ')
        motif_stats[motif_name] = [round(gc_content,4), motif_length]
    return motif_stats


#--------------------------------------------------------------------------------------------------------#
def write_motif_stats(stats, out_file):
    """Write stats to file

    Parameter:
    ----------
    stats: dict
        dictionary with motif name as key and list with stats as value
    out_file: string
        path and name of out-file
    """

    with open(out_file, "w") as f: 
        f.write("Motif\tGC content\tMotif length\n")
        for motif_name,stats in stats.items():
            stats_string = "\t".join(map(str, stats))
            f.write("{}\t{}\n".format(motif_name, stats_string))


#--------------------------------------------------------------------------------------------------------#
def clustering_motifs(similarity_matrix, threshold, method, subset_1=None, subset_2=None):
    """Clusters motif similarity matrix hierarchicaly

    Parameter:
    ----------
    similarity_matrix : DataTable
        a DataTable contaning the similarity score
    threshold : float
        clustering threshold
    method : string
        clustering method used by scipy.cluster.hierarchy.linkage
    subset_1 : DataTable
        subset of similarity matrix 
    subset_2 : DataTable
        subset of similarity matrix 
    
    Returns:
    --------
    List : ndarray, ndarray, dict, dict 
        The hierarchical clustering of rows and cols encoded as a linkage matrix.
        A dictionary containing named clusters for rows and cols.
    """

    # Clustering
    if subset_1 is not None and subset_2 is not None:
        sub_1_vector = ssd.squareform(subset_1.to_numpy())
        sub_2_vector = ssd.squareform(subset_2.to_numpy())
        col_linkage = linkage(sub_1_vector, method=method)
        row_linkage = linkage(sub_2_vector, method=method)
    else:
        vector = ssd.squareform(similarity_matrix.to_numpy())
        row_linkage = linkage(vector, method=method)
        col_linkage = row_linkage

    # Flatten cluster
    row_labels = fcluster(row_linkage, threshold, criterion="distance")
    col_labels = fcluster(col_linkage, threshold, criterion="distance")

    # Extract cluster
    if subset_1 is not None and subset_2 is not None:
        row_cluster = get_cluster(row_labels, subset_2.columns)
        col_cluster = get_cluster(col_labels, subset_1.columns)
    else:
        row_cluster = get_cluster(row_labels, similarity_matrix.index.values)
        col_cluster = get_cluster(col_labels, similarity_matrix.columns)

    return row_linkage, col_linkage, row_cluster, col_cluster


#--------------------------------------------------------------------------------------------------------#
def get_cluster(labels, label_names):
    """Map clusternumber to label names

    Parameter:
    ----------
    labels : ndarray
        An array of length n. T[i] is the flat cluster number.
    label_names : list
        List of names
    
    Return:
    -------
    cluster : dict
        A dictionary containing named clusters.
    """

    cluster = defaultdict(list)
    for i in range(len(labels)):
        cluster["Cluster_" + str(labels[i])].append(label_names[i])
    return cluster


#--------------------------------------------------------------------------------------------------------#
def scaling(axis_len):
    """Scaling of plot figure size

    Parameter:
    ----------
    axis_len : int
        Number of rows or columns

    Return:
    ----------
    scale : float
        scale for axis
    """

    return 0.5/math.log10(axis_len)


#--------------------------------------------------------------------------------------------------------#
def write_yaml(clusters, name):
    """Writes cluster dict to yaml

    Parameter:
    ----------
    clusters : dict
        dictionary containing lists as values.
    name : string
        path and prefix for outfile
    """

    yml_out = name + "_cluster.yml"
    with open(yml_out, 'w') as outfile:
        yaml.dump(dict(clusters), outfile, default_flow_style=False)


#--------------------------------------------------------------------------------------------------------#
def plot_dendrogram(label, linkage, font_size, out ,name, threshold, y, dpi, t):
    """Plot dendrogram with highlighted threshold 

    Parameter:
    ----------
    label : list
        List of labels
    linkage : ndarray
        The hierarchical clustering of rows or cols encoded as a linkage matrix.
    font_size : int
        font size
    out : String
        Output path
    name : String
        Plot title
    threshold : float
        dendrogram cluster threshold
    y : int
        yaxis size of plot
    dpi : int
        dpi of plot
    t : string
        type of plot file
    """

    plt.figure(figsize=(20, y))
    plt.title(name, fontsize=20)
    plt.axvline(x=threshold, color="red")
    dendrogram(linkage, color_threshold=threshold, labels=label, leaf_font_size=font_size, orientation="right")
    try:
        plt.savefig(out + "/" + name + "_dendrogram." + t, dpi=dpi)
    except ValueError as e:
        print("Skipped ploting of Heatmap.")
        print("Error: " + str(e))


#--------------------------------------------------------------------------------------------------------#
def generate_consensus_motifs(motif_list, cluster, similarity_dict):
    """Generates a consensus motif for each cluster

    Parameter:
    ----------
    motif_list : list
        list of all motif instances
    cluster : dict
        A dictionary containing named clusters.
    similarity_dict : dict
        Output from gimmemotif.get_all_scores

    Returns:
    --------
        dict of consensus motifs as gimmemotif motif instance
    """

    consensus_motifs = dict()
    for cluster, value in cluster.items():
        cluster_motifs = list()
        # Get motif instances
        for m in motif_list:
            if m.id in value:
                cluster_motifs.append(m)
        consensus_motifs[cluster] = generate_consensus_motif(cluster_motifs, similarity_dict)
    return consensus_motifs


#--------------------------------------------------------------------------------------------------------#
def generate_consensus_motif(motifs, similarity_dict):
    """Recursive algorithm to generate consensus motif from multiple motifs

    Parameter:
    ----------
    motifs : list
        list of gimmemotif motif instances
    similarity_dict : dict
        Output from gimmemotifs.get_all_scores

    Returns:
    --------
    gimmemotifs motif instance
        consensus motif of all motifs in motifs list
    """

    mc = MotifComparer()

    if len(motifs) == 1:
        return motifs[0]

    # score cannot be higher lower that 0 so max score is initialized with 0
    max_score = 0
    for motif_1, motif_2  in itertools.combinations(motifs, 2):
        score = similarity_dict[motif_1.id][motif_2.id]
        if score[0] >= max_score:
            max_score = score[0]
            merge_1 = motif_1
            merge_2 = motif_2

    # generate average motif
    average_motif = merge_1.average_motifs(merge_2, pos=score[1], orientation=score[2])
    average_motif.id = merge_1.id + "_" + merge_2.id

    # alter list 
    motifs.append(average_motif)
    motifs.remove(merge_1)
    motifs.remove(merge_2)

    sim_dict = mc.get_all_scores(motifs, motifs, match = "total", metric = "pcc", combine = "mean" )

    cons_motif = generate_consensus_motif(motifs, sim_dict)
    return cons_motif


#--------------------------------------------------------------------------------------------------------#
def save_motif_image(consensus_motif, prefix, ext, out):
    """Save motif to image file

    Paramater:
    ----------
    consensus_motif : gimmemotif motif instance
        The motif to be plotted
    prefix : string
        Prefix for image file
    ext : string
        File format
    out : string
        Output path
    """

    if ext == "jpg":
        ext = "png"
        warnings.warn("The 'jpg' format is not supported for motif image. Type is set tp 'png'")
    consensus_motif.to_img(os.path.join(out, prefix + "_consensus." + ext))


#--------------------------------------------------------------------------------------------------------#
def convert_motif(consensus_motif, f):
    """ Convert motif to given format

    Parameter:
    ----------
    consensus_motif : gimmemotif motif instance
        The motif to be written to file
    out : string
        Output string without file ending
    f : string
        Format of motif file

    Returns:
    --------
    String
        Motif in given format
    """

    if f == "meme":
        motif = consensus_motif.to_meme()
    elif f == "pwm":
        motif = consensus_motif.to_pwm()
    elif f == "transfac":
        motif = consensus_motif.to_transfac()
    else:
        raise ValueError("Format " + f + " is not supported")
    return motif


#--------------------------------------------------------------------------------------------------------#
def consesus_motif_out_wrapper(consensus_motifs, out_prefix, img_out, cons_format, name, typ, file_name = None):
    """ Wrapper for consensus motif output

    Parameter:
    ----------
    consensus_motifs : dict
        dictionary a consensur motif for each cluster
    out_prefix : string
        Output path + prefix
    img_out : string
        Path for images
    cons_format : string
        Format for motif file
    name : string
        Prefix for motif images
    file_name : string
        File name of clustered motifs
    typ : string
        File format of image
    """

    if file_name is not None:
        img_prefix = name + "_" + file_name
        cons_out_path = out_prefix + "_" + file_name + "_consensus_motifs." + cons_format
    else:
        img_prefix = name 
        cons_out_path = out_prefix + "_consensus_motifs." + cons_format

    cons_out = open(cons_out_path, 'w+')

    for cluster, consensus_motif in consensus_motifs.items():
        consensus_motif.id = cluster + "_" + consensus_motif.id

        # Save motif as image
        save_motif_image(consensus_motif, img_prefix + "_" + cluster, typ, img_out)

        # Save motifs to file
        motif_string = convert_motif(consensus_motif, cons_format)
        cons_out.write(motif_string + "\n\n")


#--------------------------------------------------------------------------------------------------------#
def plot_heatmap(similarity_matrix, out, x, y, col_linkage, row_linkage, dpi, x_name, y_name, color, col_clust, row_clust, zscore):
    """Plots clustered similarity matrix as a heatmap

    Parameter:
    ----------
    similarity_matrix :  DataTable
        a DataTable contaning the similarity score
    out : string
        Prefix to output file
    x : int
        width for heatmap plot
    y : int
        height for heatmap plot
    col_linkage : ndarray
        The hierarchical clustering of cols encoded as a linkage matrix.
    row_linkage : ndarray
        The hierarchical clustering of rows encoded as a linkage matrix.
    dpi : int
        dpi of plot
    x_name : string
        label for x-axis
    y_name : string
        label for y-axis
    color : string
        colorpalette
    col_clust : boolean
        if True, cluster the columns
    row_clust : boolean
        if True, cluster the rows
    """

    vmin, vmax = None, None
    if zscore == "row":
        zs = 0
    elif zscore == "col":
        zs = 1
    else:
        zs = None
        vmin, vmax = 0, 1

    try:
        plot = sns.clustermap(similarity_matrix,
            row_linkage=row_linkage,
            col_linkage=col_linkage,
            col_cluster=col_clust,
            row_cluster=row_clust,
            z_score=zs,
            cmap=color,
            vmin=vmin, vmax=vmax,
            xticklabels=similarity_matrix.columns,
            yticklabels=similarity_matrix.index.values,
            figsize=(x,y))
        plot.ax_heatmap.set_xlabel(x_name)
        plot.ax_heatmap.set_ylabel(y_name)
        plot.savefig(out, bbox_inches='tight', dpi=dpi)
    except ValueError as e:
        print("Skipped ploting of Heatmap.")
        print("Error: " + str(e))


#--------------------------------------------------------------------------------------------------------#
def run_motifclust(args):

    print("entered run_motifclust")
    print(args)
    ###### Check input arguments ######
    check_required(args, ["motifs1"]) #Check input arguments
    check_files([args.motifs1, args.motifs2]) #Check if files exist
    make_directory(args.out)
    out_prefix = os.path.join(args.out, args.name)

    ###### Create logger and write argument overview ######
    logger = TobiasLogger("MotifClust.log")
    logger.begin()
    parser = add_motifclust_arguments(argparse.ArgumentParser())
    logger.arguments_overview(parser, args)
    #logger.output_files([])

    out_prefix = os.path.join(args.out, args.name)

    #---------------------------------------- Reading motifs from file(s) -----------------------------------#

    logger.info("Handling input file/files")

    motif_list = get_motifs(args.motifs1, args.Format)
    m1_names = list(m.id for m in motif_list)

    if args.motifs2:
        m2  = get_motifs(args.motifs2, args.Format_2)
        m2_names = list(m.id for m in m2)
        motif_list = motif_list + m2

    #---------------------------------------- Generating similarity matrix ----------------------------------#
    
    logger.info("Generating similarity matrix")

    mc = MotifComparer()

    # Pairwise comparison of a set of motifs compared to reference motifs
    score_dict = mc.get_all_scores(motif_list, motif_list, match = "total", metric = "pcc", combine = "mean" )

    # Generating similarity matrix
    similarity_matrix = generate_similarity_matrix(score_dict)

    if args.motifs2 and not args.merge :
        # Subsetting matrix
        sub_matrix_1, sub_matrix_2, similarity_matrix = subset_matrix(similarity_matrix, m1_names, m2_names)
    else:
        sub_matrix_1, sub_matrix_2 = None, None

    # Safe matrix to file
    matrix_out = out_prefix + "_matrix.txt"
    logger.info("Saving matrix to the file " + str(matrix_out))
    similarity_matrix.to_csv(matrix_out, sep = '\t')

    #---------------------------------------- Motif stats ---------------------------------------------------#

    logger.info("Making matrix statistics about dissimilar motifs and GC-content")
    
    dissimilar_motifs_out = out_prefix + "_dissimilar_motifs.txt"
    full_motifs_out = out_prefix + "_stats_motifs.txt"
    dissimilar_motifs = get_dissimilar_motifs(similarity_matrix, args.d_thresh)
    motifs_stats = get_gc_content(motif_list)
    write_motif_stats(motifs_stats, full_motifs_out)

    # Get motif stats only for dissimilar motifs
    stats_dissimilar_motifs = dict((k, motifs_stats[k]) for k in dissimilar_motifs)
    write_motif_stats(stats_dissimilar_motifs, dissimilar_motifs_out)
    
    #---------------------------------------- Motif clustering ----------------------------------------------#
    
    logger.info("Clustering motifs")

    row_linkage, col_linkage, row_cluster, col_cluster = clustering_motifs(similarity_matrix, args.threshold, args.method, sub_matrix_1, sub_matrix_2)

    # Scaling for plots
    y_len = len(similarity_matrix.index.values)
    x_len = len(similarity_matrix.columns)

    if args.height:
        y = args.height
    else:
        y = y_len*scaling(y_len)
    if args.width:
        x = args.width
    else:
        x = x_len*scaling(x_len)

    filename_1 = os.path.splitext(os.path.basename(args.motifs1))[0]
    f_name_2 = filename_1

    if args.merge or not args.motifs2:
        if args.ncc or args.nrc:
            if args.motifs2:
                filename_1 = os.path.splitext(os.path.basename(args.motifs1))[0] + "_" + os.path.splitext(os.path.basename(args.motifs2))[0]
                f_name_2 = filename_1
            write_yaml(col_cluster, out_prefix + "_" + filename_1)
            plot_dendrogram(similarity_matrix.columns, col_linkage, 12, args.out, args.name + "_" + filename_1, args.threshold, x, args.dpi, args.type)

    if args.motifs2 and not args.merge:
        filename_2 = os.path.splitext(os.path.basename(args.motifs2))[0]
        f_name_2 = filename_2
        if args.nrc: # if false skip writing yaml file and dendrogram for row clustering
            write_yaml(row_cluster, out_prefix + "_" + filename_2)
            plot_dendrogram(similarity_matrix.index.values, row_linkage, 12, args.out, args.name + "_" + filename_2, args.threshold, y, args.dpi, args.type)
        if args.ncc:
            write_yaml(col_cluster, out_prefix + "_" + filename_1)
            plot_dendrogram(similarity_matrix.columns, col_linkage, 12, args.out, args.name + "_" + filename_1, args.threshold, y, args.dpi, args.type)

    #---------------------------------------- Consensus motif -----------------------------------------------#

    logger.info("Building consensus motifs")

    # Image output path
    out_cons_img = os.path.join(args.out, "consensus_motifs_img")
    # Check if out_cons_img directory exists
    if not os.path.exists(out_cons_img):
        # Create if not
        os.mkdir(out_cons_img)

    # Generate output depending on set parameters
    if args.merge or not args.motifs2: # col and row are identical
        if args.ncc or args.nrc:
            # Generate consensus motifs for each column cluster 
            cons = generate_consensus_motifs(motif_list, col_cluster, score_dict)
            # Save output
            consesus_motif_out_wrapper(cons, out_prefix, out_cons_img, args.cons_format, args.name, args.type)
    else: # col and row differ
        if args.ncc:
            # Generate consensus motifs for each column cluster 
            col_cons = generate_consensus_motifs(motif_list, col_cluster, score_dict)
            # Save output
            consesus_motif_out_wrapper(col_cons, out_prefix, out_cons_img, args.cons_format, args.name, args.type, os.path.splitext(os.path.basename(args.motifs1))[0])
        if args.nrc:
            # Generate consensus motifs for each row cluster 
            row_cons = generate_consensus_motifs(motif_list, row_cluster, score_dict)
            # Save output
            consesus_motif_out_wrapper(row_cons, out_prefix, out_cons_img, args.cons_format, args.name, args.type, os.path.splitext(os.path.basename(args.motifs2))[0])

    #---------------------------------------- Plot heatmap --------------------------------------------------#

    if args.no_heatmap:
        pdf_out = out_prefix + "_heatmap." + args.type
        logger.info("Plotting the heatmap to the file " + str(pdf_out))
        plot_heatmap(similarity_matrix, pdf_out, x, y, col_linkage, row_linkage, args.dpi, filename_1, f_name_2, args.color, args.ncc, args.nrc, args.zscore)
    
    logger.end()


#--------------------------------------------------------------------------------------------------------#
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser = add_motifclust_arguments(parser)
    args = parser.parse_args()

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit()

    run_motifclust(args)
