TOBIAS - Transcription factor Occupancy prediction By Investigation of ATAC-seq Signal 
=======================================

[![PyPI Version](https://img.shields.io/pypi/v/tobias.svg?style=plastic)](https://pypi.org/project/tobias/)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/tobias/README.html)

Introduction 
------------

ATAC-seq (Assay for Transposase-Accessible Chromatin using high-throughput sequencing) is a sequencing assay for investigating genome-wide chromatin accessibility. The assay applies a Tn5 Transposase to insert sequencing adapters into accessible chromatin, enabling mapping of regulatory regions across the genome. Additionally, the local distribution of Tn5 insertions contains information about transcription factor binding due to the visible depletion of insertions around sites bound by protein - known as _footprints_. 

**TOBIAS** is a collection of command-line bioinformatics tools for performing footprinting analysis on ATAC-seq data, and includes:

<img align="right" width=150 src="/figures/tobias.png">

- Correction of Tn5 insertion bias
- Calculation of footprint scores within regulatory regions
- Estimation of bound/unbound transcription factor binding sites
- Visualization of footprints within and across different conditions

For information on each tool, please see the [wiki](https://github.com/loosolab/TOBIAS/wiki/).

Installation
------------
TOBIAS is written as a python package and can be quickly installed via pip:
```bash
pip install tobias
```

or directly from github using:
```bash
$ git clone https://github.com/loosolab/TOBIAS
$ cd TOBIAS
$ python setup.py install
```
Please see the [installation](https://github.com/loosolab/TOBIAS/wiki/installation) page for more info.

Usage
------------
All tools are available through the command-line as ```TOBIAS <TOOLNAME>```, for example:
``` 
$ TOBIAS ATACorrect
__________________________________________________________________________________________

                                   TOBIAS ~ ATACorrect
__________________________________________________________________________________________

ATACorrect corrects the cutsite-signal from ATAC-seq with regard to the underlying
sequence preference of Tn5 transposase.

Usage:
TOBIAS ATACorrect --bam <reads.bam> --genome <genome.fa> --peaks <peaks.bed>

Output files:
- <outdir>/<prefix>_uncorrected.bw
- <outdir>/<prefix>_bias.bw
- <outdir>/<prefix>_expected.bw
- <outdir>/<prefix>_corrected.bw
- <outdir>/<prefix>_atacorrect.pdf

(...)
```

Command-line examples
-------------
These examples use the test data provided in the [TOBIAS/test_data](https://github.com/loosolab/TOBIAS/tree/master/test_data) 
directory, so please make sure you are in the upper TOBIAS/ directory when running the commands. Larger files such as .bigwig and .bam-files are stored on an S3-storage server, which you get access to by using the download script provided in [S3_downloader](https://github.com/loosolab/TOBIAS/tree/dev/S3_downloader). NOTE: This download script requires "boto3" and "pyyaml" libraries, which are available in the supplied [conda environment](https://github.com/loosolab/TOBIAS/blob/dev/conda_env.yaml) or can be downloaded using pip (```$ pip install boto3 pyyaml```). 

Run the script with the config file to obtain the data:
```
$ python S3_downloader/S3_Downloader.py --yaml S3_downloader/s3_tobias_data_config.yaml
$ mv S3_download/data-tobias-2019/test_data/* test_data/
```

**ATACorrect: Bias correction of ATAC-seq reads in open chromatin**     
```
$ TOBIAS ATACorrect --bam test_data/Bcell_chr4.bam --genome test_data/genome_chr4.fa.gz --peaks test_data/merged_peaks.bed --blacklist test_data/blacklist_chr4.bed --outdir atacorrect_test --prefix Bcell --cores 8
```
_~ 3 minutes_ 

**FootprintScores: Calculate footprint scores from corrected cutsites**
```
$ TOBIAS FootprintScores --signal test_data/Bcell_corrected.bw --regions test_data/merged_peaks.bed --output Bcell_footprints.bw --cores 8
```
_~ 1 minute_

**BINDetect: Estimation of differentially bound motifs based on scores, sequence and motifs**   
```
$ TOBIAS BINDetect --motifs test_data/example_motifs.jaspar --signals test_data/Bcell_footprints.bw test_data/Tcell_footprints.bw --genome test_data/genome_chr4.fa.gz --peaks test_data/annotated_peaks.bed --peak_header test_data/annotated_peaks_header.txt --outdir bindetect_output --cond_names Bcell Tcell --cores 8
```
_~ 1 minute (two conditions, 7500 peaks, 86 motifs)_

**PlotAggregate: Plot aggregated ATAC-seq signals in combinations of .bed/.bw to visualize footprints**  

Visualize the difference in footprints between two conditions for all accessible sites:    
```
$ TOBIAS PlotAggregate --TFBS test_data/BATFJUN_all.bed  --signals test_data/Bcell_corrected.bw test_data/Tcell_corrected.bw --output BATFJUN_footprint_comparison_all.pdf --share_y both --plot_boundaries
```
_~ 3 seconds_

Visualize the difference in footprints between two conditions exclusively for bound sites:   
```
$ TOBIAS PlotAggregate --TFBS test_data/BATFJUN_Bcell_bound.bed test_data/BATFJUN_Tcell_bound.bed --signals test_data/Bcell_corrected.bw test_data/Tcell_corrected.bw --output BATFJUN_footprint_comparison_subsets.pdf --share_y both --plot_boundaries
```
_~ 5 seconds_

Visualize the split of bound/unbound sites for one condition:   
```
$ TOBIAS PlotAggregate --TFBS test_data/IRF1_all.bed test_data/IRF1_bound.bed test_data/IRF1_unbound.bed --signals test_data/Bcell_uncorrected.bw test_data/Bcell_expected.bw test_data/Bcell_corrected.bw --output IRF1_footprint.pdf  --share_y sites --plot_boundaries
```
_~ 8 seconds_

**PlotHeatmap: Plot heatmaps and aggregates of ATAC-seq signals in combinations of .bed/.bw to visualize footprints**   
```
$ TOBIAS PlotHeatmap --TFBS test_data/BATFJUN_Bcell_bound.bed test_data/BATFJUN_Bcell_unbound.bed --TFBS test_data/BATFJUN_Tcell_bound.bed test_data/BATFJUN_Tcell_unbound.bed --signals test_data/Bcell_corrected.bw test_data/Tcell_corrected.bw --output BATFJUN_heatmap.pdf --signal_labels Bcell Tcell --share_colorbar --sort_by -1
```
_~ 6 seconds_

**PlotTracks: Plot IGV-style genomic signals such as cutsites and footprints across a selection of regions**
```
$ TOBIAS PlotTracks --bigwigs test_data/*_corrected.bw --bigwigs test_data/*_footprints.bw --regions test_data/plot_regions.bed --sites test_data/binding_sites.bed --highlight test_data/binding_sites.bed --gtf test_data/transcripts_chr4.gtf --colors red darkblue red darkblue
```
_~ 10 seconds_

**FormatMotifs: A utility to convert and join/split across different motif-file formats**    
Join individual motif files to one:    
```
$ TOBIAS FormatMotifs --input test_data/individual_motifs/* --task join --output example_motifs.txt
```
_~ < 1 second_ 

Split a motif file containing several motifs:  
```
$ TOBIAS FormatMotifs --input test_data/example_motifs.jaspar --format meme --task split --output split_motifs
```
_~ < 1 second
Filter a larger motif file using TF names:
```
$ echo 'MAFK CTCF JUNB' > TF_names.txt
$ TOBIAS FormatMotifs --input test_data/example_motifs.jaspar --output filtered_motifs.txt --filter TF_names.txt
```
_~ < 1 second_

**Cluster Motifs: Cluster motifs and create consensus motifs based on similarity**
```
$ TOBIAS ClusterMotifs --motifs test_data/example_motifs.jaspar
```
_~ 1 minute_

**Filter fragments from a .bam-file using a .bed-file of regions:**
```
$ TOBIAS FilterFragments --bam test_data/Bcell_chr4.bam --regions test_data/merged_peaks.bed
```
_~ 30 seconds_

Snakemake pipeline
------------

You can run each TOBIAS tool independently or as part of a pipeline. We provide a pre-set snakemake workflow which is found [here](https://github.molgen.mpg.de/loosolab/TOBIAS_snakemake).

Nextflow pipeline
------------

You can also run the TOBIAS tool as a nextflow pipeline. The pre-set workflow can be found [here](https://github.molgen.mpg.de/loosolab/TOBIAS-nextflow).

kubernetes/de.NBI cloud aware pipeline
------------

We also provide the TOBIAS nextflow pipeline for a cloud computing environment. There is one version utilising a [kubernetes framework](https://github.molgen.mpg.de/loosolab/TOBIAS-nextflow/tree/master/TOBIAS_over_S3), and a second version utilizing a webbased job scheduler, started automatically within a local TOBIAS run, making use of the de.NBI [cloud](https://github.molgen.mpg.de/loosolab/TOBIAS-nextflow/tree/master/TOBIAS_over_NGINX).



License
------------
This project is licensed under the [MIT license](LICENSE). 


Contact
------------
Mette Bentsen (mette.bentsen (at) mpi-bn.mpg.de)
