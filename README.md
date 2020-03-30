TOBIAS - Transcription factor Occupancy prediction By Investigation of ATAC-seq Signal 
=======================================

[![PyPI Version](https://img.shields.io/pypi/v/tobias.svg?style=plastic)](https://pypi.org/project/tobias/)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=plastic)](http://bioconda.github.io/recipes/tobias/README.html)
[![bioRxiv badge](https://img.shields.io/badge/bioRxiv-10.1101%2F869560-blue?style=plastic)](https://www.biorxiv.org/content/10.1101/869560v2)

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
$ pip install tobias
```

TOBIAS is also available as a conda package on the Bioconda channel:
```bash
$ conda install tobias -c bioconda
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

Overview and command-line examples
-------------

* [ATACorrect](https://github.com/loosolab/TOBIAS/wiki/ATACorrect): Bias correction of ATAC-seq reads in open chromatin
* [ScoreBigwig](https://github.com/loosolab/TOBIAS/wiki/ScoreBigwig): Calculate footprint scores from corrected cutsites
* [BINDetect](https://github.com/loosolab/TOBIAS/wiki/BINDetect): Estimation of differentially bound motifs based on scores, sequence and motifs
* [PlotAggregate](https://github.com/loosolab/TOBIAS/wiki/PlotAggregate): Plot aggregated ATAC-seq signals in combinations of .bed/.bw to visualize footprints
* [PlotHeatmap](https://github.com/loosolab/TOBIAS/wiki/PlotHeatmap): Plot heatmaps and aggregates of ATAC-seq signals in combinations of .bed/.bw to visualize footprints
* [PlotTracks](https://github.com/loosolab/TOBIAS/wiki/PlotTracks): Plot IGV-style genomic signals such as cutsites and footprints across a selection of regions
* [FormatMotifs](https://github.com/loosolab/TOBIAS/wiki/FormatMotifs): A utility to convert and join/split across different motif-file formats
* [ClusterMotifs](https://github.com/loosolab/TOBIAS/wiki/Additional) : Cluster motifs and create consensus motifs based on similarity
* [CreateNetwork](https://github.com/loosolab/TOBIAS/wiki/CreateNetwork): Create TF-TF binding network from annotated TFBS
* [FilterFragments](https://github.com/loosolab/TOBIAS/wiki/Additional): Filter fragments from a .bam-file using a .bed-file of regions
* [Additional utility tools](https://github.com/loosolab/TOBIAS/wiki/Additional)


Pipelines
----------------
While each TOBIAS tool can be run independently, they are developed to be run as part of an analysis pipeline. We provide ready-made pipelines for performing bias-correction, footprinting, differential binding and visualization for multiple conditions automatically.

**Snakemake pipeline**  
We provide a pre-set snakemake workflow which is found [here](https://github.molgen.mpg.de/loosolab/TOBIAS_snakemake).

**Nextflow pipeline**  
You can also run the TOBIAS tool as a nextflow pipeline. The pre-set workflow can be found [here](https://github.molgen.mpg.de/loosolab/TOBIAS-nextflow).

**Nextflow kubernetes/de.NBI cloud aware pipeline**  
We also provide the TOBIAS nextflow pipeline for a cloud computing environment. One version utilizes a [kubernetes framework](https://github.molgen.mpg.de/loosolab/TOBIAS-nextflow/tree/master/TOBIAS_MAPOKS), and a second version utilizing a webbased job scheduler, started automatically within a local TOBIAS run, making use of the de.NBI [cloud](https://github.molgen.mpg.de/loosolab/TOBIAS-nextflow/tree/master/TOBIAS_MACSEK).


License
------------
This project is licensed under the [MIT license](LICENSE). 


Contact
------------
Mette Bentsen (mette.bentsen (at) mpi-bn.mpg.de)
