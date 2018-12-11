TOBIAS - Transcription factor Occupancy prediction By Investigation of ATAC-seq Signal 
=======================================

Introduction 
------------

ATAC-seq (Assay for Transposase-Accessible Chromatin using high-throughput sequencing) is a sequencing assay for investigating genome-wide chromatin accessibility. The assay applies a Tn5 Transposase to insert sequencing adapters into accessible chromatin, enabling mapping of regulatory regions across the genome. Additionally, the local distribution of Tn5 insertions contains information about transcription factor binding due to the visible depletion of insertions around sites bound by protein - known as _footprints_. 

**TOBIAS** is a collection of command-line bioinformatics tools for performing footprinting analysis on ATAC-seq data, and includes:

<img align="right" width=150 src="/figures/tobias.png">

- Correction of Tn5 insertion bias
- Calculation of footprint scores within regulatory regions
- Estimation of bound/unbound transcription factor binding sites
- Visualization of footprints within and across different conditions

For information on each tool, please see the [wiki](https://github.molgen.mpg.de/loosolab/TOBIAS/wiki/).

Installation
------------
TOBIAS is written as a python package and can be quickly installed within a conda environment using:
```bash
$ git clone https://github.molgen.mpg.de/loosolab/TOBIAS
$ cd TOBIAS
$ conda env create -f snakemake_pipeline/environments/tobias.yaml
$ conda activate TOBIAS_ENV
$ python setup.py install
```
Please see the [installation](https://github.molgen.mpg.de/loosolab/TOBIAS/wiki/installation) page for more info.

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

Snakemake pipeline
------------

You can run each TOBIAS tool independently or as part of a pipeline using the included snakemake workflow. Simply set the paths to required data within snakemake_pipeline/TOBIAS.config and run using:
```bash
$ cd snakemake_pipeline
$ conda activate TOBIAS_ENV
$ snakemake --snakefile TOBIAS.snake --configfile TOBIAS.config --cores [number of cores] --keep-going
```
For further info on setup, configfile and output, please consult the [wiki](https://github.molgen.mpg.de/loosolab/TOBIAS/wiki/snakemake-pipeline).

License
------------
This project is licensed under the [MIT license](LICENSE). 


Contact
------------
Mette Bentsen (mette.bentsen (at) mpi-bn.mpg.de)
