This is the work-in-progress repo for the ribosome footprints and traslation efficiency analysis supporitng the manuscript _DUX4 orchestrates translational reprograming by broadly suppressing translation efficiency_. This repo hosts a [gitbook](https://FredHutch.github.io/DUX4-IFNg-ribosome-footprints), bulit with [bookdown](https://bookdown.org/), containing detailed description of the analysis and R code that reproduces the results and figures on the fly.  

## Samples and treatments
The translation efficiency was calculated by Ribo-seq and RNA-seq samples, and each of the platform consisted with triplicates of untreated, IFN-gamma, DUX4, and DUX4+INF-gamma treated.

## Softwear requirement
* R (>=4.0.3): the tidyverse project, knitr, bookdown, rmarkdown
* Bioconductor: DESeq2, goseq, GenomicAlignment, ribosomeProfilingQC, and many others

## Directories
There are five major directories:

### data 
This directory includes processed files:
* p_sites.rda: a list of _GRagnes_ instances of processed p-sites for 12 Ribo-seq samples
* dds_cds.rda: a list of _DESeqDataSet_ instances containing p-site counts on CDS regions and metadata (size factor and treatment)
* dds_.rda: DESeqDataSet instances containing p-site counts on different genomic features inclduing 5' UTR and 3' UTR


### scripts
This directory contains all the R code and shell scripts performing the bioinformatics analysis for the manuscripts. 

### gitbook
This directory includes the orginal R markdown and related files that made the gitbook. The code in the Rmd files were orgainzed, readable and can reproduce the results and figures for the manuscript.

### docs
This directory for the [gitbook](https://FredHutch.github.io/DUX4-IFNg-ribosome-footprints), made of the HTML pages randered by the files in the _gitbook_ directory and [bookdown](https://bookdown.org).


