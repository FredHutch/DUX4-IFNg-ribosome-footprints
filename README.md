To support the reproducibility and transparency of the computational work for the manuscript _DUX4 orchestrates translational reprogramming by broadly suppressing translation efficiency_, we made this repository to include our processed RNA-seq and Ribosome footprints sequencing (Ribo-seq) datasets, as well as a [__book__](https://FredHutch.github.io/DUX4-IFNg-ribosome-footprints) containing seven chapters demonstrating all the bioinformatic analysis and automatically executable R code (you can directly run the code and reproduce the statistics and figures in the manuscript).

## Samples and treatments
A total of 12 samples were prepared for RNA-seq and Ribo-seq, seperately, and each comprising four distinct treatments: untreated, DUX4-pulse, IFN-gamma, and DUX4-pulse+INF-gamma, with triplicates for each treatment. The processed data and profiling reside at the repository's [data](https://github.com/FredHutch/DUX4-IFNg-ribosome-footprints/data) directory. 

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
This directory for the [gitbook](https://FredHutch.github.io/DUX4-IFNg-ribosome-footprints), made of the HTML pages rendered by the files in the _gitbook_ directory and the [bookdown](https://bookdown.org) package.


