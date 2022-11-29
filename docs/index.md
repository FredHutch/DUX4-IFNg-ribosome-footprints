--- 
title: "DUX4 ribosome footprints profiling and translation efficiency"
author: "Chao-Jen Wong"
date: "2022-11-29"
site: bookdown::bookdown_site
documentclass: book
bibliography: [book.bib, packages.bib]
url: https://FredHutch.github.io/DUX4-IFNg-ribosome-footprints
#cover-image: images/cover.png
description: |
  This is a book providing detailed descriptions and R code to serve the purposed of transparency and reproducibility of our ribosome footprints analysis from the publication: DUX4 orchestrates translational reprograming by broadly suppressing translation.
biblio-style: apalike
csl: chicago-fullnote-bibliography.csl
---

# About

This book provides detailed descriptions and R code serving the purposes of transparency and reproducibility of the ribosome footprints analysis from our publication: _DUX4 orchestrates translational reprograming_.

## Samples and treatments
The processed data sets of the RNA-seq and Ribo-seq data supporting the manuscript reside at the [data](https://github.com/FredHutch/DUX4-IFNg-ribosome-footprints/data) directory. Both RNA-seq and Ribo-seq data sets are transformed into `DESeq2DataSet` instances and consist of triplicates of untreated, IFN-gamma, DUX4, and DUX4+INF-gamma treatments. Together, they are used find insights into the translational reprogramming of DUX4. 

## Softwear requirement
* R (>=4.0.3): the tidyverse project, knitr, bookdown, rmarkdown
* Bioconductor: DESeq2, goseq, GenomicAlignment, GenomicFeature, ribosomeProfilingQC, and etc.

## data
The [data](https://github.com/FredHutch/DUX4-IFNg-ribosome-footprints/data) directory in this repo includes the processed data sets used for our analyses:

* p_sites.rda: a list of _GRagnes_ instances of processed p-sites for 12 Ribo-seq samples    
* dds_cds.rda: a list of _DESeqDataSet_ instances containing p-site counts on CDS regions and metadata (size factor and treatment)    
* dds_.rda: DESeqDataSet instances containing p-site counts on different genomic features inclduing 5' UTR and 3' UTR    

## TxDb annotation package for gencode v35

Show how we make the annotation TxDb package here

## Additional scripts
The [scripts](https://github.com/FredHutch/DUX4-IFNg-ribosome-footprints/scripts) directory contains the R code and shell scripts performing the bioinformatics analysis for the manuscripts.
