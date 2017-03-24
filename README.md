--- 

# motifmatchr

[![Build Status](https://travis-ci.org/GreenleafLab/motifmatchr.svg?branch=master)](https://travis-ci.org/GreenleafLab/motifmatchr)

## Introduction

motifmatchr is an R package for fast motif matching, using C++ code from the MOODS library. The MOODS library was developed by Pasi Rastas, Janne Korhonen, and Petri Martinmäki. The core C++ library from MOODs version MOODS 1.9.3 code has been included in this repository. 

## Installation

Installation is easiest using the devtools package. The function `install_github` will install the package.

``` r
devtools::install_github("GreenleafLab/motifmatchr")
```

A number of needed packages are installed in this process. One of the dependencies has a system requirement for the gsl library, so if this is not installed already it may need to be installed separately.  

## match_motifs

The primary method of motifmatchr is `match_motifs`.  This method has two mandatory arguments:

1) Position weight matrices or position frequency matrices, stored in the PWMatrix, PFMatrix, PWMatrixList, or PFMatrixList objects from the TFBSTools package

2) Either a set of genomic ranges (GenomicRanges or RangedSummarizedExperiment object) or a set of sequences (either DNAStringSet, DNAString, or simple character vector)

If the second argument is a set of genomic ranges, a genome sequence is also required. If the genomic ranges include seqinfo, by default the genome specified in the seqinfo will be used (if the relevant BSgenome package is installed). Otherwise you can supply either a short string specifying the genome build if the corresponding BSgenome object is installed, a BSgenone object, a DNAStringSet object, or a FaFile object pointint to a fasta file.  

The method can return three possible outputs, depending on the `out` argument:

1) (Default, with `out = "matches"`) Boolean matrix indicating which ranges/sequences contain which motifs, stored as "matches" in assays slot of SummarizedExperiment object

2) (`out = "scores"`) Same as (1) plus two additional assays -- a matrix with the score of the high motif score within each range/sequence (score only reported if match present) and a matrix with the number of motif matches.

3) (`out = "positions"`) A GenomicRangesList with the ranges of all matches within the input ranges/sequences. 

## Quickstart

```{r}
library(motifmatchr)
library(GenomicRanges)

# load some example motifs
data(example_motifs, package = "motifmatchr") 

# Make a set of peaks
peaks <- GRanges(seqnames = c("chr1","chr2","chr2"),
                 ranges = IRanges(start = c(76585873,42772928,100183786),
                                  width = 500))

# Get motif matches for example motifs in peaks
motif_ix <- match_motifs(example_motifs, peaks, genome = "hg19") 
motif_matches(motif_ix) # Extract matches matrix from SummarizedExperiment result

# Get motif positions within peaks for example motifs in peaks 
motif_ix <- match_motifs(example_motifs, peaks, genome = "hg19",
                         out = "positions") 
```

## More information

For a more detailed overview, see [vignette](https://greenleaflab.github.io/motifmatchr/articles/motifmatchr.html).
