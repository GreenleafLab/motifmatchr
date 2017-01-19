motifmatchr
========

motifmatch is an R package for fast motif matching, using C++ code from the MOODS library. The MOODS library was developed by Pasi Rastas, Janne Korhonen, and Petri Martinmäki. The core C++ library from MOODs version MOODS 1.9.3 code has been included in this repository.   

Installation
------------

Installation is easiest using the devtools package. The function `install_github` will install the package.

``` r
devtools::install_github("GreenleafLab/motifmatchr", auth_token = "my_token")
```

The argument auth\_token takes in your github [personal access token](https://github.com/settings/tokens). This token is needed because at the moment this repository is private.

A number of needed packages are installed in this process. Note that for functions that require a genome sequence, the package [BSgenome.Hsapiens.UCSC.hg19](https://bioconductor.org/packages/release/data/annotation/html/BSgenome.Hsapiens.UCSC.hg19.html) is used as a default argument. If using another genome build, the appropraiate BSgenome object for your species should be passed to functions requiring a genome build (e.g. `match_pwms`).

Depending on your repository settings, the Bioconductor dependencies may fail to install. Use `setRepositories(graphics=F)` to see what repositories you have activated and to add the BioC software repository if need be.

Quickstart
-------------------

```{r}
library(motifmatchr)
library(GenomicRanges)

# load some example motifs
data(example_motifs, package = "motifmatchr") 

# Make a set of peaks
peaks <- GRanges(seqnames = c("chr1","chr2","chr2"),
                 ranges = IRanges(start = c(76585873,42772928,100183786),
                                  width = 500))

# Get motif matches for example motifs in peaks (using hg19 genome, the default)
motif_ix <- match_pwms(example_motifs, peaks) 

# Get motif positions within peaks for example motifs in peaks 
motif_ix <- match_pwms(example_motifs, peaks, out = "positions") 
```
