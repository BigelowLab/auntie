---
title: "Auntie"
output: github_document
---

`auntie` provides eDNA specific functionality for woirking with common eDNA related projects.  This package differes from [charlier](https://github.com/BigelowLab/charlier) which provides generic functionality for working on the charlie computer.

`auntie` is named in honor of [Auntie Edna](https://duckduckgo.com/?q=Incredibles+Auntie+Edna)

![Image of Auntie Edna](https://static.wikia.nocookie.net/pixar/images/8/8b/I2_-_Edna.png/revision/latest?cb=20180621022628)



## Requirements and Installation

 + [R v4+](https://www.r-project.org/)
 
 + [rlang](https://CRAN.R-project.org/package=rlang)

 + [dplyr](https://CRAN.R-project.org/package=dplyr)

 + [Biostrings](http://www.bioconductor.org/packages/release/bioc/html/Biostrings.html)

 + [ShortRead](http://www.bioconductor.org/packages/release/bioc/html/ShortRead.html)

 + [charlier](https://github.com/BigelowLab/charlier)
 
 
```
# for github install
remotes::install_github("BigelowLab/auntie, up = FALSE)"

# for local install (preferred for eDNA group members)
devtools::install("/mnt/storage/data/edna/packages/auntie", up = FALSE)
```

## Functionality

### Filepairs

Many operations depend upon `forward/reverse` file pairs.  We have developed some functions for listing and testing file pairings. 


### [reformat](https://github.com/BioInfoTools/BBMap/blob/a9ceda047a7c918dc090de0fdbf6f924292d4a1f/sh/reformat.sh) wrapper

We provide a convenience wrapper, `bbmap_subsample()` to subsample fastq files using `reformat.sh`

### [cutadapt]() wrapper

Not implemented yet.
