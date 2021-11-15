---
title: "Auntie"
output: github_document
---

`auntie` provides eDNA specific functionality for woirking with common eDNA related projects.  This package differes from [charlier](https://github.com/BigelowLab/charlier) which provides generic functionality for working on the charlie computer.

`auntie` is named in honor of [Auntie Edna](https://duckduckgo.com/?q=Incredibles+Auntie+Edna)

![Image of Auntie Edna](https://static.wikia.nocookie.net/pixar/images/8/8b/I2_-_Edna.png/revision/latest?cb=20180621022628){height="100"}



## Requirements and Installation

 + [R v4+](https://www.r-project.org/)
 
 + [rlang](https://CRAN.R-project.org/package=rlang)

 + [dplyr](https://CRAN.R-project.org/package=dplyr)

 + [Biostrings](http://www.bioconductor.org/packages/release/bioc/html/Biostrings.html)

 + [ShortRead](http://www.bioconductor.org/packages/release/bioc/html/ShortRead.html)

 + [charlier](https://github.com/BigelowLab/charlier)
 

### for github install
```
remotes::install_github("BigelowLab/auntie", up = FALSE)
```

### for local install (preferred for eDNA group members)
```
devtools::install("/mnt/storage/data/edna/packages/auntie", up = FALSE)
```

## Functionality

### Filepairs

Many operations depend upon `forward/reverse` file pairs.  We have developed some functions for listing and testing file pairings.


```r
library(auntie)
# create filepair listing (note use of new pipe v4.1+)
fp <- list_filepairs("/mnt/storage/data/edna/dada/example") |>
  verify_filepairs()
fp
```

```
## $forward
## [1] "/mnt/storage/data/edna/dada/example/E-EnvStd-01_S186_L001_R1_001.fastq.gz"
## [2] "/mnt/storage/data/edna/dada/example/E-EvStd_S358_L001_R1_001.fastq.gz"    
## [3] "/mnt/storage/data/edna/dada/example/E-HAB179_S370_L001_R1_001.fastq.gz"   
## 
## $reverse
## [1] "/mnt/storage/data/edna/dada/example/E-EnvStd-01_S186_L001_R2_001.fastq.gz"
## [2] "/mnt/storage/data/edna/dada/example/E-EvStd_S358_L001_R2_001.fastq.gz"    
## [3] "/mnt/storage/data/edna/dada/example/E-HAB179_S370_L001_R2_001.fastq.gz"
```

```r
interleave_filepairs(fp)
```

```
## [1] "/mnt/storage/data/edna/dada/example/E-EnvStd-01_S186_L001_R1_001.fastq.gz"
## [2] "/mnt/storage/data/edna/dada/example/E-EnvStd-01_S186_L001_R2_001.fastq.gz"
## [3] "/mnt/storage/data/edna/dada/example/E-EvStd_S358_L001_R1_001.fastq.gz"    
## [4] "/mnt/storage/data/edna/dada/example/E-EvStd_S358_L001_R2_001.fastq.gz"    
## [5] "/mnt/storage/data/edna/dada/example/E-HAB179_S370_L001_R1_001.fastq.gz"   
## [6] "/mnt/storage/data/edna/dada/example/E-HAB179_S370_L001_R2_001.fastq.gz"
```

```r
interleave_filepairs(fp) |>
  deinterleave_filepairs()
```

```
## $forward
## [1] "/mnt/storage/data/edna/dada/example/E-EnvStd-01_S186_L001_R1_001.fastq.gz"
## [2] "/mnt/storage/data/edna/dada/example/E-EvStd_S358_L001_R1_001.fastq.gz"    
## [3] "/mnt/storage/data/edna/dada/example/E-HAB179_S370_L001_R1_001.fastq.gz"   
## 
## $reverse
## [1] "/mnt/storage/data/edna/dada/example/E-EnvStd-01_S186_L001_R2_001.fastq.gz"
## [2] "/mnt/storage/data/edna/dada/example/E-EvStd_S358_L001_R2_001.fastq.gz"    
## [3] "/mnt/storage/data/edna/dada/example/E-HAB179_S370_L001_R2_001.fastq.gz"
```

```r
size_filepairs(fp, min_size = 6e6)
```

```
## # A tibble: 3 × 5
##   forward               reverse              forward_size reverse_size size_pass
##   <chr>                 <chr>                       <dbl>        <dbl> <lgl>    
## 1 /mnt/storage/data/ed… /mnt/storage/data/e…      9248335      9977645 TRUE     
## 2 /mnt/storage/data/ed… /mnt/storage/data/e…      5807710      6351454 FALSE    
## 3 /mnt/storage/data/ed… /mnt/storage/data/e…      6049072      6454063 TRUE
```


### [reformat](https://github.com/BioInfoTools/BBMap/blob/a9ceda047a7c918dc090de0fdbf6f924292d4a1f/sh/reformat.sh) wrapper

We provide a convenience wrapper, `bbmap_subsample()` to subsample fastq files using `reformat.sh`

### [cutadapt]() wrapper

Not implemented yet.
