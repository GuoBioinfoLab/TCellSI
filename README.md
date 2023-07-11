# TstateScore
TstateScore is a tool to calculate the scores of 8 T-cell status from gene expression datasets in RNA-Seq

## Installation

You can install the development version of TstateScore by:

``` r
# install.packages("devtools")
devtools::install_github("JingminYang/TstateScore")
library(TstateScore)
```

## Getting started

#sample_expression: Sample expression data frame in TPM format by log2-transformed RNA-seq data.

Tstate_calcScore(sample_expression)
```
#output
```
The output of the function is a matrix, where each row corresponds to a score name and each column represents a sample name. 
