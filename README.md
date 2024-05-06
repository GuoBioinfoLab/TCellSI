# TCSS
TCSS is a tool to assess the degree of eight distinct T-cell states from gene expression profiles. The major algorithm of TCSS is shown as follows:
![image](https://github.com/JingminYang/TstateScore/blob/main/TCSSalgorithm.jpg)
## Installation

You can install the development version of TCSS by:

``` r
# install.packages("devtools")
devtools::install_github("JingminYang/TstateScore")
library(TstateScore)
```

## Getting started

#sample_expression: Sample expression data frame in TPM format by log2-transformed RNA-seq data.
```
Tstate_calcScore(sample_expression)
```
#output

The output of the function is a matrix, where each row corresponds to a score name and each column represents a sample name. 
