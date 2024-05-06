# TCSS
TCSS is a tool to assess the degree of eight distinct T-cell states including Quiescence, Regulating, Proliferation, Helper, Cytotoxicity, Progenitor exhaustion, Terminal exhaustion, and Senescence based on the degree of resting, activation, and suppression using specific gene sets and a compiled reference spectrum from transcriptomic data. The major algorithm of TCSS is shown as follows:

![image](https://github.com/JingminYang/TstateScore/blob/main/TCSSalgorithm.jpg)
## Installation

You can install the development version of TCSS by:

``` r
# install.packages("devtools")
devtools::install_github("VyvyanYjm/TCSS")

library(TCSS)
```

## Example tutorial

#sample_expression: Sample expression data frame in TPM format by log2-transformed RNA-seq data.
```
data(exampleSample)
ResultScores <- Tstate_calcScore(exampleSample)
``` 
#output

The output of the function is a matrix, where each row corresponds to a score name and each column represents a sample name. 
