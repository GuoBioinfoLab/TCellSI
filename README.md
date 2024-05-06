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
ResultScores <- Tstate_calcScore(exampleSample)
``` 
#output:
The output of the function is a matrix, where each row corresponds to a score name and each column represents a sample name.  
```
# ResultScores
#                      ERR3502705 ERR3502706 ERR3502712
#Quiescence             0.7721079  0.7722078  0.7482268
#Terminal_exhaustion    0.6513540  0.6003654  0.6041851
#Progenitor_exhaustion  0.6793562  0.6450539  0.6275441
#Senescence             0.6846261  0.6487488  0.6340548
#Cytotoxicity           0.6529064  0.5928662  0.6044746
#Regulating             0.5308267  0.4070450  0.4134950
#Proliferation          0.5340890  0.4552266  0.4388102
#Helper                 0.6267208  0.5927787  0.5699727
```
 

