# TCellSI 
T cell state identifier (TCellSI) is a tool to access eight distinct T cell states including Quiescence, Regulating, Proliferation, Helper, Cytotoxicity, Progenitor exhaustion, Terminal exhaustion, and Senescence. TCellSI provides T cell state scores (TCSS) for samples using specific marker gene sets and a compiled reference spectrum of T cell states from transcriptomic data. The major algorithm of TCellSI is shown as follows: 

![image](https://github.com/VyvyanYjm/TCellSI/blob/main/algorithm.jpg)
## Installation

You can install the development version of TCellSI by:

``` r
# install.packages("devtools")
devtools::install_github("GuoBioinfoLab/TCellSI")

library(TCellSI)
```

## Example tutorial

#sample_expression: Complete gene expression data.frame in TPM format by log2-transformed RNA-seq data.

```
sample_expression <- TCellSI::exampleSample
# sample_expression
#            SRR5088825  SRR5088828  SRR5088830
# 5_8S_rRNA    0.000000    0.000000   0.0000000
# 5S_rRNA      0.000000    0.000000   0.4346853
# 7SK          0.0c0000    0.000000   0.0000000
# A1BG         1.907599    3.284418   4.0821248
# A1BG-AS1     3.083914    2.021501   4.5169002
ResultScores <- TCSS_Calculate(sample_expression) 
```
#output:
The output of the function is a data.frame with TCSS metrics, where each row corresponds to a T cell state and each column represents a sample name.  

```
# ResultScores
#                        SRR5088825  SRR5088828  SRR5088830
# Quiescence              0.6422158   0.6153140   0.5105856
# Regulating              0.4404173   0.3974391   0.1670356
# Proliferation           0.5780474   0.5761891   0.3994404
# Helper                  0.5647604   0.4063729   0.2121070
# Cytotoxicity            0.6482720   0.5826073   0.1829111
# Progenitor_exhaustion   0.5905679   0.5308304   0.2482848
# Terminal_exhaustion     0.6624943   0.6222188   0.3311532
# Senescence              0.5921543   0.5725185   0.2896508
```

#If you want to apply this method to other cell states that interest you, you should compile a reference spectrum and prepare specific marker gene sets of your cell states. You can then calculate the scores for your cell states using the following function.

```
OtherScores <- CSS_Calculate‎(sample_expression, reference = XXX, markers = XXX)
```
#Forms of reference and markers look like:

```
#reference
#The self-constructed reference should contain log2-transformed, TPM-normalized gene expression data from RNA-seq or scRNA-seq. 
#             cell_state1  cell_state2  cell_state3
# DDX11L1      0.32323232   0.54567463   0.32456323
# WASH7P       0.82670591   1.89565638   1.40492732
# MIR6859-1    0.02172025   0.03816506   0.52313432
# MIR1302-2HG  0.00000000   0.00000000   0.00032302
# MIR1302-2    0.00000000   0.00000000   0.00002132
```
```
#markers: A list of multiple cell states containing specific gene sets
#The number of marker genes per cell state can vary.
#$cell_state1
#[1] "XXX"  "XXX"  "XXX" ...
#$cell_state2
#[1] "XXX"  "XXX"  "XXX" ...
#$cell_state3
#[1] "XXX"  "XXX"  "XXX" ...
```
![image](https://github.com/VyvyanYjm/TCellSI/blob/main/Logo.jpg)
