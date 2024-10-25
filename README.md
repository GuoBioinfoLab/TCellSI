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

Sample_expression: Complete gene expression data.frame in TPM format by log2-transformed RNA-seq data.

```
sample_expression <- TCellSI::exampleSample
# sample_expression
#            SRR5088825  SRR5088828  SRR5088830  ...
# 5_8S_rRNA    0.000000    0.000000   0.0000000
# 5S_rRNA      0.000000    0.000000   0.4346853
# 7SK          0.0c0000    0.000000   0.0000000
# A1BG         1.907599    3.284418   4.0821248
# A1BG-AS1     3.083914    2.021501   4.5169002
# ...
ResultScores <- TCellSI::TCSS_Calculate(sample_expression, ref = TRUE) (Optionally, ref = FALSE)
```
Output:
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

If you want to apply this method to other states that interest you, you should compile a reference spectrum and prepare specific marker gene sets of your cell states. You can then calculate the scores for your cell states using the following function. If you choose not to provide a reference expression spectrum but to do the calculation directly, you can use the parameter ref=FALSE and do not need to provide the reference parameter.

```
OtherScores <- TCellSI::CSS_Calculate‎(sample_expression, ref=TRUE, reference = XXX, markers = XXX)
```

Forms of reference and markers look like:

```
#reference
#The self-constructed reference should contain log2-transformed, TPM-normalized gene expression data from RNA-seq or scRNA-seq. 
#             cell_state1  cell_state2  cell_state3  ...
# DDX11L1      0.32323232   0.54567463   0.32456323
# WASH7P       0.82670591   1.89565638   1.40492732
# MIR6859-1    0.02172025   0.03816506   0.52313432
# MIR1302-2HG  0.00000000   0.00000000   0.00032302
# MIR1302-2    0.00000000   0.00000000   0.00002132
#              ...
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
## How to use TCellSI for scRNA-seq data?
TCellSI still shows excellent results in the calculation of single-cell data and can assist in `single-cell annotation`.

In terms of operation. First, you should extract the count expression of single-cell data by reading the count file of single-cell data directly or seurat_obj@assays$RNA@counts in the seurat object. You can use the count data directly to do the calculations of TCellSI. If the results are unsatisfactory, we suggested to convert it to the log(TPM +1) format or use other scaled method. Then you can use TCellSI to perform calculations of the states scores for each cell of the single-cell data. 
```
scRNA_scores <- TCellSI::TCSS_scRNAseqCalculate(sample_scRNA, core= XXX, ref = TRUE) # core: default value is 4; Optionally, ref = FALSE
```
Then you can add the score value of the result of the calculation into the metadata data box of the seurat object.
```
FeaturePlot(object = seurat_object, features = "TCSS")  #viewing the distribution of scores in a umap; "TCSS" is the name of the column in which the categorical value is added to the metadata object
```
In addition, if you have an single-cell population annotation, you can create pseudobulk samples and then calculate the state scores for each samples, which can reduce the problem of drop-out in the single-cell data that leads to less accurate results. The creation of the pseudobulk is as follows:
## Pseudobulk creation tutorial for single-cell data analysis
How to create pseudobulk samples from single cell data ? If you want to do this, you should prepare an expression data, which should be either log2(TPM+1) or normalized single-cell data. In this data, each row represents a gene and each column represents a cell ID (see example as follows). Also, you should prepare a single-cell annotation file, which includes columns of cell annotation and cell ID in expression file (see example as follows). 
```
# expression data
#             NP710.20180123  NP711.20180123  NP71.20180123 ...
#A1BG         0.070079488      0.216835131     6.269805313
#NAT2         0.002001509      0.003654851     0.003190016
#ADA          0.085464008      0.088970085     0.057264107
#...
# single-cell annotation file
#             UniqueCell_ID   annotation
#             NTH5.20180123   CD4_C01_CCR7
#             NTH64.20180123  CD4_C01_CCR7
#             NTR57.20180123  CD4_C01_CCR7 
#             ...
```
Then, you can use the following function to get pseudobulk samples.
```
pseudo_bulk <- TCellSI::create_pseudo_bulk(
  annotation_data = XXX, 
  expression_data = XXX, 
  cluster_col = "annotation", # the column names of annotation in single-cell annotation file
  cell_id_col = "UniqueCell_ID", # the column names of Cell_ID in single-cell annotation file
  n_clusters = 18, # number of cell types annotated
  factor = 5, # number of samples for downsampling, default is 5
  sampling_rate = 0.6 # percentage of cells downsampled, default is 0.6
) 
# see examples of the result
# pseudo_bulk, each column represents a newly pseudobulk samples, each row represents a gene.
#         CD4_C01_CCR7_bulk   CD4_C01_CCR7_bulk.1  CD4_C01_CCR7_bulk.2
#A1BG      0.495165739          0.67542360           0.737122107
#NAT2      0.006033183          0.00337272           0.007104438
#ADA       0.855647562          1.06058830           0.898952625
Result <- TCSS_Calculate(pseudo_bulk)
```
## Credit
Please cite our paper if TCellSI is helpful.

* Yang, Jing‐Min, NanZhang, Tao Luo, Mei Yang, Wen‐Kang Shen,Zhen‐Lin Tan, Yun Xia, et al. 2024. "TCellSI: A Novel Method for T Cell State Assessment and Its Applications in Immune Environment Prediction."iMeta e231. https://doi.org/10.1002/imt2.231


![image](https://github.com/VyvyanYjm/TCellSI/blob/main/Logo.jpg)
