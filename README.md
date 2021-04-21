# BronchiCellPred
To predict bronchial cell proportions in bulk RNA-seq datasets

![introduction](https://github.com/CancanQi/BronchiCellPred/blob/ca0893d545d6adcac28c81b3ff6fd9d911f31ef4/introduction.png)

## Introduction

`BronchiCellPred` is an r package containing models for predicting the proportions of bronchial epithelial cells using bulk gene expression data from bronchial biopsies. <br>
Gene signature matrix used for cell type deconvolution were generated using a single-cell RNAseq dataset of bronchial biopsies.<br>
Three cell type deconvolution methods were available in this package, including nnls, MuSiC (https://github.com/xuranw/MuSiC) and BSeq-sc(https://github.com/shenorrLab/bseqsc). Users can choose the methods according to their own preference.  <br>
For additional details on methods and results please go our [manuscript](link to be updated).

## Installation

library(devtools)
install_github("CancanQi/BronchiCellPred") (to be updated)

## Quick start

#### Load the library and dependencies.

```R
library(MuSiC)
library(bseqsc)
library(xbioc)
library(reshape2)
library(tidyverse)
library(Biobase)
```

#### Prepare input data

```R
## read bulk gene expression count table, with genes in row and samples in column
count.table<-read.table("your file name")
## convert count table to expression set
bulk.eset <- Biobase::ExpressionSet(assayData = data.matrix(count.table))
```

#### Predict bronchial cell proportions

```R
## use the bulk expression set as the input, chose the methods from nnls, MuSiC and bseq
est_prop<-BronCell.prop(bulk.eset,method="MuSiC")
## the output is a dataframe of cell proporion for each cell type (column) and each sample (row)

## plot the results
BronCell.plot(est_prop)

```
![result](https://github.com/CancanQi/BronchiCellPred/blob/ca0893d545d6adcac28c81b3ff6fd9d911f31ef4/main_result.png)


