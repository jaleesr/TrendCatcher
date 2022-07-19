# TrendCatcher

# Introduction

TrendCatcher is a versatile R package for identifying dynamic differentially expressed genes (DDEGs) in RNA-seq longitudinal studies. A time course experiment is a widely used design in the study of cellular processes such as cell differentiation or response to external stimuli. Temporal changes to the gene expression, such as mRNA, is the key to characterizing such biological processes. Here, we present a versatile R package named TrendCatcher to identify the dynamic differentially expressed genes along the biological process time course. 


# Installation

* Install latest development version from GitHub (requires [devtools](https://github.com/hadley/devtools) package):

```r
if (!require("devtools")) {
  install.packages("devtools")
}
devtools::install_github("jaleesr/TrendCatcher", dependencies = TRUE, build_vignettes = FALSE)
```

# TrendCatcher Framework Overview

TrendCatcher requires 2 main inputs: the raw count table C of a temporal study with a dimension of m × n, where m denotes the number of genes and n denotes the number of samples, and a user-defined baseline time variable T, such as “0 hour”. Since samples may have different sequencing depths and batch effect, TrendCatcher integrates with limma and provides preprocessing steps, such as batch correction and normalization. For scRNA-Seq data sets, TrendCatcher extracts cells for each cell type annotated in the meta data slot of Seurat object and converts it into a cell type–specific “pseudobulk” time course RNA library. Based on a user-specified threshold, genes of relatively low abundance are removed from the count table, reads are normalized, and batch effects are removed. TrendCatcher’s core algorithm is composed of 5 main steps: (a) baseline fluctuation confidence interval estimation, (b) model dynamic longitudinal count, (c) time point dynamic P value calculation, (d) gene-wise dynamic P value calculation, and (e) break point screening and gene-wise dynamic pattern assignment. Mathematical details will be expanded in the following sections. 

For the output of TrendCatcher, there are mainly 2 components: a master table and a set of functions for versatile visualization purposes. The master table contains all the dynamic details of each single gene, including its dynamic P value, its break point location time, and its dynamic trajectory pattern. In addition to the master table, TrendCatcher produces 5 main types of visualizations: (a) a figure showing the observed counts and fitted splines of each gene, (b) genes trajectories from each trajectory pattern group, (c) a hierarchical pie chart that represents trajectory pattern composition, (d) a TimeHeatmap to infer trajectory dynamics of top dynamic biological pathways, and (e) a 2-sided bar plot to show the top most positively and negatively changed (averaged accumulative log2FC) biological pathways.

![plot](./man/figures/TrendCatcherWorkFlow.png)

