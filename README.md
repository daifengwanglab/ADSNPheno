README
================

`{r setup, include=FALSE} knitr::opts_chunk$set(echo = TRUE)`

ADSNPheno: A computational pipeline of multi-omics analysis to predict gene regulatory networks from disease variants to phenotypes with applications in Alzheimer’s disease and SARS-CoV-2
===========================================================================================================================================================================================

Summary
-------

Genome-wide association studies (GWAS) have found many genetic risk
variants associated with Alzheimer’s disease (AD). However, how these
risk variants affect deeper phenotypes such as disease progression and
immune response remains elusive. Also, our understanding of cellular and
molecular mechanisms from disease risk variants to various disease
phenotypes is still limited. To address these problems, we developed a
computational pipeline of integrated multi-omics analysis from genotype,
transcriptomics, epigenomics to phenotypes for revealing gene regulatory
mechanisms from disease variants to phenotypes.

Method
------

This pipeline, ADSNPheno, aims to predict gene regulatory networks of AD
risk Single-Nucleotide Polymorphisms (SNPs) to different AD phenotypes.
In particular, ADSNPheno first clusters gene co-expression networks and
identifies the gene modules for various AD phenotypes. ADSNPheno further
predicts the transcription factors (TFs) that significantly regulate the
genes in each module, as well as the AD SNPs interrupting the TF binding
sites on the regulatory elements. Finally, ADSNPheno constructs a full
gene regulatory network linking SNPs, interrupted TFs, and regulatory
elements to target genes for each phenotype. This network thus provides
mechanistic insights of gene regulation from disease risk variants to AD
phenotypes.

Please note our pipeline, ADSNPheno:

<p align="center">
  <img width="1500" src="adsnpheno.png">
</p>


## Hardware Requirements

Please note that this analysis is based on R 4.0. You will only need a standard computer with enough RAM to support the operations. For predicting gene regulatory networks, a *Linux* system with 32 GB RAM and 32GB storage would be enough to support.

## Software Requirements



This is an R Markdown format used for publishing markdown documents to
GitHub. When you click the **Knit** button all R code chunks are run and
a markdown file (.md) suitable for publishing to GitHub is generated.


Including Code
--------------

You can include R code in the document as follows:

`{r cars} summary(cars)`

Including Plots
---------------

You can also embed plots, for example:

`{r pressure, echo=FALSE} plot(pressure)`

Note that the `echo = FALSE` parameter was added to the code chunk to
prevent printing of the R code that generated the plot.
