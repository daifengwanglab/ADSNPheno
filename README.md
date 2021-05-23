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

## Software Requirements!



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




## Support

Please note that this code was developed by Saniya Khullar and Daifeng Wang, Ph.D.
If you experience any issues with the code or steps in the pipeline, please reach out to Saniya Khullar at skhullar2@wisc.edu and she will be very happy to help you!

YouTube tutorials will also be available on ADSNPheno on Saniya's YouTube channel to guide you every step of the way:  https://www.youtube.com/channel/UCNhVAcIdarXzTCWZ27N1EmQ.  We also have presented our computational pipeline at an ISCB-SC webinar https://www.youtube.com/watch?v=ITwEzqhQnZU. 


   
## Demo

This demo applies our general, open-source computational pipeline, **ADSNPheno**, in the context of Alzheimer's disease (AD).  Here, we aim to reveal underlying gene regulatory mechanisms of AD risk variants to different AD phenotypes. Particularly, ADSNPheno first identifies the gene co-expression modules for various AD phenotypes via clustering gene co-expression networks. ADSNPheno further predicts the transcription factors (TFs) that significantly regulate the genes in each module, as well as the AD SNPs interrupting the TF binding sites on the regulatory elements. Finally,  ADSNPheno constructs a full gene regulatory network linking SNPs, TFs, and regulatory elements to target genes for each phenotype. This network thus provides mechanistic insights of gene regulation from disease risk variants to the phenotype in AD.

Below, please note that we provide an example application of ADNSPheno to the Hippocampus Ca1 gene expression data set GSE1297: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE1297. 

**Step 1: Disease gene expression data at the population level**

Since this is microarray data, a lot of pre-processing was needed. 

We followed "RMA normalization for microarray data" (https://felixfan.github.io/RMA-Normalization-Microarray/) for the Demo (that was specifically tailored for our same GSE1297 dataset). 

Ultimately, after pre-processing, our final dataset for the Hippocampus is:


