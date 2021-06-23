README
================

ADSNPheno: A computational pipeline of multi-omics analysis to predict gene regulatory networks from disease variants to phenotypes with applications in Alzheimer’s disease and SARS-CoV-2
===========================================================================================================================================================================================
Paper
-------

Please note that this code corresponds to this current pre-print, available on Biorxiv: [Integrative analysis of multi-omics reveals gene regulatory networks across brain regions from risk variants to phenotypes of Alzheimer's disease and Covid-19](https://www.biorxiv.org/content/10.1101/2021.06.21.449165v1.full.pdf+html)
by Saniya Khullar and Daifeng Wang, Ph.D.


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

Please note that this analysis is based on R 4.0 and Python. You will only need a standard computer with enough RAM to support the operations. For predicting gene regulatory networks, a *Linux* system with 32 GB RAM and 32GB storage would be enough to support.

## Software Requirements

Please note that you will need to have R installed as well as Python (Version 3.8.5 or above).  In Python, we will be installing packages such as numpy, pandas, math, datetime, sklearn, and fsspec. 
Please note that you can find more information on R packages that are needed under code $\rightarrow$ setupScripts $\rightarrow$ packagesNeeded.R.  Please note that many of the R packages do require Bioconductor. 


## Support

Please note that this code was developed by Saniya Khullar and Daifeng Wang, Ph.D.
If you experience any issues with the code or steps in the pipeline, please reach out to Saniya Khullar at skhullar2@wisc.edu and she will be very happy to help you!

YouTube tutorials will also be available on ADSNPheno on [Saniya's YouTube channel](https://www.youtube.com/channel/UCNhVAcIdarXzTCWZ27N1EmQ) to guide you every step of the way. We also have presented our computational pipeline at an [ISCB-SC webinar](https://www.youtube.com/watch?v=ITwEzqhQnZU). 

Please visit our [lab website](https://daifengwanglab.org/) at the University of Wisconsin - Madison, to learn more about our team and our work! 
   
## Demo

This demo applies our general, open-source computational pipeline, **ADSNPheno**, in the context of Alzheimer's disease (AD).  Here, we aim to reveal underlying gene regulatory mechanisms of AD risk variants to different AD phenotypes. Particularly, ADSNPheno first identifies the gene co-expression modules for various AD phenotypes via clustering gene co-expression networks. ADSNPheno further predicts the transcription factors (TFs) that significantly regulate the genes in each module, as well as the AD SNPs interrupting the TF binding sites on the regulatory elements. Finally, ADSNPheno constructs a full gene regulatory network linking SNPs, TFs, and regulatory elements to target genes for each phenotype. This network thus provides mechanistic insights of gene regulation from disease risk variants to the phenotype in AD.

Please note that we provide an example application of ADNSPheno to the Lateral Temporal Lobe gene expression data set [GSE159699](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE159699).  For the demo, we randomly selected 200 genes from the gene expression data, and around 1,222 rows of Enhancer epigenomics data for the Lateral Temporal Lobe. 


Tutorial
-------

Please note that a detailed Tutorial of this code will be available here: [ADSNPheno Book](https://saniyakhullar.github.io/SNPheno/).

Please let us know if you have any questions!  Thank you!


```
