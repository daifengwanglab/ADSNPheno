README
================

Integrative analysis of multi-omics reveals gene regulatory networks across brain regions from risk variants to phenotypes of Alzheimer's disease and Covid-19
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
This integrative analysis, named as ADSNPheno, aims to predict gene regulatory networks of AD risk Single-Nucleotide Polymorphisms (SNPs) to different AD phenotypes. In particular, we first cluster gene co-expression networks and
identifies the gene modules for various AD phenotypes. Next, we further
predict the transcription factors (TFs) that significantly regulate the
genes in each module, as well as the AD SNPs interrupting the TF binding
sites on the regulatory elements. Finally, we construct a full
gene regulatory network linking SNPs, interrupted TFs, and regulatory
elements to target genes for each phenotype. This network thus provides
mechanistic insights of gene regulation from disease risk variants to AD
phenotypes.

Please note these steps in our integrative analysis (that can be summarized as a pipeline):

<p align="center">
  <img width="1500" src="adsnpheno.png">
</p>


## Hardware Requirements

Please note that this integrative analysis is based on R 4.0 and Python. You will only need a standard computer with enough RAM to support the operations. For predicting gene regulatory networks, a *Linux* system with 32 GB RAM and 32GB storage would be enough to support.

## Software Requirements

Please note that you will need to have R installed as well as Python (Version 3.8.5 or above).  In Python, we will be installing packages such as numpy, pandas, math, datetime, sklearn, and fsspec. 
Please note that you can find more information on R packages that are needed under [packages.R](https://github.com/daifengwanglab/ADSNPheno/blob/master/code/setupScripts/packagesNeeded.R).  Please note that many of the R packages do require Bioconductor. 

Please install the following main R packages prior to using ADSNPheno from an R terminal:
```{r}
c('plyr', 'dplyr', 'doRNG', 'reshape2', 'stringr','reticulate','hash', 'lubridate',
  'rvest', 'ggplot2', 'glmnet', 'data.table','doParallel', 'foreach', 'Seurat', 'Rmagic')

```
We also need to install BiocManager so we can download packages from Bioconductor to run in R.  Please run this command to install BiocManager if this is your first time using it: 
```{r}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
```

Then, please install these packages:
```{r}
BiocManager::install(c('JASPAR2018','S4Vectors', 'biomaRt','SNPlocs.Hsapiens.dbSNP142.GRCh37','GenomicRanges','TFBSTools','motifmatchr','GenomicFeatures','GenomeInfoDb','IRanges', 'aliases2entrez', 'biomaRt', 'GENIE3', 'GenomicInteractions', 'motifbreakR', 'rentrez', 'trena', 'WGCNA', 'clusterProfiler'
'DOSE', 'meshes', 'MotifDb', 'msigdbr', 'RedeR', 'RTN', 'rvest', 'DESeq2',
'TxDb.Hsapiens.UCSC.hg19.knownGene', 'MeSH.Hsa.eg.db', 'org.Hs.eg.db', 'BSgenome.Hsapiens.UCSC.hg19', 'SNPlocs.Hsapiens.dbSNP144.GRCh37'))
```

We will be calling a Python script [BuildingFullAndSubNetworksPythonCode.py](https://github.com/daifengwanglab/ADSNPheno/blob/master/code/setupScripts/BuildingFullAndSubNetworksPythonCode.py) from R itself, using the [reticulate](https://rstudio.github.io/reticulate/) package.  Please note we used Python version 3.8.5 for our code.  In [packages.R](https://github.com/daifengwanglab/ADSNPheno/blob/master/code/setupScripts/packagesNeeded.R), we create a [virtual environment](https://docs.python.org/3/tutorial/venv.html) (that we call "r-reticulate") where we download these 5 packages in python version 3.8.5:  [scipy](https://www.scipy.org/), [pandas](https://pandas.pydata.org/), [numpy](https://numpy.org/), [datetime](https://docs.python.org/3/library/datetime.html) and [fsspec](https://filesystem-spec.readthedocs.io/en/latest/).  

For our Covid-19 Severity Machine Learning task, we will use ipython notebooks like [AD-Covid GRN for LTL.ipynb](https://github.com/daifengwanglab/ADSNPheno/blob/master/code/otherScripts/Covid19_Analysis/Part2_MachineLearningPredictionOfCovid19SeverityFromADCovidGRNs/AD-Covid%20GRN%20for%20LTL.ipynb) and you will need to have sklearn, scipy, matplotlib, and seaborn (optional, but for visualizing confusion matrix) [installed on your computer](https://packaging.python.org/tutorials/installing-packages/) by typing in each of the following commands on terminal or command line:

```{r}
pip install sklearn
pip install scipy
pip install matplotlib
pip install seaborn
```

## Support

Please note that this code was developed by **Saniya Khullar** and **Daifeng Wang, Ph.D**.
If you experience any issues with the code or steps, please reach out to Daifeng Wang, Ph.D. at daifeng.wang@wisc.edu for more assistance.

YouTube tutorials will also be available on ADSNPheno on [Saniya's YouTube channel](https://www.youtube.com/channel/UCNhVAcIdarXzTCWZ27N1EmQ) to guide you every step of the way. We also have presented our computational pipeline at an [ISCB-SC webinar](https://www.youtube.com/watch?v=ITwEzqhQnZU). 

Please visit our [lab website](https://daifengwanglab.org/) at the University of Wisconsin - Madison, to learn more about our team and our work! 
   
## Demo
This demo applies our integrative analysis, ADSNPheno, in the context of Alzheimer's disease (AD).  Here, we aim to reveal underlying gene regulatory mechanisms of AD risk variants to different AD phenotypes. Particularly, ADSNPheno first identifies the gene co-expression modules for various AD phenotypes via clustering gene co-expression networks. ADSNPheno further predicts the transcription factors (TFs) that significantly regulate the genes in each module, as well as the AD SNPs interrupting the TF binding sites on the regulatory elements. Finally, ADSNPheno constructs a full gene regulatory network linking SNPs, TFs, and regulatory elements to target genes for each phenotype. This network thus provides mechanistic insights of gene regulation from disease risk variants to the phenotype in AD.


We provide a demo of our integrative analysis to the Lateral Temporal Lobe (LTL) region for predicting its gene regulatory network and linking SNPs to AD phenotypes. The gene expression data set and Enhancers data set are from [GSE159699](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE159699). For the demo, we randomly selected 200 genes from the gene expression data, and around 1,222 rows of Enhancer epigenomics data for the Lateral Temporal Lobe. Please see all codes for running this demo at [inputsNeededByUser_DemoLTL.R](https://github.com/daifengwanglab/ADSNPheno/blob/master/code/userInputs/inputsNeededByUser_DemoLTL.R). 

## Code for General Usage

In addition, all codes for our integrative analysis for general usage is provided in a single file [adsnphenoCodesToRun.R](https://github.com/daifengwanglab/ADSNPheno/blob/master/code/adsnphenoCodesToRun.R). Also, we provide the codes for each step of our analysis as follows (where you can edit the parameters directly within the file itself, at the top, and then run the respective code below in that file):

*  [Identifying gene co-expression modules and modular enrichments](https://github.com/daifengwanglab/ADSNPheno/tree/master/code/otherScripts/Part1_Identifying%20Gene%20Co-Expression%20Modules%20and%20Module%20Enrichments)
*  [Predicting gene regulatory network linking TFs, regulatory elements to target genes](https://github.com/daifengwanglab/ADSNPheno/tree/master/code/otherScripts/Part2_PredictingGeneRegulatoryNetworksLinkingTFsAndRegElementsToTarget%20Genes)
*  [Linking AD SNPs to interrupted TFBSs, target genes to phenotypes](https://github.com/daifengwanglab/ADSNPheno/tree/master/code/otherScripts/Part3_LinkingSNPsToInterruptedTFBindingSitesAndDysregulatedTGsAndPhenotypes)
*  [Differential expressions analysis of Covid-19 severity](https://github.com/daifengwanglab/ADSNPheno/tree/master/code/otherScripts/Covid19_Analysis/Part1_DifferentialExpressionAnalysis)
*  [Machine learning predicting Covid-19 severity from AD-Covid gene regulatory networks](https://github.com/daifengwanglab/ADSNPheno/tree/master/code/otherScripts/Covid19_Analysis/Part2_MachineLearningPredictionOfCovid19SeverityFromADCovidGRNs) 
*  [Decision curve analysis for Covid-19 severity using AD-Covid genes and benchmarking genes](https://github.com/daifengwanglab/ADSNPheno/tree/master/code/otherScripts/Covid19_Analysis/Part3_DecisionCurveAnalysisForCovid19SeverityPredictions)


You can run the analysis for the entire Hippocampus Ca1 region (by running [inputsNeededByUser_Hippo.R](https://github.com/daifengwanglab/ADSNPheno/blob/master/code/userInputs/inputsNeededByUser_Hippo.R)) and the entire Lateral Temporal Lobe region (by running [inputsNeededByUser_LatTempLobe.R](https://github.com/daifengwanglab/ADSNPheno/blob/master/code/userInputs/inputsNeededByUser_LatTempLobe.R)). Please ensure that you adjust the parameters (especially the file paths) so that you are able to run those files on your computer/server. 


Tutorial
-------

Please note that a detailed Tutorial of this code will be available soon.
Please let us know if you have any questions!  Thank you!


```
