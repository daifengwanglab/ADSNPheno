README
================

#`{r setup, include=FALSE} knitr::opts_chunk$set(echo = TRUE)`

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

Please note that this analysis is based on R 4.0 and Python. You will only need a standard computer with enough RAM to support the operations. For predicting gene regulatory networks, a *Linux* system with 32 GB RAM and 32GB storage would be enough to support.

## Software Requirements



This is an R Markdown format used for publishing markdown documents to
GitHub. When you click the **Knit** button all R code chunks are run and
a markdown file (.md) suitable for publishing to GitHub is generated.



## Support

Please note that this code was developed by Saniya Khullar and Daifeng Wang, Ph.D.
If you experience any issues with the code or steps in the pipeline, please reach out to Saniya Khullar at skhullar2@wisc.edu and she will be very happy to help you!

YouTube tutorials will also be available on ADSNPheno on [Saniya's YouTube channel](https://www.youtube.com/channel/UCNhVAcIdarXzTCWZ27N1EmQ) to guide you every step of the way. We also have presented our computational pipeline at an [ISCB-SC webinar](https://www.youtube.com/watch?v=ITwEzqhQnZU). 

Please visit our [lab website](https://daifengwanglab.org/) at the University of Wisconsin - Madison, to learn more about our team and our work! 
   
## Demo

This demo applies our general, open-source computational pipeline, **ADSNPheno**, in the context of Alzheimer's disease (AD).  Here, we aim to reveal underlying gene regulatory mechanisms of AD risk variants to different AD phenotypes. Particularly, ADSNPheno first identifies the gene co-expression modules for various AD phenotypes via clustering gene co-expression networks. ADSNPheno further predicts the transcription factors (TFs) that significantly regulate the genes in each module, as well as the AD SNPs interrupting the TF binding sites on the regulatory elements. Finally, ADSNPheno constructs a full gene regulatory network linking SNPs, TFs, and regulatory elements to target genes for each phenotype. This network thus provides mechanistic insights of gene regulation from disease risk variants to the phenotype in AD.

Below, please note that we provide an example application of ADNSPheno to the Hippocampus Ca1 gene expression data set [GSE1297](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE1297). 



**Step 0: Please setup the appropriate workspace environment**

Please note that it is important to download our ADSNPheno code onto your local computer or server, where you can access it. It is recommended that you work with the layout of our existing folders (with our data files and code), and then make adjustments as needed, such as:

* adding your respective data files to the appropriate folders (please go to data folder)

* changing the parameters in the code based on your needs and requirements (by going to code > userInputForEachStep)

* adapting our code to further suit your specific needs

This is to help you adapt on our code as well, better understand the tools involved, and be more user-friendly.

When you download our pipeline ADSNPheno, as a zip file("ADSNPheno-master.zip"), it hopefully will come with the sub-folders shown above (such as the "code" and "data" folders). Please unzip the folder and then note the complete file path location of the unzipped folder, which will be the variable adsnphenoDirectory. For instance, on Saniya's computer, this ADSNPheno folder is located here: "F://organizedAlzheimers//ADSNPheno//", so adsnphenoDirectory = "F://organizedAlzheimers//ADSNPheno//".

We also have shared 3 data sources that we have used for our analysis, which we believe may be general-purpose and helpful references for you.  These are also found in the dataDefaults folder (data > dataDefaults) and are information on:

*  1. Additional Gene Regulatory Information:  Transcription Factor (TF) to Target Gene (TG) relationships based on TRRUST 2.0, which may serve as an additional source of TF-TG interactions.  Information is also provided on some relationships and if the TF activates or represses the TG. 



Please note an example of some rows of our "data//dataDefaults//additionalGeneRegulatoryInfo//Trrust2.csv". Here, CombinedName basically uses the gene symbol names and 2 pipes to represent the relationship:  "TF || RegulatedGene".  Thus, the 1st row has "AATF"" as the TF and "CDKN1A"" as the RegulatedGene, so CombinedName will be "AATF || CDKN1A"
```{r}
dim(trrust2DF) # 9,395 rows of unique TF-TG relationships and these 5 columns (at least) are needed.  

head(trrust2DF) # the first 5 rows of the TRRUST2.0 data source.
#     TF RegulatedGene  Source                             Info   CombinedName
# 1 AATF        CDKN1A Trrust2    TRRUST2 (Regulation: Unknown) AATF || CDKN1A
# 2 AATF          KLK3 Trrust2    TRRUST2 (Regulation: Unknown)   AATF || KLK3
# 3 AATF           MYC Trrust2 TRRUST2 (Regulation: Activation)    AATF || MYC
# 4 AATF          TP53 Trrust2    TRRUST2 (Regulation: Unknown)   AATF || TP53
# 5 ABL1           BAX Trrust2 TRRUST2 (Regulation: Activation)    ABL1 || BAX
# 6 ABL1          BCL2 Trrust2 TRRUST2 (Regulation: Repression)   ABL1 || BCL2
```

*  2. Entrez Mapping Information:  To assist you with your gene set, please note that we have tried to build an extensive resource on some of the Entrez IDs for various genes so that you can map many genes to Entrez IDs. Please beware that Excel may distort some gene names like SEPT7 and MARC1 to numeric date values by mistake (we have tried to troubleshoot for some of those situations).  It is not exhaustive, so please check additional sources such as NCBI for any genes you cannot find and make adjustments to this default data source, as needed, if you find any discrepancies. And, please email Saniya at skhullar2@wisc.edu, so she can update the files as needed (thanks for your help!). For many genes that begin with "LOC", the entrez ID is often the numeric part of the gene name.   This mapping will be used as a lookup in the pipeline to help with merging relationships, from disparate sources, reliably!

Please note an example of some rows of our "data//dataDefaults//entrezMappingInfo//infoDFMappingsGeneSymbolAndID.csv":
```{r}
someRows = c(1:5, 200:205)# 10 different rows of the CSV as an example
dim(infoDF) # 59,307 rows and 2 columns

infoDF[someRows,]  # geneSymbolName mapping to respective entrezIDs. Please double-check.
#     geneSymbolName  entrezID
# 1           15.212      3805
# 2            47.11      3802
# 3             3635    201626
# 4           Jul-60     11054
# 5           19-Feb     60343
# 200        ABHD14B     84836
# 201         ABHD15    116236
# 202     ABHD15-AS1 104355133
# 203        ABHD16A      7920
# 204        ABHD16B    140701
# 205        ABHD17A     81926
```

* 3.  Transcription Factor Information:  We shared a combined list of Transcription Factors by Lambert and Jasper (as well as EntrezIDs for several of them), and had 1,668 unique TFs.  

```{r}
# dim(tfsDF) # 1,668 rows (TFs) and 2 columns: TFs, entrezID

# head(tfsDF) # the 1st 5 rows
#   TFs entrezID
# 1  RUNX1      861
# 2 TFAP2A     7020
# 3  NR2F1     7025
# 4  CREB1     1385
# 5   E2F1     1869
# 6  NFIL3     4783

```


To adjust data sources, please note that you can just go to the "data > dataDefaults" folder. 

Please navigate to code > userInputForEachStep and then go to "step0_loadingInDefaultDatasets.R" to make adjustments to any of these parameters below.  

Below, please note the main parameters to adjust (and please see step0_loadingInDefaultDatasets.R for additional parameters that may be varied) 

| Variable Input Needed by You | Description | Default Value |
|------|------|------|
| adsnphenoDirectory | The region of the body that your gene expression data is for. In our example, it is the Hippocampus Ca1 Region | adsnphenoDirectory = "F://organizedAlzheimers//ADSNPheno//"
| tfsUsed | The sources of information for the Transcription Factors (TFs).  Here, we use TFs found  in Lambert and/or Jaspar as our default. | tfsUsed = "LambertAndJaspar" |
| transcriptionFactorsFilePath | In the dataDefaults > transcriptionFactorInfo directory, this is the name of the CSV file with the list of TFs you want to use. Here, our list is for Lambert and/or Jaspar TFs. | transcriptionFactorsFilePath = "LambertAndJasparTFs.csv" |
| geneToEntrezIDMappingPath | In the dataDefaults > entrezMappingInfo directory, this is the name of the CSV file with 2 columns:  [Symbols, entrezID].  Here, "Symbols" is the geneName (or # or numeric date), and "entrezID" is what we typically found as the entrezID.   | geneToEntrezIDMappingPath = "infoDFMappingsGeneSymbolAndID.csv" | |
| additionalGeneRegulatoryNetworkFilePath |  In the dataDefaults > additionalGeneRegulatoryInfo directory, this is the name of the CSV file with external gold standard gene regulatory network information you would like to include in the final analysis. Here, we used TRRUST2.0.  Our data file needs to contain at least these 5 columns: [TF, RegulatedGene, Source, Info, CombinedName].  CombinedName is basically a column based on the gene symbols and is separated with 2 pipes | additionalGeneRegulatoryNetworkFilePath = "Trrust2.csv" | |
|------|------|------|
|numberOfRoundingDigits| Please note that this is the # of digits after the decimal place for rounding any calculations performed during ADSNPheno. | numberOfRoundingDigits = 3|
|includeTimeStamp| Please note that if this Boolean variable is set to TRUE, then additional time-stamped data files will be saved as the pipeline is run. | includeTimeStamp = TRUE|
|disease| Please note that this is where you would mention the disease or condition or relevant biological context.  Here, we are analyzing the potential role of non-coding SNPs in Alzheimer's disease, so we set disease = "Alzheimers". | disease = "Alzheimers"|
|outputPath| Please note that this is where the outputs of ADSNPheno would be stored. By default, ADSNPheno will build an output folder and store the corresponding datafiles and results in sub-folders there.| outputPath = "ADSNPheno//ADSNPhenoOutputs"|


**Step 1: Disease gene expression data at the population level**

Since this is microarray data, a lot of pre-processing was needed. The dataset had total RNA expression values for 22,283 HG-U133 Affymetrix Human Genome U133 Plus 2.0 Microarray Identifier probes for 31 individual postmortem samples.  These individual samples include 9 control samples (no AD), 7 samples at the initial stage, 8 samples at the moderate stage, and 7 samples at the severe stage. 

Please note that we utilized Bioconductor packages, such as: GEOquery, hgu133a.db, hgu133acdf, and Affy R packages to download the raw data. Then, we performed Robust Multichip Average (RMA) normalization [31] to account for background and technical variations among these samples.  In fact, we followed [RMA normalization for microarray data](https://felixfan.github.io/RMA-Normalization-Microarray/) for the Demo (that was specifically tailored for our same GSE1297 dataset). 


Next, we mapped those probes of microarray data to genes, averaging values that mapped to the same gene Entrez ID.  We removed probes that did not map to any known genes.  Please note that some helpful tools for obtaining the status of certain genes (such as updated gene names or status on whether they are withdrawn or discontinued) is: https://www.genenames.org/tools/multi-symbol-checker/.  Some gene names can map to the same Entrez ID.  In addition, Entrez IDs can be more robust (and thereby better when using different data tools with various gene name databases), especially for some gene-set enrichment sources.


In addition, we consulted the aliases2entrez package to retrieve EntrezID information on the genes.  To assist you with your gene set, please note that we have tried to build an extensive resource on some of the Entrez IDs for various genes, which we have stored here: https://github.com/daifengwanglab/ADSNPheno/blob/master/data/entrezMappingInfo/infoDFMappingsGeneSymbolAndID.csv  [infoDFMappingsGeneSymbolAndID.csv]("https://github.com/daifengwanglab/ADSNPheno/blob/master/data/entrezMappingInfo/infoDFMappingsGeneSymbolAndID.csv") (Please see: data > entrezMappingInfo > infoDFMappingsGeneSymbolAndID.csv). 


Ultimately, after pre-processing, our final dataset for the Hippocampus had 13,073 unique genes for these 31 samples. (Please see: data > originalGeneExpressionDataSets > originalHippocampalCa1RegionGeneExpressionData.csv). 

| Variable Input Needed by You | Description | Demo Value |
|------|------|------|
| bodyRegion | The region of the body that your gene expression data is for. In our example, it is the Hippocampus Ca1 Region | bodyRegion = "HippocampusCa1"
| gene | Gene name | |
| weight | The weight of the variant to gene expression with respect to effect allele counts | |
| ref_allele | Reference allele | |
| eff_allele | Effect allele | |
| bd | The gene coordinate | |


```{r}
inputGeneExpressionDataFilePath = "ADSNPheno\\data\\originalGeneExpressionDataSets\\originalHippocampalCa1RegionGeneExpressionData.csv"




```
