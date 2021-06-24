# :) Please note that this R file by Saniya has all the parameters we need
# Step 1:  Populate gene expression data for Hippocampus
#  how do we organize that input data


#### STEP 0: DATA DEFAULTS:
#  BY: SANIYA KHULLAR
# Please note that we have loaded in our defaults for:
# entrezMappingInfo (gene name to Entrez ID mapping) in infoDF
# Transcription Factors in TFs
# Gene Regulatory Network information from trrust2DF

# step 0:
# please set the defaults.

adsnphenoDirectory = "F://organizedAlzheimers//ADSNPheno//"  # please supply full directory path
setwd(adsnphenoDirectory) # please set this folder as working directory for code and other resources
codeFolder = paste0(adsnphenoDirectory, "code//") # the main folder where code would be within ADSNPheno

dataFolder = paste0(adsnphenoDirectory, "data//") # the main folder where data would be within ADSNPheno
# We use a combined list by Lambert and Jaspar for our TFs
transcriptionFactorsFilePath = "LambertAndJasparTFs.csv"

tfsUsed = "LambertAndJaspar"

# please note the location of the Gene Name to Entrez ID mapping
geneToEntrezIDMappingPath = paste0(dataFolder, "dataDefaults//entrezMappingInfo//infoDFMappingsGeneSymbolAndID.csv")

# our background gene regulatory data came from TRRUST2, so we input:
# please note that this data source should have columns like this: "TF", "RegulatedGene", "Source", "Info", "CombinedName"
additionalGeneRegulatoryNetworkFilePath = "Trrust2.csv"

tfSource = "TRRUST2"

numberOfRoundingDigits = 3
includeTimeStamp = "TRUE"

outputPathNameADSNPhenoOutputs = paste0("F://organizedAlzheimers//ADSNPhenoOutputs_NEW//")




getwd() # please note this is our current working directory
setwd(codeFolder) # please set this folder as working directory


disease = "Alzheimers"
diseasePetName = "AD" # or diseaseNickName


pathForPythonCode = paste0(codeFolder, "setupScripts//BuildingFullAndSubNetworksPythonCode.py")


#################################################################################################

tissueName = "Brain"
bodyRegion = "Mini_LTL_Demo" 

inputGeneExpressionDataFilePath = paste0(dataFolder, "originalGeneExpressionDataSets//LateralTemporalLobeRegionGeneExpressionData_miniDemo_200genes.csv")

computerType = "Windows" # please note that this may be helpful given that some code may work differently depending on Operating System.  ex. a Mac may differ from PC

useTFsFoundByADSNPheno = TRUE # please note that by default we use  combined list of Lambert and Jaspar TFs.
filePathOfTFsToUse = NULL # please note that if this is NULL, we will use combined list of Lambert and Jaspar TFs.  Else, please specify csv filepath of TFs, where the TFs are in 1 column with the title called "tfName", with names that are found exactly in gene expression dataset

tfsUsed = "LambertAndJaspar" # please adjust if ye are using another dataset

log2transformInputData = "TRUE" # (should we apply a log2(x+1) transform on data
scaleInputData = "FALSE"  # should we apply a scale() on data
phenotypesFilePath = "F://organizedAlzheimers//Setup//InputDataFiles//phenotypes//AlzheimersLateralTemporalLobePhenotypesUpdated.csv"


########################################################################################################################
### Weighted Gene Co-Expression Network Analysis (WGCNA) with K-Means:
# Step 2:  How to construct co-expression network and modules (WGCNA and K-means too).  Define parameters

saveOutputs = "TRUE"
netType = "signed"  # # network type:but please check WGCNA for additional parameters.  We retained the defaults

# Soft-Thresholding Power to be used:
maxNumPowers = 50 # when selecting a Soft-Thresholding Power, what is the maximum Power, beta, that you want to consider?

ourPower = NULL # if we have a recommended power in mind please put that there (this can come from analyzing the Soft-Thresholding Analysis results)
useRecommendedPower = "TRUE" #FALSE" # otherwise, please use the default that WGCNA selects


# module trait and gene-trait associations:
positiveCorrelationsOnly = TRUE #   for module-trait associations and gene-trait associations
pValueCutOffForSignificance = 0.05 # please note that associations with a p-value < pValueCutOffForSignificance will be considered statistically significant
# co-expression network part:
percentilesForCoexpressionNetworkVec = c(0, 0.25, 0.50, 0.75, 0.85, 0.90, 0.95, 1)
minTOMThresholdValue = 0.0001

# when WGCNA is performing initial step of creating modules from the TOM, what is the minimum module size?
minModuleSize = 30;


# tree cut height ??  we kept it at 0, but individuals could adjust as needed to help merge modules?
treeCutHeightForMergingCloseModules = 0  # we usually do not merge close modules

# K-Means step performed after WGCNA
performAdditionalKMeansStep = TRUE # or can be FALSE if you want original WGCNA results

# please note source of this code:
# https://github.com/juanbot/km2gcn
min.genes.for.grey=0 # originally, grey module corresponds a group of unclustered genes.  How many minimum genes do we want in grey module after k-means?
# stopping conditions
nIterations = 201 # what is the maximum # of iterations of k-means that we are willing to perform?
meg = 0 # Minimum number of genes exchanged: if before the number of iterations, the number of genes exchanged in the current iteration is less than this parameter, k-means algorithm stops.



#############################################################################################################################
####### GENE MODULE ANNOTATION:
# Step 3:  Annotate modular functions by enrichment analysis
enrichmentsForKeyModulesOnly = FALSE # please note that this says if we only want the enrichments for modules that are associated with at least 1 key phenotype.  If positiveCorrelationsOnly is TRUE then this is only for modules positively associated with at least 1  phenotype

performMESHEnrichment = TRUE
MESHCategoryVec = c("A", "B", "C", "D", "E", "F", "G", "H")
meshDatabaseNamesVec = c("gendoo", "gene2pubmed", "RBBH")

performMolSigDBEnrichment = TRUE
molSigDBCategoryVec = c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "H")


keggPathwaysOfInterest <- c("hsa05010", "hsa05171") # what are the IDs of the kegg pathways of interest?



########################################
# Single Nucleotide Polymorphisms (SNPs)
# organizing the Genome-Wide Association Study Data
# Requirements for SNPs and P-value Columns:
# 1.  Please ensure that your dataset has labeled the SNPs column as any of these: SNP, MarkerName, MARKERNAME,
# snp, snpid, SNPID, SNPId, RS, RSID, RS_NUMBER, or RS_NUMBERS
# 2.   Please ensure that your dataset has labeled the P-vlaue column as any of these: P, PVALUE, P_VALUE, PVAL,
# P_VAL, GC_PVALUE, P-VALUE, P-VAL, GC-PVALUE, SNP-P, or SNP-PVALUE
gwasDataFile = "F:\\organizedAlzheimers\\SNPs\\allADGWAS_SNPS.csv"
pValThreshForSNPs = 5e-10 # please note we only consider the SNPs with a p-value below this threshold :)
contextForGWASDataSet = "Jansen et. al 2019 GWAS on Alzheimer's Disease SNPs" # please enter some context to help ye better understand the GWAS dataset used for disease/trait

# Finding the SNPs that disrupt the TF Binding Sites in key regions 
snpsDataBase <- SNPlocs.Hsapiens.dbSNP144.GRCh37 # SNPlocs.Hsapiens.dbSNP151.GRCh38
versionUsed = "SNPlocs.Hsapiens.dbSNP144.GRCh37" # please note that this should be the same as snpsDataBase but should have quotes
# if ye need to use another snpsDataBase (different from this default), please look at:
# BiocManager::install("SNPlocs.Hsapiens.dbSNP144.GRCh37")
# library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
# also Kaviar et. al could help with mapping some SNPs to rsIDs, for example.
searchGenomeName = BSgenome.Hsapiens.UCSC.hg19
genomeUsed = "BSgenome.Hsapiens.UCSC.hg19" # please note that this should be the same as searchGenomeName, but have quotes

# if ye need to use another genomeName (different from this default), please look at:
# BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
# library(BSgenome.Hsapiens.UCSC.hg19)
motifbreakrIntervalLength = 200 # please note this is the interval length for each snp object we consider

motifBreakrThreshold = 1e-3
numWorkersMotifbreakR = 32 # please note that the code can be run in parallel with many workers.  More workers = faster computing, but may be based on your laptop

#threshold = motifBreakrThreshold #, #1e-4,
methodToUseMotifBreakr = "default" # please note that other options include ic for information content, ...
motifBackgroundProbabilitiesVec = c(A=0.25, C=0.25, G=0.25, T=0.25)  # please adjust as needed, but please ensure they add to 1
# please note that if BPPARAM =  BiocParallel::SnowParam(numWorkersMotifbreakR) does not work, please adjust the code as needed
#https://www.bioconductor.org/packages/devel/bioc/vignettes/BiocParallel/inst/doc/Introduction_To_BiocParallel.pdf
# please note that according to page 5 top, SnowParam: Supported on all platforms.
motifBreakRParallelComputing_BPPARAM = BiocParallel::SnowParam(numWorkersMotifbreakR)

chromosomeVec = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6",
                  "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
                  "chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
                  "chr19", "chr20", "chr21", "chr22", "chrX")#, "chrY")

##############################################################################################################################

# please perhaps also have an R data object that stores locations to the files (file paths)

#############################################################################################################################
# Step 5:
# Chromatin Interaction Part:
otherInfo = "" # please add in any additional information that may help for your reference, such as histone data, etc.
haveChromatinInteractionRegulatoryNetwork = FALSE # usually, we assume ye do not have a Chromatin Interaction Regulatory Network from before.  Thus, we will need to construct TF to TG using scgrnom
chromatinInteractionRegulatoryNetworkFilePath = NULL # please fill this in with file path for Chromatin Interaction Regulatory Network, only if you have one
haveEnhancerAndPromoterInteractionsData = FALSE # if ye already have a dataframe (or can make a dataframe in this format) in this format: gene, gene_chr, promoter_start, promoter_end,
# enh_chr,enh_start,enh_end
enhAndPromoterInteractionsFilePath = NULL # please fill this in with file path for interactions dataframe in the form:  gene, gene_chr, promoter_start, promoter_end,
# enh_chr,enh_start,enh_end


haveEnhancersAndPromotersHICData = FALSE# TRUE # do ye have HI-C data on Enhancers and/or Promoters
# this may be optional chromatinInteractionDataFilePath if we do not have hi-c on Chromatin and Enhancer Interactions
chromatinInteractionDataFilePath = "F:\\organizedAlzheimers\\Setup\\InputDataFiles\\chromatinInteractionData\\chromatinInteractionData_Hippocampus.csv"


haveEnhancersInteractionData = TRUE # please note that if haveEnhancersAndPromotersHICData is TRUE, this should also be true. Otherwise, do ye have data on Enhancers: (chrom, start, end) for region of interest?
# please note that we need to at least have this dataframe:

# NECESSARY:
enhancerChromatinFilePath = "C://Users//saniy//Documents//alzheimersResearchChanu//SNPheno//data//chromatinInteractionNetworkData//LateralTemporalLobe//miniDemoEnhancerInteractionsLateralTemporalLobe.csv"

# if we need to get corresponding promoter regions for our enhancer DF, please note that we will run this command
pleaseOrganizeEnhancerRegions = TRUE  # if this is true, we need to organize the enhancer regions
# (ex. from chr1:752537-759037) to separate columns with chrom = 1, enhancer_start = 752537, enhancer_end = 759037
enhancerRegionCol = NULL  # please note that if pleaseOrganizeEnhancerRegions is TRUE, we should specify column for enhancer regions, otherwise, function will try to determine it based on : and - location in row
maxBasesUpstreamForPromoterLength = 5000 #  Promoters given a # of bases of interest for the promoter length.
geneIDColumnNameEnhancerDF = "entrezID"

numOfCoresToUseForScgrnom = 6


#################################################################################################
# Gene Regulatory Networks from Gene Expression Data:
minNumSourcesGeneGRN = 2
hasFinalFullNetworkForChromosome = "FALSE"
filePathForFullNetwork = ""

numOfCoresToUseForGenie3 = 3 # please note that this tells how many cores of computer power we want for running GENIE3

pleasePrintProgressAfterNumIterations_GetTF = 1000 # when using getTF, after how many iterations should we print the progress
includeTheIndividualTFsInGroup = TRUE  # should we also include the individual components of group TFs?  For example: should EWSR1-FLI1 TF be made into 3 elements (TRUE): EWSR1-FLI1, EWSR1, and FLI1?  Or, just left as 1 element (FALSE): EWSR1-FLI1?

numRTNPermutations = 1000 # please input how many permutations ye would like to perform for RTN
weightThresholdGenie3 =  0.0025  #instead of 0.001?

# Trena: please note that there are ways to speed this up:
# default Ensemble solvers: "lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman","xgboost"
numIterationsAfterForSavingTrenaOutput = 1000 # please note that this file will save all input so far from Trena for these number of iterations.  This is to be safe in case R dies or something else happens.
numIterationsToPrintProgressForTrenaOutput = 25 # please note that this is the # of iterations after which ye want a progress update on the # of target genes (in gene expression data) so far that Trena has done.
numOfCoresToUseForSlowestLassoSqrtSolverTrena = 4 # please note this corresponds to TrenaEnsemble parameter of :nCores.sqrt An integer denoting the number of computational cores to devote to the square
#root LASSO solver, which is the slowest of the solvers (default = 4)

speedUpTrenaEnsembleByEliminatingLassoRelatedModels = TRUE #FALSE # please note that if this is true, we do not perform the computationally intensive and slower lasso and lassopv solvers.  Instead, only these solvers are used: "pearson", "randomForest", "ridge", "spearman","xgboost"
useVarianceFilterToReduceCandidateTFsForTrenaToSpeedUpProcess = TRUE #FALSE  # if true, please note that this will apply a variance filter to each target gene to help find best TFs as candidate regulators:
# https://bioconductor.org/packages/release/bioc/vignettes/trena/inst/doc/TReNA_Vignette.html: For instance, we can create a VarianceFilter and use it to find all transcription factors with variance within 50% of the target gene's variance. This will return a named list with both the names of the transcription factors and their variances.
varianceSizeForTFsAndTargetGeneForTrena = 0.25# NULL # please note that if useVarianceFilterToReduceCandidateTFsForTrenaToSpeedUpProcess is TRUE, then this should be a number between 0 and 1



########################################
# Single Nucleotide Polymorphisms (SNPs)
# organizing the Genome-Wide Association Study Data
# Requirements for SNPs and P-value Columns:
# 1.  Please ensure that your dataset has labeled the SNPs column as any of these: SNP, MarkerName, MARKERNAME,
# snp, snpid, SNPID, SNPId, RS, RSID, RS_NUMBER, or RS_NUMBERS
# 2.   Please ensure that your dataset has labeled the P-vlaue column as any of these: P, PVALUE, P_VALUE, PVAL,
# P_VAL, GC_PVALUE, P-VALUE, P-VAL, GC-PVALUE, SNP-P, or SNP-PVALUE
gwasDataFile = "F:\\organizedAlzheimers\\SNPs\\allADGWAS_SNPS.csv"
pValThreshForSNPs = 5e-10# 5e-3 # please note we only consider the SNPs with a p-value below this threshold :)
contextForGWASDataSet = "Jansen et. al 2019 GWAS on Alzheimer's Disease SNPs" # please enter some context to help ye better understand the GWAS dataset used for disease/trait

# Finding the SNPs that disrupt the TF Binding Sites in key regions
snpsDataBase <- SNPlocs.Hsapiens.dbSNP144.GRCh37 # SNPlocs.Hsapiens.dbSNP151.GRCh38
versionUsed = "SNPlocs.Hsapiens.dbSNP144.GRCh37" # please note that this should be the same as snpsDataBase but should have quotes
# if ye need to use another snpsDataBase (different from this default), please look at:
# BiocManager::install("SNPlocs.Hsapiens.dbSNP144.GRCh37")
# library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
# also Kaviar et. al could help with mapping some SNPs to rsIDs, for example.
searchGenomeName = BSgenome.Hsapiens.UCSC.hg19
genomeUsed = "BSgenome.Hsapiens.UCSC.hg19" # please note that this should be the same as searchGenomeName, but have quotes

# if ye need to use another genomeName (different from this default), please look at:
# BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
# library(BSgenome.Hsapiens.UCSC.hg19)
motifbreakrIntervalLength = 200 # please note this is the interval length for each snp object we consider

motifBreakrThreshold = 1e-3
numWorkersMotifbreakR = 32 # please note that the code can be run in parallel with many workers.  More workers = faster computing, but may be based on your laptop

#threshold = motifBreakrThreshold #, #1e-4,
methodToUseMotifBreakr = "default" # please note that other options include ic for information content, ...
motifBackgroundProbabilitiesVec = c(A=0.25, C=0.25, G=0.25, T=0.25)  # please adjust as needed, but please ensure they add to 1
# please note that if BPPARAM =  BiocParallel::SnowParam(numWorkersMotifbreakR) does not work, please adjust the code as needed
#https://www.bioconductor.org/packages/devel/bioc/vignettes/BiocParallel/inst/doc/Introduction_To_BiocParallel.pdf
# please note that according to page 5 top, SnowParam: Supported on all platforms.
motifBreakRParallelComputing_BPPARAM = BiocParallel::SnowParam(numWorkersMotifbreakR)


