
#########################################################################################################################################
### INFORMATION:  Please note that the code below was written by Saniya Khullar and Daifeng Wang, Ph.D.
#########################################################################################################################################
### USER INPUTS:
###############################################################################################################################################################################
######## SNPs breaking TF Binding Sites:

########################################
# Single Nucleotide Polymorphisms (SNPs)
outputPathNameADSNPhenoOutputs = paste0("F://organizedAlzheimers//ADSNPhenoOutputs_NEW//")
disease = "Alzheimers"
tissueName = "Brain"
bodyRegion = "MiniFinalDemo_LTL1"



adsnphenoDirectory = "F://organizedAlzheimers//ADSNPheno//"  # please supply full directory path
codeFolder = paste0(adsnphenoDirectory, "code//") # the main folder where code would be within ADSNPheno
dataFolder = paste0(adsnphenoDirectory, "data//") # the main folder where data would be within ADSNPheno


packagesNeededSourceCode = paste0(codeFolder, "setupScripts//packagesNeeded.R")
source(packagesNeededSourceCode)
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



########################################################################################################################################
### PLEASE RUN:
########################################################### running other things
generalDataForPipeline = paste0(dataFolder, "generalDataForPipeline.RData")
load(generalDataForPipeline)
pathForPythonCode = paste0(codeFolder, "setupScripts//BuildingFullAndSubNetworksPythonCode.py")
#filePathDerivationsSourceCode = paste0(codeFolder, "setupScripts//filePathDerivations.R")
functionsNeededSourceCode = paste0(codeFolder, "setupScripts//functionsNeeded.R")
userInputsSourceCode = paste0(codeFolder, "pythonCode//BuildingFullAndSubNetworksPythonCode.py")
generalDataForPipeline = paste0(dataFolder, "generalDataForPipeline.RData")

source(functionsNeededSourceCode)
load(generalDataForPipeline)
#source(filePathDerivationsSourceCode)
source_python(pathForPythonCode)

print(paste0(":) Please note: packagesNeededSourceCode: ", packagesNeededSourceCode))
print(paste0(":) Please note: functionsNeededSourceCode: ", functionsNeededSourceCode))
print(paste0(":) Please note: generalDataForPipeline: ", generalDataForPipeline))
#print(paste0(":) Please note: filePathDerivationsSourceCode: ", filePathDerivationsSourceCode))
print(paste0(":) Please note: pathForPythonCode: ", pathForPythonCode))


snpsOutputPathRDataObjects = pleaseCreateSNPsMotifBreakrFileNames(outputPathNameADSNPhenoOutputs, pValThreshForSNPs, disease, tissueName, bodyRegion)

load(snpsOutputPathRDataObjects)
motifBreakRSNPsDF = findingSNPsBreakingTFBindingSitesFromGWASDataFunction(gwasDataFile,
                                                                          pValThreshForSNPs,
                                                                          contextForGWASDataSet,
                                                                          updatedGWASDataSetFileNameALL_CSV,
                                                                          updatedGWASDataSetFileNameALL_RData,
                                                                          updatedGWASDataSetFileNameSNPandP_CSV,
                                                                          updatedGWASDataSetFileNameSNPandP_RData,
                                                                          motifbreakrIntervalLength,
                                                                          filePathOfMotifBreakR_Step1,
                                                                          filePathOfMotifBreakR_Step2,
                                                                          versionUsed, genomeUsed,fullDiseaseName,
                                                                          motifBreakrThreshold,
                                                                          methodToUseMotifBreakr,
                                                                          motifBackgroundProbabilitiesVec,
                                                                          motifBreakRParallelComputing_BPPARAM,
                                                                          snpDFPath,
                                                                          snpDFPathRData)
