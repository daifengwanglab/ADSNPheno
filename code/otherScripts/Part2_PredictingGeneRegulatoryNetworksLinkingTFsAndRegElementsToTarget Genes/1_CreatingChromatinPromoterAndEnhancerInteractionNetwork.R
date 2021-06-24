#########################################################################################################################################
### INFORMATION:  Please note that the code below was written by Saniya Khullar and Daifeng Wang, Ph.D.
#########################################################################################################################################
### USER INPUTS:

adsnphenoDirectory = "F://organizedAlzheimers//ADSNPheno//"  # please supply full directory path

disease = "Alzheimers"
tissueName = "Brain"
bodyRegion = "MiniFinalDemo_LTL1"
codeFolder = paste0(adsnphenoDirectory, "code//") # the main folder where code would be within ADSNPheno
dataFolder = paste0(adsnphenoDirectory, "data//") # the main folder where data would be within ADSNPheno

outputPathNameADSNPhenoOutputs = paste0("F://organizedAlzheimers//ADSNPhenoOutputs_NEW//")

chromosomeVec = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6",
                  "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
                  "chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
                  "chr19", "chr20", "chr21", "chr22", "chrX")#, "chrY")
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
chromatinInteractionDataFilePath = NULL #

haveEnhancersInteractionData = TRUE # please note that if haveEnhancersAndPromotersHICData is TRUE, this should also be true. Otherwise, do ye have data on Enhancers: (chrom, start, end) for region of interest?
# please note that we need to at least have this dataframe:



# NECESSARY:
enhancerChromatinFilePath = paste0(dataFolder, "chromatinInteractionNetworkData//miniDemoEnhancerInteractionsLateralTemporalLobe.csv")

# if we need to get corresponding promoter regions for our enhancer DF, please note that we will run this command
pleaseOrganizeEnhancerRegions = TRUE  # if this is true, we need to organize the enhancer regions
# (ex. from chr1:752537-759037) to separate columns with chrom = 1, enhancer_start = 752537, enhancer_end = 759037
enhancerRegionCol = NULL  # please note that if pleaseOrganizeEnhancerRegions is TRUE, we should specify column for enhancer regions, otherwise, function will try to determine it based on : and - location in row
maxBasesUpstreamForPromoterLength = 5000 #  Promoters given a # of bases of interest for the promoter length.
geneIDColumnNameEnhancerDF = "entrezID"

includeTheIndividualTFsInGroup = TRUE  # should we also include the individual components of group TFs?  For example: should EWSR1-FLI1 TF be made into 3 elements (TRUE): EWSR1-FLI1, EWSR1, and FLI1?  Or, just left as 1 element (FALSE): EWSR1-FLI1?
pleasePrintProgressAfterNumIterations_GetTF = 1000 # when using getTF, after how many iterations should we print the progress

numOfCoresToUseForScgrnom = 6


#####################################################################################
# 5.  Please predict gene regulatory networks for genes and modules
# Step 5:
################################################################################
# Chromatin Interaction Part:
# Please note that since the chromosomes are usually smaller in size the higher the number of the chromosome,
# we start with the later chromosomes to ensure those work :)
# # Step2: please get TFs for each promoter and enhancer
# # https://github.com/daifengwanglab/scGRN

##########################################################################################
# PLEASE NOTE THIS:
# user can have the option to create an interactionsDF or to run scgrnom to get it
# caseInfo
########################################################################################################################################
### PLEASE RUN:
########################################################### running other things:
generalDataForPipeline = paste0(dataFolder, "generalDataForPipeline.RData")
load(generalDataForPipeline)
pathForPythonCode = paste0(codeFolder, "setupScripts//BuildingFullAndSubNetworksPythonCode.py")
#filePathDerivationsSourceCode = paste0(codeFolder, "setupScripts//filePathDerivations.R")
packagesNeededSourceCode = paste0(codeFolder, "setupScripts//packagesNeeded.R")
functionsNeededSourceCode = paste0(codeFolder, "setupScripts//functionsNeeded.R")
#userInputsSourceCode = paste0(codeFolder, "pythonCode//BuildingFullAndSubNetworksPythonCode.py")
generalDataForPipeline = paste0(dataFolder, "generalDataForPipeline.RData")

source(packagesNeededSourceCode)
source(functionsNeededSourceCode)
load(generalDataForPipeline)
#source(filePathDerivationsSourceCode)
source_python(pathForPythonCode)

print(paste0(":) Please note: packagesNeededSourceCode: ", packagesNeededSourceCode))
print(paste0(":) Please note: functionsNeededSourceCode: ", functionsNeededSourceCode))
print(paste0(":) Please note: generalDataForPipeline: ", generalDataForPipeline))
#print(paste0(":) Please note: filePathDerivationsSourceCode: ", filePathDerivationsSourceCode))
print(paste0(":) Please note: pathForPythonCode: ", pathForPythonCode))



chromatinInteractionFilesRData = pleaseCreateChromatinRegulatoryNetworkFileNames(disease, tissueName, bodyRegion, outputPathNameADSNPhenoOutputs, haveEnhancersAndPromotersHICData)

load(chromatinInteractionFilesRData)

scgrnGetTFsForEachChromosomesChromatinNetworkVec = pleaseBuildChromatinGeneRegulatoryNetwork(numOfCoresToUseForScgrnom,
                                                                                             haveChromatinInteractionRegulatoryNetwork,
                                                                                             haveEnhancerAndPromoterInteractionsData,
                                                                                             haveEnhancersAndPromotersHICData,
                                                                                             bodyRegion, otherInfo, folderName,
                                                                                             chromatinInteraction_Step2_FilePath,
                                                                                             chromatinInteraction_Step2_FilePath_InitialRDataObjects,
                                                                                             chromatinInteraction_Step2_FilePath_CleanedRDataObjects,
                                                                                             chromatinInteraction_Step2_FilePath_CleanedCSVObjects,
                                                                                             includeTheIndividualTFsInGroup,
                                                                                             pleasePrintProgressAfterNumIterations_GetTF,
                                                                                             chromatinInteraction_Step1ADSNPheno_FilePath,
                                                                                             maxBasesUpstreamForPromoterLength, enhancerChromatinFilePath,
                                                                                             geneIDColumnNameEnhancerDF, enhancerRegionCol,
                                                                                             pleaseOrganizeEnhancerRegions,
                                                                                             chromatinInteractionRegulatoryNetworkFilePath,
                                                                                             enhAndPromoterInteractionsFilePath,
                                                                                             chromatinInteractionDataFilePath,
                                                                                             enhancerAndPromoterInteraction)



numOfCoresToUseForScgrnom
haveChromatinInteractionRegulatoryNetwork
haveEnhancerAndPromoterInteractionsData
haveEnhancersAndPromotersHICData
bodyRegion
otherInfo
folderName
chromatinInteraction_Step2_FilePath
chromatinInteraction_Step2_FilePath_InitialRDataObjects
chromatinInteraction_Step2_FilePath_CleanedRDataObjects
chromatinInteraction_Step2_FilePath_CleanedCSVObjects
includeTheIndividualTFsInGroup
pleasePrintProgressAfterNumIterations_GetTF
chromatinInteraction_Step1ADSNPheno_FilePath
maxBasesUpstreamForPromoterLength
enhancerChromatinFilePath
geneIDColumnNameEnhancerDF
enhancerRegionCol
pleaseOrganizeEnhancerRegions
chromatinInteractionRegulatoryNetworkFilePath
enhAndPromoterInteractionsFilePath
chromatinInteractionDataFilePath
enhancerAndPromoterInteraction
