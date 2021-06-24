#########################################################################################################################################
### INFORMATION:  Please note that the code below was written by Saniya Khullar and Daifeng Wang, Ph.D.

#########################################################################
## BUILDING FULL GENE REGULATORY NETWORK BASED ON GENE EXPRESSION REGULATORY NETWORK AND CHROMATIN REGULATORY ENTWORK
##################################################################
#########################################################################################################################################
### USER INPUTS:
adsnphenoDirectory = "F://organizedAlzheimers//ADSNPheno//"  # please supply full directory path
codeFolder = paste0(adsnphenoDirectory, "code//") # the main folder where code would be within ADSNPheno
dataFolder = paste0(adsnphenoDirectory, "data//") # the main folder where data would be within ADSNPheno

outputPathNameADSNPhenoOutputs = paste0("F://organizedAlzheimers//ADSNPhenoOutputs_NEW//")

inputGeneExpressionDataFilePath = paste0(dataFolder, "originalGeneExpressionDataSets//LateralTemporalLobeRegionGeneExpressionData_miniDemo_200genes.csv") #"F://organizedAlzheimers//Setup//InputDataFiles//originalGeneExpressionDataSets//LateralTemporalLobeRegionGeneExpressionData_miniDemo.csv" #"F://organizedAlzheimers//Setup//InputDataFiles//originalGeneExpressionDataSets//originalLateralTemporalLobeRegionGeneExpressionData.csv"

disease = "Alzheimers"
tissueName = "Brain"
bodyRegion = "MiniFinalDemo_LTL1"

tfsUsed = "LambertAndJaspar"
numberOfRoundingDigits = 3
# K-Means step performed after WGCNA
performAdditionalKMeansStep = TRUE # or can be FALSE if you want original WGCNA results

phenotypesFilePath = paste0(dataFolder, "phenotypes//AlzheimersLateralTemporalLobePhenotypesUpdated.csv")

log2transformInputData = "TRUE" # (should we apply a log2(x+1) transform on data
scaleInputData = "FALSE"  # should we apply a scale() on data

includeTimeStamp = "TRUE"
netType = "signed"  # # network type:but please check WGCNA for additional parameters.  We retained the defaults

####### GENE REGULATORY NETWORKS CODE:

#################################################################################################
# Gene Regulatory Networks from Gene Expression Data:
minNumSourcesGeneGRN = 2
hasFinalFullNetworkForChromosome = "FALSE"
filePathForFullNetwork = ""


includeTheIndividualTFsInGroup = TRUE  # should we also include the individual components of group TFs?  For example: should EWSR1-FLI1 TF be made into 3 elements (TRUE): EWSR1-FLI1, EWSR1, and FLI1?  Or, just left as 1 element (FALSE): EWSR1-FLI1?


# RTN Gene Regulatory Network Inputs:
enrichmentsForKeyModulesOnly = FALSE # please note that this says if we only want the enrichments for modules that are associated with at least 1 key phenotype.  If positiveCorrelationsOnly is TRUE then this is only for modules positively associated with at least 1  phenotype

# GENIE3 Gene Regulatory Network Inputs:
weightThresholdGenie3 =  0.0025

progressBarPythonIterations = 5000



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
includeTheIndividualTFsInGroup_Python = "Yes"  # should we also include the individual components of group TFs?  For example: should EWSR1-FLI1 TF be made into 3 elements (TRUE): EWSR1-FLI1, EWSR1, and FLI1?  Or, just left as 1 element (FALSE): EWSR1-FLI1?

pleasePrintProgressAfterNumIterations_GetTF = 1000 # when using getTF, after how many iterations should we print the progress

numOfCoresToUseForScgrnom = 6

enhancer_buffer_kbp = 10 # kilobase pairs to add to the enhancer on each side
promoter_buffer_kbp = 2 # kilobase pairs to add to the promoter


########################################################################################################################################
### PLEASE RUN:
########################################################### running other things:
generalDataForPipeline = paste0(dataFolder, "generalDataForPipeline.RData")
load(generalDataForPipeline)
pathForPythonCode = paste0(codeFolder, "setupScripts//BuildingFullAndSubNetworksPythonCode.py")
#filePathDerivationsSourceCode = paste0(codeFolder, "setupScripts//filePathDerivations.R")
packagesNeededSourceCode = paste0(codeFolder, "setupScripts//packagesNeeded.R")
functionsNeededSourceCode = paste0(codeFolder, "setupScripts//functionsNeeded.R")
userInputsSourceCode = paste0(codeFolder, "pythonCode//BuildingFullAndSubNetworksPythonCode.py")
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


initialRDataFilePath_FilePathDerivationsAndInfoObjectsForWGCNA = pleaseGetInitialFilePathDerivationsAndInfoObjectsForWGCNA(outputPathNameADSNPhenoOutputs, inputGeneExpressionDataFilePath,
                                                                                                                           disease, tissueName, bodyRegion, tfsUsed, performAdditionalKMeansStep, phenotypesFilePath, log2transformInputData, scaleInputData, includeTimeStamp, numberOfRoundingDigits, netType)

load(initialRDataFilePath_FilePathDerivationsAndInfoObjectsForWGCNA)
load(wgcnaWithKmeansAllObjectsExceptTOM_SimpleFileName) # please note that after we load initialRDataFilePath_FilePathDerivationsAndInfoObjectsForWGCNA, we have names for many file paths and can load this to get WGCNA results :) YAY!


geneRegulatoryPathwayFileNamesRData = pleaseCreateGeneExpressionRegulatoryNetworkFileNames(disease, tissueName, bodyRegion, outputPathNameADSNPhenoOutputs, dataScalingOutputMini, outputAddOn, tfsUsed, weightThresholdGenie3)
load(geneRegulatoryPathwayFileNamesRData)

###########################
# ## Combining Gene Regulatory Networks:
combinedGeneExpressGRNDF = pleaseCreateCombined_GeneRegNetwork_From_RTN_Genie3_TrenaEnsembleSolver_and_Trrust2(trrust2DF,
                                                                                                               RTN_rDataPart,
                                                                                                               genie3RDataPartFinal,
                                                                                                               trenaOutpathALLRData_Final,
                                                                                                               grnResultsCombinedFileNameCSV,
                                                                                                               grnResultsCombinedFileNameRData,
                                                                                                               pathForPythonCode,
                                                                                                               filePathOfAll4GRN, disease,
                                                                                                               bodyRegion,
                                                                                                               dataScalingOutputMini,
                                                                                                               geneToEntrezIDMappingPath,
                                                                                                               minNumSourcesGeneGRN,
                                                                                                               progressBarPythonIterations)
parentName = paste0(transcriptionalRegulatoryNetworkOutputPath, "//")# #'F://organizedAlzheimers//FinalGeneRegulatoryNetworks//Alzheimers_Brain_Hippocampus//GeneRegulatNetwork//'

chromatinInteractionFilesRData = pleaseCreateChromatinRegulatoryNetworkFileNames(disease, tissueName, bodyRegion, outputPathNameADSNPhenoOutputs, haveEnhancersAndPromotersHICData)
load(chromatinInteractionFilesRData)



if (haveChromatinInteractionRegulatoryNetwork){
  hasFinalFullNetworkForChromosome = "TRUE"
} else {
  hasFinalFullNetworkForChromosome = "FALSE"
}

scgrnGetTFsForEachChromosomesChromatinNetworkVec = paste0(chromatinInteraction_Step2_FilePath_CleanedCSVObjects,
                                                          list.files(chromatinInteraction_Step2_FilePath_CleanedCSVObjects, pattern = "scgrnomGetTFs_"))
scgrnGetTFsForEachChromosomesChromatinNetworkVec
chromatinRegNetFileNameList = scgrnGetTFsForEachChromosomesChromatinNetworkVec
load(powerRelatedFilePathsRData)
source_python(pathForPythonCode)
finalGRN = pleaseGetOrganizedFinalGeneRegulatoryNetwork(geneAssignmentsModAssignmentsPhenosWithKmeansOutputNameCSV,
                                                        finalFullNetworkChromatinAndGeneExpressCSV,
                                                        parentName, disease, bodyRegion,
                                                        dataScalingOutputMini, minNumSourcesGeneGRN,
                                                        chromatinRegNetFileNameList, grnResultsCombinedFileNameCSV,
                                                        geneToEntrezIDMappingPath,
                                                        enhancer_buffer_kbp, promoter_buffer_kbp,
                                                        powerEstimate,
                                                        progressBarPythonIterations,
                                                        hasFinalFullNetworkForChromosome,
                                                        filePathForFullNetwork)

finalGRN




# fullNetworkDFForChromFunction(scgrnomOutputFileName, grnResultsCombinedFileNameCSV, parentName, disease, bodyRegion, dataScalingOutputMini, geneToEntrezIDMappingPath, minNumSourcesGeneGRN, enhancer_buffer_kbp, promoter_buffer_kbp, progressBarPythonIterations = 1000)
#
# pleaseGetOrganizedFinalGeneRegulatoryNetwork(geneAssignmentsModAssignmentsPhenosWithKmeansOutputNameCSV, finalFullNetworkChromatinAndGeneExpressCSV, parentName, disease, bodyRegion, dataScalingOutputMini, minNumSourcesGeneGRN, chromatinRegNetFileNameList, grnResultsCombinedFileNameCSV, geneToEntrezIDMappingPath,
#                                              enhancer_buffer_kbp, promoter_buffer_kbp, powerEstimate, progressBarPythonIterations, hasFinalFullNetworkForChromosome = "FALSE", filePathForFullNetwork = "")




