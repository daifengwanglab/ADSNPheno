#########################################################################################################################################
### INFORMATION:  Please note that the code below was written by Saniya Khullar and Daifeng Wang, Ph.D.

#########################################################################
## GENE EXPRESSION GENE REGULATORY NETWORKS:
##################################################################
#########################################################################################################################################
### USER INPUTS:
adsnphenoDirectory = "F://organizedAlzheimers//ADSNPheno//"  # please supply full directory path
codeFolder = paste0(adsnphenoDirectory, "code//") # the main folder where code would be within ADSNPheno
dataFolder = paste0(adsnphenoDirectory, "data//") # the main folder where data would be within ADSNPheno

outputPathNameADSNPhenoOutputs = paste0("F://organizedAlzheimers//ADSNPhenoOutputs_NEW//")

inputGeneExpressionDataFilePath = paste0(dataFolder, "//originalGeneExpressionDataSets//LateralTemporalLobeRegionGeneExpressionData_miniDemo_200genes.csv") #"F://organizedAlzheimers//Setup//InputDataFiles//originalGeneExpressionDataSets//LateralTemporalLobeRegionGeneExpressionData_miniDemo.csv" #"F://organizedAlzheimers//Setup//InputDataFiles//originalGeneExpressionDataSets//originalLateralTemporalLobeRegionGeneExpressionData.csv"

disease = "Alzheimers"
tissueName = "Brain"
bodyRegion = "MiniFinalDemo_LTL1"

tfsUsed = "LambertAndJaspar"
numberOfRoundingDigits = 3
# K-Means step performed after WGCNA
performAdditionalKMeansStep = TRUE # or can be FALSE if you want original WGCNA results

phenotypesFilePath = paste0(dataFolder, "//phenotypes//AlzheimersLateralTemporalLobePhenotypesUpdated.csv")

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

numOfCoresToUseForGenie3 = 3 # please note that this tells how many cores of computer power we want for running GENIE3

pleasePrintProgressAfterNumIterations_GetTF = 1000 # when using getTF, after how many iterations should we print the progress
includeTheIndividualTFsInGroup = TRUE  # should we also include the individual components of group TFs?  For example: should EWSR1-FLI1 TF be made into 3 elements (TRUE): EWSR1-FLI1, EWSR1, and FLI1?  Or, just left as 1 element (FALSE): EWSR1-FLI1?


# RTN Gene Regulatory Network Inputs:
numRTNPermutations = 1000 # please input how many permutations ye would like to perform for RTN

# GENIE3 Gene Regulatory Network Inputs:
weightThresholdGenie3 =  0.0025


# TReNA Ensemble Solver Gene Regulatory Network Inputs:
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
print(paste0(":) Please note: filePathDerivationsSourceCode: ", filePathDerivationsSourceCode))
print(paste0(":) Please note: pathForPythonCode: ", pathForPythonCode))


initialRDataFilePath_FilePathDerivationsAndInfoObjectsForWGCNA = pleaseGetInitialFilePathDerivationsAndInfoObjectsForWGCNA(outputPathNameADSNPhenoOutputs, inputGeneExpressionDataFilePath,
                                                                                                                           disease, tissueName, bodyRegion, tfsUsed, performAdditionalKMeansStep, phenotypesFilePath, log2transformInputData, scaleInputData, includeTimeStamp, numberOfRoundingDigits, netType)

load(initialRDataFilePath_FilePathDerivationsAndInfoObjectsForWGCNA)
load(wgcnaWithKmeansAllObjectsExceptTOM_SimpleFileName) # please note that after we load initialRDataFilePath_FilePathDerivationsAndInfoObjectsForWGCNA, we have names for many file paths and can load this to get WGCNA results :) YAY!



print(paste0(":) please note body region = ", bodyRegion))

#########################################################################
## GENE EXPRESSION GENE REGULATORY NETWORKS:
##################################################################
load(powerRelatedFilePathsRData)
#load(wgcnaWithKmeansAllObjectsExceptTOM)


grnInfoNeeded = gettingInformationForGeneExpressionGeneRegulatoryNetworks(tfsInGeneExpressionDataRData,
                                                                          entrezMappingFilePath,
                                                                          tfsInGeneExpressionDataCSV,
                                                                          filePathOfGeneExpressionDataUsedByWGCNA,
                                                                          inputGeneExpressionDataFilePath,
                                                                          geneToEntrezIDPathForGeneExpressionDataCSV,
                                                                          geneToEntrezIDPathForGeneExpressionDataRData, tfsDF)


tfsToUse = grnInfoNeeded$tfsToUse
allsamples = grnInfoNeeded$allsamples
geneExpressionSamplesVec = grnInfoNeeded$geneExpressionSamplesVec
tfsInGeneExpressionDataRData = grnInfoNeeded$tfsInGeneExpressionDataRData

colnamesVec = c("TF", "RegulatedGene", "Source", "Info", "CombinedName")
tissueVec = rep(tissueName, ncol(allsamples))


#####################################################################
############## RTN Gene Transcriptional Regulatory Network:

rtnRegulonsDF = pleaseGet_RTN_GeneRegulatoryNetwork(tfsInGeneExpressionDataRData,
                                                    numRTNPermutations,
                                                    oldRDataPart, newRDataPart,
                                                    numberOfRoundingDigits,
                                                    RTN_rDataPart, csvPartRTNFinal,
                                                    csvPartRTNFinal_NoDate)


print(paste0(":) please note that next we can find the TFs that significantly regulate the gene modules:"))

load(RTN_rDataPart)
if (enrichmentsForKeyModulesOnly){
  # we only load data for the key modules then
  load(keyGeneAssignmentsWithKmeansOutputNameRData)
} else {
  # we load info on genes in all modules
  load(geneAssignmentsWithKmeansOutputNameRData)
}

moduleNamesVec = unique(clusterInfo$module)

tni.regulon.summary(rtni)

##############################################################
# Finding Transcriptional Regulatory Network for TFs Regulating Gene Modules
mraAllModulesDF = findingTFsThatRegulateGeneModules(clusterInfo, folderName,
                                                    geneModMembershipOutputNameCSV,
                                                    filePathOfTFtoModule, RTN_rDataPart,
                                                    geneAndEntrezIDMappingDF, tfsDF,
                                                    wgcnaWithKmeansAllObjectsExceptTOM,
                                                    enrichmentsForKeyModulesOnly)


###############################################################
# Building Transcriptional Regulatory Network Using GENIE3

# GENIE3:
###

weightThreshold = weightThresholdGenie3

load(tfsInGeneExpressionDataRData) # loading info on TFsToUse

genie3df = buildGRNUsingGenie3(allsamples, tfsToUse, weightThresholdGenie3,
                               genie3_rDataPart, genie3RDataPartFinal,
                               numOfCoresToUseForGenie3)



###############################################################
# Building Transcriptional Regulatory Network Using TRENA Ensemble Solver
## TRENA:

trenaEnsembleDF = buildGRNUsingTrenaEnsembleSolver(allsamples, tfsToUse, powerEstimate, tfsInGeneExpressionDataRData,
                                                   speedUpTrenaEnsembleByEliminatingLassoRelatedModels,
                                                   useVarianceFilterToReduceCandidateTFsForTrenaToSpeedUpProcess,
                                                   varianceSizeForTFsAndTargetGeneForTrena,
                                                   filePathOfTrena,
                                                   trenaOutpathALLCSV,
                                                   trenaOutpathALLRData,
                                                   trenaOutpathALLCSV_Final,
                                                   trenaOutpathALLRData_Final,
                                                   filePathOfGeneExpressionDataUsedByWGCNA,
                                                   numIterationsAfterForSavingTrenaOutput,
                                                   numIterationsToPrintProgressForTrenaOutput)
