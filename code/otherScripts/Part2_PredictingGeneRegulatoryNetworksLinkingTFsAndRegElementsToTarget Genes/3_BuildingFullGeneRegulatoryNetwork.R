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


