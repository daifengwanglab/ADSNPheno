# File Path Derivations




pleaseGetInitialFilePathDerivationsAndInfoObjectsForWGCNA <- function(outputPathNameADSNPhenoOutputs,
disease, tissueName,
bodyRegion, tfsUsed,
performAdditionalKMeansStep,
phenotypesFilePath, log2transformInputData,
scaleInputData, includeTimeStamp){
print(":) Please note that this function takes in initial user input and helps us store the file paths and other key variables that we may need for our analysis.") 

pipeline = "ADSNPheno"
fullDiseaseName = paste0(disease, "Disease")
print(paste0(":) please note body region = ", bodyRegion))
library(stringr)

print(paste0(":) please note the key output directory (folder) is: outputPathNameADSNPhenoOutputs = ", outputPathNameADSNPhenoOutputs))
# please create a folder within ADSNPhenoOutputs, which will have the results of pipeline for Hippocampus Ca1 here
folderName = str_replace_all(paste0(disease, "_", tissueName, "_", bodyRegion), " ", "")
dir.create(folderName)
print(paste0("outputPathNameADSNPhenoOutputs: ", outputPathNameADSNPhenoOutputs))

mainOutputFolder = paste0(
  outputPathNameADSNPhenoOutputs, folderName)
dir.create(mainOutputFolder, showWarnings = FALSE)

print(paste0(":) Please note that a directory (folder) has been created here: ",
             mainOutputFolder, " for all of the output of ADSNPheno"))



filePathInfoDF_FileNameRData = paste0(mainOutputFolder, "//fileInfoDF_", pipeline, "_",folderName, ".RData")
filePathInfoDF_FileNameCSV = paste0(mainOutputFolder, "//fileInfoDF_", pipeline, "_",folderName, ".csv")
tfsInGeneExpressionDataCSV = paste0(mainOutputFolder, "//tfsInGeneExpressionDataset_",tfsUsed, "_", folderName, ".csv")
tfsInGeneExpressionDataRData = paste0(mainOutputFolder, "//tfsInGeneExpressionDataset_",tfsUsed, "_", folderName, ".RData")

if(performAdditionalKMeansStep){
  wgcnaAndKMeansOutputPath = paste0(mainOutputFolder, "//", "WGCNAandKMeansAfter")
  print(paste0("Please note that since performAdditionalKMeansStep = ", performAdditionalKMeansStep,
               " we will store the results of the Co-Expression Network by WGCNA with an additional K-Means step after",
               " in this directory (folder): ", wgcnaAndKMeansOutputPath))

  wgcnaToAdd = "//WGCNAWithKMeans"

} else {
  wgcnaAndKMeansOutputPath = paste0(mainOutputFolder, "//", "WGCNA_Only_Outputs")
  print(paste0("Please note that since performAdditionalKMeansStep = ", performAdditionalKMeansStep,
               " we will store the results of the Co-Expression Network by WGCNA only ",
               " in this directory (folder): ", wgcnaAndKMeansOutputPath))

  wgcnaToAdd = "//WGCNA_Only_"

}
dir.create(wgcnaAndKMeansOutputPath, showWarnings = FALSE)
datTraitsRDataFilePath = paste0(mainOutputFolder, "//datTraitsPhenotypeInfo_", folderName, ".RData")


# file path for storing filePaths related to the power used for WGCNA
powerRelatedFilePathsRData = paste0(wgcnaAndKMeansOutputPath, "//rDataOfTheFilePathInfoBasedOnWGCNAPower_", folderName, ".RData")
print(paste0(":) Please note that we will store the information on the file paths for the RData related to the WGCNA power selected for ",
             folderName, " here: ", powerRelatedFilePathsRData))
			 
# file path for the gene expression data

geneExpressionPreparationFilePath = paste0(wgcnaAndKMeansOutputPath, "//GeneExpressionPreparation")
print(paste0("geneExpressionPreparationFilePath: ", geneExpressionPreparationFilePath))
dir.create(geneExpressionPreparationFilePath, showWarnings = FALSE)

#final gene expression:
finalGeneExpressionPreparationFilePath = paste0(geneExpressionPreparationFilePath, "//finalGeneExpressionData")
print(paste0("finalGeneExpressionPreparationFilePath: ", finalGeneExpressionPreparationFilePath))
dir.create(finalGeneExpressionPreparationFilePath, showWarnings = FALSE)


# file path for the entrez ID mapping info
entrezMappingFilePath = paste0(wgcnaAndKMeansOutputPath, "//EntrezMappingInfo")
print(paste0("entrezMappingFilePath: ", entrezMappingFilePath))
dir.create(entrezMappingFilePath, showWarnings = FALSE)

# file path for the soft thresholding power info
softThresholdPowerFilePath = paste0(wgcnaAndKMeansOutputPath, "//SoftThresholdPowerAnalysis")
print(paste0("softThresholdPowerFilePath: ", softThresholdPowerFilePath))
dir.create(softThresholdPowerFilePath, showWarnings = FALSE)


# file path for the Topological Overlap Matrix Info
tomFilePath = paste0(wgcnaAndKMeansOutputPath, "//TopologicalOverlapMatrix")
print(paste0("tomFilePath: ", tomFilePath))
dir.create(tomFilePath, showWarnings = FALSE)


# file path for WGCNA only
wgcnaFilePath = paste0(wgcnaAndKMeansOutputPath, "//WGCNAOnly")
print(paste0("wgcnaFilePath: ", wgcnaFilePath))
dir.create(wgcnaFilePath, showWarnings = FALSE)	


# file path for WGCNA and Kmeans after
if (performAdditionalKMeansStep){
  wgcnaWithKMeansAfterFilePath = paste0(wgcnaAndKMeansOutputPath, "//WGCNAThenKMeans")
  pathForWGCNA = wgcnaWithKMeansAfterFilePath
  print(paste0("wgcnaWithKMeansAfterFilePath: ", wgcnaWithKMeansAfterFilePath))
  dir.create(wgcnaWithKMeansAfterFilePath, showWarnings = FALSE)
  print(paste0(":) please note that the outputs for ", bodyRegion,
               " WGCNA Gene Co-Expression Analysis with K-Means after will be written here: ", wgcnaWithKMeansAfterFilePath))

} else {
  print(paste0(":) since performAdditionalKMeansStep is False, no K-means step will be performed afterwards.  Instead, please note that the outputs for ", bodyRegion,
               " WGCNA Gene Co-Expression Analysis Only after will be written here: ", wgcnaFilePath))
  pathForWGCNA = wgcnaFilePath

}

MEsFilePath = paste0(pathForWGCNA, "//ModuleEigengenes")
moduleAssignmentsFilePath = paste0(pathForWGCNA, "//ModuleAssignments")
geneModuleMembershipFilePath = paste0(pathForWGCNA, "//GeneModuleMembership")

allObjectsFilePath = paste0(pathForWGCNA, "//AllObjects")
dir.create(allObjectsFilePath, showWarnings = FALSE)


phenotypeAssociationFilePath = paste0(pathForWGCNA, "//PhenoCorrs")
dir.create(phenotypeAssociationFilePath, showWarnings = FALSE)

modulePhenoAssociationFilePath = paste0(phenotypeAssociationFilePath, "//ModulePhenoCorrs")
genePhenoAssociationFilePath = paste0(phenotypeAssociationFilePath, "//GenePhenoCorrs")
dir.create(modulePhenoAssociationFilePath, showWarnings = FALSE)
dir.create(genePhenoAssociationFilePath, showWarnings = FALSE)



print(paste0(":) MEsFilePath: ", MEsFilePath))
print(paste0(":) moduleAssignmentsFilePath: ", moduleAssignmentsFilePath))
print(paste0(":) geneModuleMembershipFilePath: ", geneModuleMembershipFilePath))

print(paste0(":) phenotypeAssociationFilePath: ", phenotypeAssociationFilePath))
print(paste0(":) modulePhenoAssociationFilePath: ", modulePhenoAssociationFilePath))
print(paste0(":) genePhenoAssociationFilePath: ", genePhenoAssociationFilePath))



dir.create(MEsFilePath, showWarnings = FALSE)
dir.create(moduleAssignmentsFilePath, showWarnings = FALSE)
dir.create(geneModuleMembershipFilePath, showWarnings = FALSE)

# Step 2:  How to construct co-expression network and modules (WGCNA and K-means too).  Define parameters
#wgcnaAndKMeansOutputPath = "F://organizedAlzheimers//GeneExpressionPreparation//outputs"

# :) Please note that this R file by Saniya has all the parameters we need

diseaseTraits = read.csv(phenotypesFilePath,
                         header = TRUE)

newDate = str_replace_all(Sys.Date(), "-", "_")



######### Please note the function
if (log2transformInputData == "TRUE"){
  if (scaleInputData == "TRUE"){
    outputAddOn = "_log2AndScaledInput"
    print(":) Please note we will apply a Log2(x+1) transformation on the input gene expression data and then will scale that.")
    print(":) This final data set will be our input for WGCNA")
    dataScaling = "log2TransformedAndScaledGeneExpression" #  <-- will come from WGCNA
    dataScalingOutputMini = "Log2AndScaledGeneExpress"

  } else {
    outputAddOn = "_log2InputOnly"
    print(":) Please note we will apply a Log2(x+1) transformation on the input gene expression data.")
    print(":) This final data set will be our input for WGCNA")
    dataScaling = "log2TransGeneExpOnly" #  <-- will come from WGCNA
    dataScalingOutputMini = "Log2TransformedGeneExpress"
  }
} else if (scaleInputData == "TRUE"){
  outputAddOn = "_ScaledInputOnly"
  print(":) Please note we will apply a Scale transformation on the input gene expression data.")
  print(":) This final data set will be our input for WGCNA")
  dataScaling = "ScaledGeneExpOnly" #  <-- will come from WGCNA
  dataScalingOutputMini = "ScaledGeneExpressInput"
} else {
  outputAddOn = "_OriginalInputData"
  print(":) Please note we will use our original input gene expression dataset for WGCNA.")
  dataScaling = "originalGeneExpression" #  <-- will come from WGCNA
  dataScalingOutputMini = "OriginalGeneExpress"


}


geneToEntrezIDPathForGeneExpressionDataRData = paste0(entrezMappingFilePath, "//geneToEntrezIDMapping_",
                                                      disease, "_", bodyRegion, ".RData")


geneToEntrezIDPathForGeneExpressionDataCSV = paste0(entrezMappingFilePath, "//geneToEntrezIDMapping_",
                                                    disease, "_", bodyRegion, ".csv")


originalGeneFilteredCSV = paste0(geneExpressionPreparationFilePath, "//geneExpressionData_",
                                 disease, "_", bodyRegion, "afterFilteringZeroGenes.csv")

originalGeneFilteredRData = paste0(geneExpressionPreparationFilePath, "//geneExpressionData_",
                                   disease, "_", bodyRegion, "afterFilteringZeroGenes.RData")


updatedOutputAddOn = paste0(finalGeneExpressionPreparationFilePath, "//finalGeneExpressionDataAndDatExpr_",
                            disease, "_", bodyRegion, outputAddOn, ".csv" )

filePathOfGeneExpressionDataUsedByWGCNA = updatedOutputAddOn #updatedOutputAddOnAndDatExpr

updatedOutputAddOnAndDatExpr = paste0(finalGeneExpressionPreparationFilePath, "//finalGeneExpressionDataDatExpr_",
                                      disease, "_", bodyRegion, outputAddOn, ".csv" )


sftThresholdingOutputNameCSV = paste0(softThresholdPowerFilePath, "//softThresholdPowerDF_",netType, "_network_recommendedPower_",disease, "_", bodyRegion, outputAddOn, ".csv")
sftThresholdingOutputNameRData = paste0(softThresholdPowerFilePath, "//softThresholdPowerDF_",netType, "_network_recommendedPower_",disease, "_", bodyRegion, outputAddOn, ".RData")


filePathInfoDF = data.frame(c("geneToEntrezIDPathForGeneExpressionDataRData:"), c(geneToEntrezIDPathForGeneExpressionDataRData),
                            c("RData"), c(paste0("Please note this contains information on the gene to the Entrez IDs mapping for ", folderName)))


colnames(filePathInfoDF) = c("variableName", "filePath", "fileType", "description")


# maybe also have a parameter for overall grouping?  such as soft-threshold power for WGCNA
if (includeTimeStamp == "TRUE"){
  outputAddOnWithDate = paste0(outputAddOn, returnCurrentDateToFilePath())
  updatedOutputAddOn_withDate = paste0(finalGeneExpressionPreparationFilePath, "//finalGeneExpressionDataExpr_",
                                       disease, "_", bodyRegion, outputAddOnWithDate, ".csv" )

  updatedOutputAddOnAndDatExpr_withDate = paste0(finalGeneExpressionPreparationFilePath, "//finalGeneExpressionDatExpr_",
                                                 disease, "_", bodyRegion, outputAddOnWithDate, ".csv" )

  geneExpressionOutputNameCSV_withDate = str_replace_all(updatedOutputAddOn_withDate, " ", "")

}

# please get rid of any possible spaces in the filePath
inputGeneExpressionDataFilePath = str_replace_all(inputGeneExpressionDataFilePath,  " ", "")
geneToEntrezIDPathForGeneExpressionDataCSV = str_replace_all(geneToEntrezIDPathForGeneExpressionDataCSV,  " ", "")
geneToEntrezIDPathForGeneExpressionDataRData = str_replace_all(geneToEntrezIDPathForGeneExpressionDataRData,  " ", "")
geneExpressionOutputNameCSV_withDate = str_replace_all(geneExpressionOutputNameCSV_withDate,  " ", "")
geneExpressionOutputNameCSV = str_replace_all(updatedOutputAddOn, " ", "")
geneExpressionOutputNameRData = str_replace_all(geneExpressionOutputNameCSV, ".csv", ".RData") #paste0("_", diseaseName, "_", bodyRegion, outputAddOn, ".RData")

geneExpressionOutputNameCSV = str_replace_all(geneExpressionOutputNameCSV,  " ", "")
originalGeneFilteredCSV = str_replace_all(originalGeneFilteredCSV,  " ", "")
originalGeneFilteredRData = str_replace_all(originalGeneFilteredRData,  " ", "")
geneExpressionOutputNameRData = str_replace_all(geneExpressionOutputNameRData,  " ", "")
updatedOutputAddOnAndDatExpr_withDate = str_replace_all(updatedOutputAddOnAndDatExpr_withDate,  " ", "")
updatedOutputAddOnAndDatExpr = str_replace_all(updatedOutputAddOnAndDatExpr,  " ", "")
geneExpressionOutputNameRData_withDate = str_replace_all(geneExpressionOutputNameCSV_withDate, ".csv", ".RData") #paste0("_", diseaseName, "_", bodyRegion, outputAddOn, ".RData")
print(paste0("geneExpressionOutputNameCSV_withDate: ", geneExpressionOutputNameCSV_withDate))
print(paste0("geneExpressionOutputNameRData_withDate: ", geneExpressionOutputNameRData_withDate))

print(paste0("geneExpressionOutputNameCSV: ", geneExpressionOutputNameCSV))
print(paste0("geneExpressionOutputNameRData: ", geneExpressionOutputNameRData))


initialRDataFilePath_FilePathDerivationsAndInfoObjectsForWGCNA = paste0(mainOutputFolder, "//initialRDataFilePath_FilePathDerivationsAndInfoObjectsForWGCNA.RData")


# please note that we also save the data objects from WGCNA here so we can call them in future functions without needing to know the power :)
wgcnaWithKmeansAllObjectsExceptTOM_SimpleFileName = paste0(mainOutputFolder, wgcnaToAdd, bodyRegion, "_wgcnaSimplePathALLKeyObjectsExceptTOM",outputAddOn,".RData")



if (includeTimeStamp == "TRUE"){		 
save(fullDiseaseName, folderName, mainOutputFolder,
wgcnaWithKmeansAllObjectsExceptTOM_SimpleFileName,
filePathInfoDF_FileNameRData,
filePathInfoDF_FileNameCSV,
tfsInGeneExpressionDataCSV,
tfsInGeneExpressionDataRData,
wgcnaAndKMeansOutputPath,
wgcnaToAdd,
datTraitsRDataFilePath,
powerRelatedFilePathsRData,
geneExpressionPreparationFilePath,
finalGeneExpressionPreparationFilePath,
entrezMappingFilePath,
softThresholdPowerFilePath,
tomFilePath, wgcnaFilePath, wgcnaWithKMeansAfterFilePath,
pathForWGCNA, MEsFilePath,
moduleAssignmentsFilePath,
geneModuleMembershipFilePath,
allObjectsFilePath, phenotypeAssociationFilePath,
modulePhenoAssociationFilePath,
genePhenoAssociationFilePath,
diseaseTraits, newDate, outputAddOn,
dataScaling, dataScalingOutputMini,
geneToEntrezIDPathForGeneExpressionDataRData,
geneToEntrezIDPathForGeneExpressionDataCSV,
originalGeneFilteredCSV,
originalGeneFilteredRData,
updatedOutputAddOn, filePathOfGeneExpressionDataUsedByWGCNA,
updatedOutputAddOnAndDatExpr,
sftThresholdingOutputNameCSV,
sftThresholdingOutputNameRData,
filePathInfoDF, inputGeneExpressionDataFilePath,
geneExpressionOutputNameCSV_withDate,
geneExpressionOutputNameCSV,
geneExpressionOutputNameRData,
updatedOutputAddOnAndDatExpr_withDate,
geneExpressionOutputNameRData_withDate,
outputAddOnWithDate, updatedOutputAddOn_withDate,
updatedOutputAddOnAndDatExpr_withDate,
geneExpressionOutputNameCSV_withDate,

file = initialRDataFilePath_FilePathDerivationsAndInfoObjectsForWGCNA)
} else {
save(fullDiseaseName, folderName, mainOutputFolder,
wgcnaWithKmeansAllObjectsExceptTOM_SimpleFileName,
filePathInfoDF_FileNameRData,
filePathInfoDF_FileNameCSV,
tfsInGeneExpressionDataCSV,
tfsInGeneExpressionDataRData,
wgcnaAndKMeansOutputPath,
wgcnaToAdd,
datTraitsRDataFilePath,
powerRelatedFilePathsRData,
geneExpressionPreparationFilePath,
finalGeneExpressionPreparationFilePath,
entrezMappingFilePath,
softThresholdPowerFilePath,
tomFilePath, wgcnaFilePath,
wgcnaWithKMeansAfterFilePath,
pathForWGCNA, MEsFilePath,
moduleAssignmentsFilePath,
geneModuleMembershipFilePath,
allObjectsFilePath,
phenotypeAssociationFilePath,
modulePhenoAssociationFilePath,
genePhenoAssociationFilePath,
diseaseTraits, newDate, outputAddOn, dataScaling,
dataScalingOutputMini, geneToEntrezIDPathForGeneExpressionDataRData,
geneToEntrezIDPathForGeneExpressionDataCSV,
originalGeneFilteredCSV,
originalGeneFilteredRData, updatedOutputAddOn,
filePathOfGeneExpressionDataUsedByWGCNA,
updatedOutputAddOnAndDatExpr,
sftThresholdingOutputNameCSV,
sftThresholdingOutputNameRData,
filePathInfoDF, inputGeneExpressionDataFilePath,
geneExpressionOutputNameCSV,
geneExpressionOutputNameRData,
file = initialRDataFilePath_FilePathDerivationsAndInfoObjectsForWGCNA)


}


print(paste0(":) Please note that we wrote a lot of the key files for our WGCNA Gene Co-Expression Network Analysis to an RData File :)"))
print(paste0("initialRDataFilePath_FilePathDerivationsAndInfoObjectsForWGCNA:", initialRDataFilePath_FilePathDerivationsAndInfoObjectsForWGCNA))
return(initialRDataFilePath_FilePathDerivationsAndInfoObjectsForWGCNA)

}





initialRDataFilePath_FilePathDerivationsAndInfoObjectsForWGCNA = pleaseGetInitialFilePathDerivationsAndInfoObjectsForWGCNA(outputPathNameADSNPhenoOutputs,
disease, tissueName, bodyRegion, tfsUsed, performAdditionalKMeansStep, phenotypesFilePath, log2transformInputData, scaleInputData, includeTimeStamp)








######################################






# Step 3:  Annotate modular functions by enrichment analysis


pleaseCreateModuleEnrichmentsFileNames <- function(regionName, pathForWGCNA, performMESHEnrichment, performMolSigDBEnrichment){
print(paste0(":) Please note that this function creates file names and file paths for the Module Enrichments :)"))
library(DOSE)
library(meshes)
library(MeSH.Hsa.eg.db)
library(clusterProfiler)
# "NeuroscienceProjectAxons"
#regionName = folderName

outputNamesList = list() # please create a list to return everything
outputNamesList$moduleEnrichmentsFilePath = moduleEnrichmentsFilePath
moduleEnrichmentsFilePath = paste0(pathForWGCNA, "//ModuleEnrichments")
dir.create(moduleEnrichmentsFilePath, showWarnings = FALSE)


#print(paste0(":) please note that there are ", length(keyModules), " gene modules that are associated with at least 1 of the phenotypes"))


if (performMESHEnrichment){
  print(paste0(":) Please note that since performMESHEnrichment is", performMESHEnrichment, " we are performing Medical Subject Enrichment Analysis! :)"))
  meshEnrichmentsFilePath = paste0(moduleEnrichmentsFilePath, "//MeSH")
  dir.create(meshEnrichmentsFilePath, showWarnings = FALSE)
  outputNamesList$meshEnrichmentsFilePath = meshEnrichmentsFilePath

}
library(clusterProfiler)
library(msigdbr)

if(performMolSigDBEnrichment){
  performMolSigDBEnrichment = TRUE

  molSigDBEnrichmentsFilePath = paste0(moduleEnrichmentsFilePath, "//MolSigDB//")
  dir.create(molSigDBEnrichmentsFilePath, showWarnings = FALSE)
	outputNamesList$molSigDBEnrichmentsFilePath = molSigDBEnrichmentsFilePath

}

return(outputNamesList)



}



disease
tissueName
bodyRegion
outputPathNameADSNPhenoOutputs
haveEnhancersAndPromotersHICData

pleaseCreateChromatinRegulatoryNetworkFileNames <- function(disease, tissueName, bodyRegion, outputPathNameADSNPhenoOutputs, haveEnhancersAndPromotersHICData){
print(paste0(":) Please note that this function creates file names and file paths for the Chromatin (Promoter and Enhancer) Regulatory Network (Epigenomics) :)"))

  fullDiseaseName = paste0(disease, "Disease")
  folderName = str_replace_all(paste0(disease, "_", tissueName, "_", bodyRegion), " ", "")
  dir.create(folderName)
  print(paste0("outputPathNameADSNPhenoOutputs: ", outputPathNameADSNPhenoOutputs))

  mainOutputFolder = paste0(
    outputPathNameADSNPhenoOutputs, folderName)
  dir.create(mainOutputFolder, showWarnings = FALSE)
  
  
  
  
  

# 5.  Please predict gene regulatory networks for genes and modules

enhancerAndPromoterInteraction = paste0(mainOutputFolder, "//PromAndEnhInteractions//")
dir.create(enhancerAndPromoterInteraction, showWarnings = FALSE)



chromatinInteraction_Step2_FilePath = paste0(enhancerAndPromoterInteraction, "Chromatin_scgrnom_getTF_Step2//")
chromatinInteraction_Step2_FilePath_InitialRDataObjects = paste0(enhancerAndPromoterInteraction, "Chromatin_scgrnom_getTF_Step2//getTF_InitialRData//")
chromatinInteraction_Step2_FilePath_CleanedRDataObjects = paste0(enhancerAndPromoterInteraction, "Chromatin_scgrnom_getTF_Step2//getTF_FinalRData//")
chromatinInteraction_Step2_FilePath_CleanedCSVObjects = paste0(enhancerAndPromoterInteraction, "Chromatin_scgrnom_getTF_Step2//getTF_finalCSVData//")

#chromatinInteraction_Step2_FilePath_CSVWithEntrezID =

dir.create(chromatinInteraction_Step2_FilePath)
dir.create(chromatinInteraction_Step2_FilePath_InitialRDataObjects)
dir.create(chromatinInteraction_Step2_FilePath_CleanedRDataObjects)
dir.create(chromatinInteraction_Step2_FilePath_CleanedCSVObjects)

chromatinInteractionFilesRData = paste0(enhancerAndPromoterInteraction, folderName,"_chromatinInteractionFilePathNames.RData") 

if (haveEnhancersAndPromotersHICData){
  chromatinInteraction_Step1_FilePath = paste0(enhancerAndPromoterInteraction, "Chromatin_scgrnom_Step1//")
  dir.create(chromatinInteraction_Step1_FilePath)
  
  save(fullDiseaseName, mainOutputFolder, enhancerAndPromoterInteraction, 
chromatinInteraction_Step2_FilePath, 
chromatinInteraction_Step2_FilePath_InitialRDataObjects, 
chromatinInteraction_Step2_FilePath_CleanedRDataObjects,
chromatinInteraction_Step2_FilePath_CleanedCSVObjects,
chromatinInteraction_Step1_FilePath,
file = chromatinInteractionFilesRData)
} else { #if (haveEnhancersInteractionData}{
  chromatinInteraction_Step1ADSNPheno_FilePath = paste0(enhancerAndPromoterInteraction, "Chromatin_adsnpheno_Step1//")
  dir.create(chromatinInteraction_Step1ADSNPheno_FilePath)
   save(fullDiseaseName, mainOutputFolder, enhancerAndPromoterInteraction, 
chromatinInteraction_Step2_FilePath, 
chromatinInteraction_Step2_FilePath_InitialRDataObjects, 
chromatinInteraction_Step2_FilePath_CleanedRDataObjects,
chromatinInteraction_Step2_FilePath_CleanedCSVObjects,
chromatinInteraction_Step1ADSNPheno_FilePath,
file = chromatinInteractionFilesRData)
  
}


print(paste0(":) Please note that the information on the Chromatin Interaction Regulatory Network File Path Names is stored in this RData Object: chromatinInteractionFilesRData = ", chromatinInteractionFilesRData))
return(chromatinInteractionFilesRData)
}







pleaseCreateGeneExpressionRegulatoryNetworkFileNames <- function(disease, tissueName, bodyRegion, outputPathNameADSNPhenoOutputs, dataScalingOutputMini, outputAddOn, tfsUsed, weightThresholdGenie3){
print(paste0(":) Please note that this function creates file names and file paths for the Gene Expression Regulatory Network (RTN, GENIE3, TReNA Ensemble) :)"))

  fullDiseaseName = paste0(disease, "Disease")
  folderName = str_replace_all(paste0(disease, "_", tissueName, "_", bodyRegion), " ", "")
  dir.create(folderName)
  print(paste0("outputPathNameADSNPhenoOutputs: ", outputPathNameADSNPhenoOutputs))

  mainOutputFolder = paste0(
    outputPathNameADSNPhenoOutputs, folderName)
  dir.create(mainOutputFolder, showWarnings = FALSE)
  
  

transcriptionalRegulatoryNetworkOutputPath = paste0(mainOutputFolder, "//", "GeneRegulatNetwork")
filePathOfRTN  = paste0(transcriptionalRegulatoryNetworkOutputPath, "//", "RTN")
filePathOfGenie3  = paste0(transcriptionalRegulatoryNetworkOutputPath, "//", "Genie3")
filePathOfTrena  = paste0(transcriptionalRegulatoryNetworkOutputPath, "//", "TrenaEnsemble")
dir.create(transcriptionalRegulatoryNetworkOutputPath)
dir.create(filePathOfRTN)
dir.create(filePathOfGenie3)
dir.create(filePathOfTrena)

# please note that this file path will combine the RTN, Genie3, Trena GRNs from gene expression data along with Trrust2 database results
filePathOfAll4GRNs  = paste0(transcriptionalRegulatoryNetworkOutputPath, "//", "CombinedGRN_4Sources//")

if (minNumSourcesGeneGRN > 1){
        grnResultsCombinedFileNameCSV = paste0(filePathOfAll4GRNs,"finalComboGeneExpRegNet_GenieRTNTrenaTrrust2_TFtoGene_", disease, "_", bodyRegion, "_",
		dataScalingOutputMini, "_", (minNumSourcesGeneGRN), "_minSources.csv")

		 grnResultsCombinedFileNameRData = paste0(filePathOfAll4GRNs,"finalComboGeneExpRegNet_GenieRTNTrenaTrrust2_TFtoGene_", disease, "_", bodyRegion, "_",
		dataScalingOutputMini, "_", (minNumSourcesGeneGRN), "_minSources.RData")
} else {
		grnResultsCombinedFileNameCSV = paste0(filePathOfAll4GRNs,"finalComboGeneExpRegNet_GenieRTNTrenaTrrust2_TFtoGene_", disease, "_", bodyRegion, "_",
		dataScalingOutputMini, ".csv")

		 grnResultsCombinedFileNameRData = paste0(filePathOfAll4GRNs,"finalComboGeneExpRegNet_GenieRTNTrenaTrrust2_TFtoGene_", disease, "_", bodyRegion, "_",
		dataScalingOutputMini, ".RData")
}


grnResultsComboForPythonInputFileNameCSV = grnResultsCombinedFileNameCSV

dir.create(filePathOfAll4GRNs)

newDate = str_replace_all(Sys.Date(), "-", "_")
csvPart = paste0("RTN_", tfsUsed, "TFs_", disease, "_", bodyRegion,  outputAddOn,  newDate, "Final.csv")
rDataPart = paste0("RTN_", tfsUsed, "TFs_", disease, "_", bodyRegion, outputAddOn, "_", newDate, "Final.RData")
oldRDataPart = paste0(filePathOfRTN, "RTN_", tfsUsed, "TFs_", disease, "_", bodyRegion, outputAddOn, "_", newDate, "_Original.RData")
RTN_rDataPart = paste0(filePathOfRTN, "//RTN_", tfsUsed, "TFs_", disease, "_", bodyRegion, outputAddOn,  "_FinalDataObjects.RData")

csvPartRTNFinal = paste0(filePathOfRTN, "//RTN_", tfsUsed, "TFs_", disease, "_", bodyRegion,  outputAddOn,  newDate, "FinalRTN_rtnRegulonsDF.csv")
csvPartRTNFinal_NoDate = paste0(filePathOfRTN, "//RTN_", tfsUsed, "TFs_", disease, "_", bodyRegion,  outputAddOn, "FinalRTN_rtnRegulonsDF.csv")


newCsvPart = paste0(filePathOfRTN, "//", csvPart)
newRDataPart = paste0(filePathOfRTN,"//",  rDataPart)
print(paste0("RTN newCsvPart: ", newCsvPart))
print(paste0("RTN rDataPart: ", newRDataPart))
rm(csvPart)
rm(rDataPart)
print(newCsvPart)
print(newRDataPart)

genie3CsvPart = paste0(filePathOfGenie3,"//", tfsUsed, "TFs_", disease, "_", bodyRegion, "_", dataScaling, "_", newDate, "Final.csv")
genie3RDataPart = paste0(filePathOfGenie3,"//", tfsUsed, "TFs_", disease, "_", bodyRegion, "_", dataScaling,  "_", newDate, "Final.RData")
oldGenie3RDataPart = paste0(filePathOfGenie3, "//",tfsUsed, "TFs_", disease, "_", bodyRegion, "_", dataScaling, "_", newDate, "_Original.RData")
genie3_rDataPart = paste0(filePathOfGenie3, "//",tfsUsed, "TFs_", disease, "_", bodyRegion, "_", dataScaling, "_FinalDataObjects_Genie3All.RData")

genie3RDataPartFinal = paste0(filePathOfGenie3, "//",tfsUsed, "TFs_", disease, "_", bodyRegion, "_", dataScaling,"_minWeightThresh_", weightThresholdGenie3, "_FinalGenie3DF.RData")
genie3CsvPartFinal = paste0(filePathOfGenie3, "//",tfsUsed, "TFs_", disease, "_", bodyRegion, "_", dataScaling, "_minWeightThresh_", weightThresholdGenie3, "_FinalGenie3DF.csv")


filePathOfTFtoModule  = paste0(transcriptionalRegulatoryNetworkOutputPath, "//", "TFtoModuleRTN//")
dir.create(filePathOfTFtoModule)

trenaOutpathALLCSV = paste0(filePathOfTrena, "//trenaEnsemble_", bodyRegion, "_ALL_",
                            disease, "_", tissueName, ".csv")


trenaOutpathALLRData = paste0(filePathOfTrena, "//trenaEnsemble_", bodyRegion, "_ALL_",
                              disease, "_", tissueName, ".RData")


trenaOutpathALLCSV_Final = paste0(filePathOfTrena, "//trenaEnsemble_", bodyRegion, "_ALL_",
                                  disease, "_", tissueName, "_final_trenaEnsembleDF.csv")


trenaOutpathALLRData_Final = paste0(filePathOfTrena, "//trenaEnsemble_", bodyRegion, "_ALL_",
                                    disease, "_", tissueName, "_final_trenaEnsembleDF.RData")
print(paste0("trenaOutpathALLCSV: ", trenaOutpathALLCSV))
print(paste0("trenaOutpathALLRData: ", trenaOutpathALLRData))

finalFullGeneRegulatoryNetworkOutputPath = paste0(mainOutputFolder, "//", "FinalFullGRN//")

fullNetALLOrganizedCSV = paste0(finalFullGeneRegulatoryNetworkOutputPath,
                                "fullNetwork_chromAndGeneExpr_TFs_ALL_", folderName, "_FINAL.csv")
fullNetALLOrganizedCSV_Done = paste0(finalFullGeneRegulatoryNetworkOutputPath,
                                     "final_fullNetwork_chromAndGeneExpr_TFs_ALL_", folderName, "_FINAL.csv")



finalFullNetworkChromatinAndGeneExpressCSV = paste0(finalFullGeneRegulatoryNetworkOutputPath,
                                     "FINAL_FULL_GRN_chromAndGeneExpr_TFs_ALL_", folderName, ".csv")



geneRegulatoryPathwayFileNamesRData = paste0(finalFullGeneRegulatoryNetworkOutputPath, "finalFullGRNRDataObjectsForFilePaths.RData")

save(finalFullNetworkChromatinAndGeneExpressCSV,
fullNetALLOrganizedCSV_Done,
fullNetALLOrganizedCSV,
finalFullGeneRegulatoryNetworkOutputPath,
trenaOutpathALLCSV,
trenaOutpathALLRData,
filePathOfTrena,
trenaOutpathALLCSV_Final,
filePathOfTFtoModule,
transcriptionalRegulatoryNetworkOutputPath,
genie3CsvPartFinal,
genie3RDataPartFinal,
genie3_rDataPart,
oldGenie3RDataPart,
genie3CsvPart,
genie3RDataPart,
newRDataPart,
newCsvPart,
csvPartRTNFinal_NoDate,
csvPartRTNFinal,
RTN_rDataPart,
oldRDataPart,
rDataPart,
csvPart,
filePathOfAll4GRNs,
grnResultsCombinedFileNameCSV,
grnResultsComboForPythonInputFileNameCSV,
grnResultsCombinedFileNameRData,
filePathOfGenie3,
filePathOfRTN, file = geneRegulatoryPathwayFileNamesRData)

print(paste0(":) Please note that the information on the Gene Expression Regulatory Network File Path Names is stored in this RData Object: geneRegulatoryPathwayFileNamesRData = ", geneRegulatoryPathwayFileNamesRData))
return(geneRegulatoryPathwayFileNamesRData)

}
######################
# SNPs GWAS:

pleaseCreateSNPsMotifBreakrFileNames <- function(outputPathNameADSNPhenoOutputs, pValThreshForSNPs, disease, tissueName, bodyRegion){

print(paste0(":) Please note that this function creates file names and file paths for the SNP DF and the MotifBreakR Results :)"))

fullDiseaseName = paste0(disease, "Disease")
pValThreshForSNPs
folderName = str_replace_all(paste0(disease, "_", tissueName, "_", bodyRegion), " ", "")
dir.create(folderName)
print(paste0("outputPathNameADSNPhenoOutputs: ", outputPathNameADSNPhenoOutputs))

mainOutputFolder = paste0(
  outputPathNameADSNPhenoOutputs, folderName)
dir.create(mainOutputFolder, showWarnings = FALSE)

print(paste0(":) Please note that a directory (folder) has been created here: ",
             mainOutputFolder, " for all of the output of ADSNPheno"))
			 
			 
			 


snpsOutputPath = paste0(mainOutputFolder, "//", "SNPs_Analysis")
filePathOfGWAS  = paste0(snpsOutputPath, "//", "GWAS_Data")
filePathOfMotifBreakR_Step1  = paste0(snpsOutputPath, "//", "MotifbreakR_Part1//")
filePathOfMotifBreakR_Step2  = paste0(snpsOutputPath, "//", "MotifbreakR_Part2_BrokenTFBS//")
filePathOfMotifBreakR_Final  = paste0(snpsOutputPath, "//", "MotifbreakR_FinalResults//")

pValName = str_replace(pValThreshForSNPs, "-", "Minus")

updatedGWASDataSetFileNameALL_CSV = paste0(filePathOfGWAS, "//", "allGWAS_Data_", fullDiseaseName, "_SNPsWithPvalBelow_", pValName,".csv")
updatedGWASDataSetFileNameSNPandP_CSV = paste0(filePathOfGWAS, "//", "SNPsAndP_GWAS_", fullDiseaseName, "_SNPsWithPBelow_", pValName,".csv")
updatedGWASDataSetFileNameALL_RData = paste0(filePathOfGWAS, "//", "allGWAS_Data_", fullDiseaseName, "_SNPsWithPvalBelow_", pValName,".RData")
updatedGWASDataSetFileNameSNPandP_RData = paste0(filePathOfGWAS, "//", "SNPsAndP_GWAS_", fullDiseaseName, "_SNPsWithPBelow_", pValName,".RData")


snpsOutputPathRDataObjects = paste0(snpsOutputPath, "//", "SNPs_Analysis_RDataObjectsFilePaths.RData")

save(fullDiseaseName, mainOutputFolder, snpsOutputPath, filePathOfGWAS, filePathOfMotifBreakR_Step1,
 filePathOfMotifBreakR_Step2, filePathOfMotifBreakR_Final, pValName, 
 updatedGWASDataSetFileNameALL_CSV, updatedGWASDataSetFileNameSNPandP_CSV, 
 updatedGWASDataSetFileNameALL_RData, updatedGWASDataSetFileNameSNPandP_RData,
 file = snpsOutputPathRDataObjects)

print(paste0(":) Please note that the SNPs and MotifbreakR file paths info (file path names) is stored here: ", snpsOutputPathRDataObjects))
return(snpsOutputPathRDataObjects)

}







#filePathOfTrena  = paste0(transcriptionalRegulatoryNetworkOutputPath, "//", "TrenaEnsemble")
dir.create(snpsOutputPath)
dir.create(filePathOfGWAS)
dir.create(filePathOfMotifBreakR_Step1)
dir.create(filePathOfMotifBreakR_Step2)
dir.create(filePathOfMotifBreakR_Final)


snpDFPath = paste0(filePathOfMotifBreakR_Final, "motifbreakr_FinalResultsForSNPsWithP_", pValName,".csv")
snpDFPathRData = paste0(filePathOfMotifBreakR_Final, "motifbreakr_FinalResultsForSNPsWithP_", pValName,".RData")

#################
# SNP Regulatory Network (Sub-Networks:)
snpRegulatNetworkFolder = paste0(mainOutputFolder, "//SNPRegulatNetwork//")
#dir.create(snpRegulatNetworkFolder)
#\\snpRegulatoryPathway_"
SNPOutputNameStem = paste0(snpRegulatNetworkFolder, "snpRegulatoryPathway_")
# SNPOutputNameStem = "F:\\organizedAlzheimers\\ADSNPhenoOutputs_NEW\\Alzheimers_Brain_Hippocampus\\SNPRegulatNetwork\\snpRegulatoryPathway_"
