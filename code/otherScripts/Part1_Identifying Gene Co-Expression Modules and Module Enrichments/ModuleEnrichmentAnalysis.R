#########################################################################################################################################
### INFORMATION:  Please note that the code below was written by Saniya Khullar and Daifeng Wang, Ph.D.
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

phenotypesFilePath = paste0(dataFolder, "phenotypes//AlzheimersLateralTemporalLobePhenotypesUpdated.csv") 

log2transformInputData = "TRUE" # (should we apply a log2(x+1) transform on data
scaleInputData = "FALSE"  # should we apply a scale() on data

includeTimeStamp = "TRUE"
netType = "signed"  # # network type:but please check WGCNA for additional parameters.  We retained the defaults

####### GENE MODULE ANNOTATION:
# Step 3:  Annotate modular functions by enrichment analysis
enrichmentsForKeyModulesOnly = FALSE # please note that this says if we only want the enrichments for modules that are associated with at least 1 key phenotype.  If positiveCorrelationsOnly is TRUE then this is only for modules positively associated with at least 1  phenotype

performMESHEnrichment = TRUE
MESHCategoryVec = c("A", "B", "C", "D", "E", "F", "G", "H")
meshDatabaseNamesVec = c("gendoo", "gene2pubmed", "RBBH")

performMolSigDBEnrichment = TRUE
molSigDBCategoryVec = c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "H")


keggPathwaysOfInterest <- c("hsa05010", "hsa05171") # what are the IDs of the kegg pathways of interest?


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

# Step 2:  How to construct co-expression network and modules (WGCNA and K-means too).  Define parameters

newDate = str_replace_all(Sys.Date(), "-", "_")



# Step 3:  Annotate modular functions by enrichment analysis
load(powerRelatedFilePathsRData)
load(wgcnaWithKmeansAllObjectsExceptTOM)

powerVal = powerEstimate

if (enrichmentsForKeyModulesOnly){
  # we only load data for the key modules then
  load(keyGeneAssignmentsWithKmeansOutputNameRData)
  #clusterInfo = clusterInfo[which(clusterInfo[,colNum] %in% keyModules),]

} else {
  # we load info on genes in all modules
  load(geneAssignmentsWithKmeansOutputNameRData)
}
head(clusterInfo)
colnames(clusterInfo)
dim(clusterInfo)
entrezIDCol = grep("entrez", colnames(clusterInfo))
print(paste0(":) Please note the entrezIDCol is: ", entrezIDCol))
modColNum = grep("module", colnames(clusterInfo)) # module column number

print(paste0(":) please note that there are ", length(keyModules), " gene modules that are associated with at least 1 of the phenotypes"))

moduleEnrichmentsInfoList = pleaseCreateModuleEnrichmentsFileNames(regionName, pathForWGCNA, performMESHEnrichment, performMolSigDBEnrichment)
if (performMESHEnrichment){
  print(paste0(":) Please note that since performMESHEnrichment is", performMESHEnrichment, " we are performing Medical Subject Enrichment Analysis! :)"))
  meshEnrichmentsFilePath =  moduleEnrichmentsInfoList$meshEnrichmentsFilePath

  meshCategoriesList = runAllMeSHEnrichments(categoryVec = MESHCategoryVec, #c("A", "B", "C", "D", "E", "F", "G", "H"),
                                             powerVal = powerEstimate,
                                             region = regionName,
                                             geneIdCol = entrezIDCol,
                                             geneClusterDFToUse = clusterInfo,
                                             parentName = meshEnrichmentsFilePath,
                                             colNum = modColNum,
                                             meshDatabaseNames = meshDatabaseNamesVec)
  print(paste0(":) Please note that we are done with Medical Subject Headings Enrichments!"))
}

colNum = modColNum

if(performMolSigDBEnrichment){
  #performMolSigDBEnrichment = TRUE
  molSigDBEnrichmentsFilePath = moduleEnrichmentsInfoList$molSigDBEnrichmentsFilePath
  print(paste0(":) Please note that since performMolSigDBEnrichment is", performMolSigDBEnrichment, " we are performing Molecular Signatures Database Enrichment Analysis! :)"))

  molsigDBDF = data.frame()
  for (i in 1:length(molSigDBCategoryVec)){
    categoryMolSigDB = molSigDBCategoryVec[i]
    print(paste0(":) Please note that we are on i = ", i, " (category = ", categoryMolSigDB,
                 " of the Molecular Signatures Database :)"))
    molSidDFForCategory = molSigDBFunctionsEnrichment(categoryName = categoryMolSigDB, powerVal = powerEstimate,
                                                      geneIdCol = entrezIDCol, colNum,
                                                      region = regionName,
                                                      clusterInfo,
                                                      parentName = molSigDBEnrichmentsFilePath,
                                                      kMeans = performAdditionalKMeansStep)

    molsigDBDF = unique(rbind(molsigDBDF, molSidDFForCategory))
  }
  print(paste0(":) Please note that we are done with Molecular Signatures Database (MolSigDB) Enrichments!"))

}


