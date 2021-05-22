# :) Please note that this R file by Saniya has all the parameters we need
source("D:\\organizedAlzheimers\\Setup\\packagesNeeded.R")
source("D:\\organizedAlzheimers\\Setup\\parametersNeededByUser.R")
source("D:\\organizedAlzheimers\\Setup\\functionsNeeded.R")
diseaseTraits = read.csv(phenotypesFilePath,
                           header = TRUE)


tfsDF <- read.csv(transcriptionFactorsFilePath, 
                  header = TRUE)[,c(1,2)]
tfs = tfsDF[,1]


infoDF = read.csv(geneToEntrezIDMappingPath, header = TRUE)
dim(infoDF)
rowsToRemove = c(which(infoDF[,2] == "Withdrawn"), which(infoDF[,2] == "unknown"))
infoDF = unique(infoDF[-rowsToRemove,])
colnames(infoDF) = c("geneSymbolName", "entrezID")
dim(infoDF)
entrezIDMappingVec = as.vector(infoDF[,2])
names(entrezIDMappingVec) = infoDF[,1]

#entrezIDMappingHashTable = hash(infoDF)
entrezIDMappingHashTable = hash() #infoDF)
#.set(entrezIDMappingHashTable, keys = infoDF[,1], values = infoDF[,2])
#h <- hash()
#.set( h, keys=letters, values=1:26 )
entrezIDMappingHashTable[["ALAD"]]

newDate = str_replace_all(Sys.Date(), "-", "_")

# entrezIDMappingHashTable = hash(as.vector(infoDF[,1]),
#                                 as.vector(infoDF[,2]))
numOfRoundingDigits = numberOfRoundingDigits
for (i in 1:length(bodyRegionVec)){
  mainDir = wgcnaAndKMeansOutputPath
  subDir = bodyRegionVec[i]
  outputPath = file.path(mainDir, subDir)
  print(paste0(":) Please note the WGCNA output for ", subDir, 
               " will be written to this folder: ", 
               outputPath))
  dir.create(outputPath, showWarnings = FALSE)
  dir.create(file.path(mainDir, subDir, "CoexpressionNetwork_WGCNAAndKmeans"), showWarnings = FALSE)
}

wgcnaAndKMeansPathsForOutputs = paste0(mainDir, "//", bodyRegion)
print(paste0(":) please note that the outputs for ", bodyRegion, 
             " will be written here: ", wgcnaAndKMeansPathsForOutputs))



powerVal = powerEstimate# "Unknown"
#dataScaling = "log2TransformedAndScaled" #  <-- will come from WGCNA
fullDiseaseName = paste0(disease, "Disease")

######### Please note the function
if (log2transformInputData == "TRUE"){
  if (scaleInputData == "TRUE"){
    outputAddOn = "_log2TransformedAndThenScaledInput"
    print(":) Please note we will apply a Log2(x+1) transformation on the input gene expression data and then will scale that.")
    print(":) This final data set will be our input for WGCNA")
    dataScaling = "log2TransformedAndScaledGeneExpression" #  <-- will come from WGCNA
    
  } else {
    outputAddOn = "_log2TransformedInputOnly"
    print(":) Please note we will apply a Log2(x+1) transformation on the input gene expression data.")
    print(":) This final data set will be our input for WGCNA")
    dataScaling = "log2TransformedGeneExpressionOnly" #  <-- will come from WGCNA
    
  }
} else if (scaleInputData == "TRUE"){
  outputAddOn = "_ScaledInputOnly"
  print(":) Please note we will apply a Scale transformation on the input gene expression data.")
  print(":) This final data set will be our input for WGCNA")
  dataScaling = "ScaledGeneExpressionOnly" #  <-- will come from WGCNA
  
} else {
  outputAddOn = "_OriginalInputData"
  print(":) Please note we will use our original input gene expression dataset for WGCNA.")
  dataScaling = "originalGeneExpression" #  <-- will come from WGCNA
  
}


if (includeTimeStamp == "TRUE"){
  outputAddOnWithDate = paste0(outputAddOn, returnCurrentDateToFilePath())
}


geneToEntrezIDPathForGeneExpressionDataRData = paste0(wgcnaAndKMeansPathsForOutputs, "//geneToEntrezIDMapping_", 
                                                    fullDiseaseName, "_", bodyRegion, ".RData")

geneToEntrezIDPathForGeneExpressionDataCSV = paste0(wgcnaAndKMeansPathsForOutputs, "//geneToEntrezIDMapping_", 
                                                 fullDiseaseName, "_", bodyRegion, ".csv")

originalGeneFilteredCSV = paste0(wgcnaAndKMeansPathsForOutputs, "//geneExpressionData_", 
                                     fullDiseaseName, "_", bodyRegion, "afterFilteringZeroGenes.csv")

originalGeneFilteredRData = paste0(wgcnaAndKMeansPathsForOutputs, "//geneExpressionData_", 
                                 fullDiseaseName, "_", bodyRegion, "afterFilteringZeroGenes.RData")

updatedOutputAddOn_withDate = paste0(wgcnaAndKMeansPathsForOutputs, "//finalGeneExpressionDataAndDatExpr_", 
                     fullDiseaseName, "_", bodyRegion, outputAddOnWithDate, ".csv" )

updatedOutputAddOn = paste0(wgcnaAndKMeansPathsForOutputs, "//finalGeneExpressionDataAndDatExpr_", 
                                     fullDiseaseName, "_", bodyRegion, outputAddOn, ".csv" )


updatedOutputAddOnAndDatExpr_withDate = paste0(wgcnaAndKMeansPathsForOutputs, "//finalGeneExpressionDataDatExpr_", 
                                     fullDiseaseName, "_", bodyRegion, outputAddOnWithDate, ".csv" )

updatedOutputAddOnAndDatExpr = paste0(wgcnaAndKMeansPathsForOutputs, "//finalGeneExpressionDataDatExpr_", 
                            fullDiseaseName, "_", bodyRegion, outputAddOn, ".csv" )

geneExpressionOutputNameCSV_withDate = str_replace_all(updatedOutputAddOn_withDate, " ", "")



inputGeneExpressionDataFilePath = str_replace_all(inputGeneExpressionDataFilePath,  " ", "")
geneToEntrezIDPathForGeneExpressionDataCSV = str_replace_all(geneToEntrezIDPathForGeneExpressionDataCSV,  " ", "")
geneToEntrezIDPathForGeneExpressionDataRData = str_replace_all(geneToEntrezIDPathForGeneExpressionDataRData,  " ", "")
geneExpressionOutputNameCSV_withDate = str_replace_all(geneExpressionOutputNameCSV_withDate,  " ", "")
geneExpressionOutputNameCSV = str_replace_all(geneExpressionOutputNameCSV,  " ", "")
originalGeneFilteredCSV = str_replace_all(originalGeneFilteredCSV,  " ", "")
originalGeneFilteredRData = str_replace_all(originalGeneFilteredRData,  " ", "")
geneExpressionOutputNameRData = str_replace_all(geneExpressionOutputNameRData,  " ", "")
#geneExpressionOutputNameRData_withDate = str_replace_all(geneExpressionOutputNameRData_withDate,  " ", "")
updatedOutputAddOnAndDatExpr_withDate = str_replace_all(updatedOutputAddOnAndDatExpr_withDate,  " ", "")
updatedOutputAddOnAndDatExpr = str_replace_all(updatedOutputAddOnAndDatExpr,  " ", "")


geneExpressionOutputNameRData_withDate = str_replace_all(geneExpressionOutputNameCSV_withDate, ".csv", ".RData") #paste0("_", diseaseName, "_", bodyRegion, outputAddOn, ".RData")
print(paste0("geneExpressionOutputNameCSV_withDate: ", geneExpressionOutputNameCSV_withDate))
print(paste0("geneExpressionOutputNameRData_withDate: ", geneExpressionOutputNameRData_withDate))


geneExpressionOutputNameCSV = str_replace_all(updatedOutputAddOn, " ", "")

geneExpressionOutputNameRData = str_replace_all(geneExpressionOutputNameCSV, ".csv", ".RData") #paste0("_", diseaseName, "_", bodyRegion, outputAddOn, ".RData")
print(paste0("geneExpressionOutputNameCSV: ", geneExpressionOutputNameCSV))
print(paste0("geneExpressionOutputNameRData: ", geneExpressionOutputNameRData))


sftThresholdingOutputNameCSV = paste0(wgcnaAndKMeansPathsForOutputs, "//softThresholdPowerDF_",netType, "_network_recommendedPower_",fullDiseaseName, "_", bodyRegion, outputAddOn, ".csv")
sftThresholdingOutputNameRData = paste0(wgcnaAndKMeansPathsForOutputs, "//softThresholdPowerDF_",netType, "_network_recommendedPower_",fullDiseaseName, "_", bodyRegion, outputAddOn, ".RData")

