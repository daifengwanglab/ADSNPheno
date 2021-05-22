# :) Please note that this R file by Saniya has all the parameters we need

tissueName = "Brain"
disease = "Alzheimers"
alzheimerTraits = read.csv("D://organizedAlzheimers//GeneExpressionPreparation//LateralTemporalLobePhenotypesUpdated.csv",
                           header = TRUE)
parentPath = "D://organizedAlzheimers//GeneExpressionPreparation//"
outputPath = "D://organizedAlzheimers//GeneExpressionPreparation//outputs"
minModuleSize = 30;
maxNumPowers = 50
min.genes.for.grey=0
nIterations = 201 
meg = 0
netType="signed"
numRoundingDigits = 3
pValueCutOffForSignificance = 0.05
includeTimeStamp = "TRUE"
bodyRegion = "HippocampusCa1"
useRecommendedPower = "TRUE"
ourPower = NULL 
log2transformInputData = "TRUE" # (should we apply a log2(x+1) transform on data
scaleInputData = "TRUE" #TRUE" #TRUE" # should we apply a scale on data
saveOutputs = "TRUE"
positiveCorrelationsOnly = TRUE





powerVal = "Unknown"
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



tfsUsed = "LambertAndJaspar"

# WGCNAwithSpecificPower
# WGCNA with recommended power

#transcriptionalFilePath = 

newDate = str_replace_all(Sys.Date(), "-", "_")
csvPart = paste0(tfsUsed, "TFs_", disease, "_", bodyRegion, "_", dataScaling, "_power", powerVal, "_", newDate, "Final.csv")
rDataPart = paste0(tfsUsed, "TFs_", disease, "_", bodyRegion, "_", dataScaling, "_power", powerVal, "_", newDate, "Final.RData")
oldRDataPart = paste0(filePathOfRTN, tfsUsed, "TFs_", disease, "_", bodyRegion, "_", dataScaling, "_power", powerVal, "_", newDate, "_Original.RData")
RTN_rDataPart = paste0(tfsUsed, "TFs_", disease, "_", bodyRegion, "_", dataScaling, "_power", powerVal, "_FinalDataObjects.RData")



genie3CsvPart = paste0(filePathOfGenie3, tfsUsed, "TFs_", disease, "_", bodyRegion, "_", dataScaling, "_power", powerVal, "_", newDate, "Final.csv")
genie3RDataPart = paste0(filePathOfGenie3, tfsUsed, "TFs_", disease, "_", bodyRegion, "_", dataScaling, "_power", powerVal, "_", newDate, "Final.RData")
oldGenie3RDataPart = paste0(filePathOfGenie3, tfsUsed, "TFs_", disease, "_", bodyRegion, "_", dataScaling, "_power", powerVal, "_", newDate, "_Original.RData")
genie3_rDataPart = paste0(filePathOfGenie3, tfsUsed, "TFs_", disease, "_", bodyRegion, "_", dataScaling, "_power", powerVal, "_FinalDataObjects.RData")

# printing out a text file or csv with information on location of each element we have,
# along with what it means

# SNP threshold, etc.
# key KEGG Pathways to visualize
# WGCNA userListEnrichment





bodyRegion = "Lateral Temporal Lobe"
tissueName = "Brain"
disease = "Alzheimers"
powerVal = "10"
dataScaling = "log2TransformedOnly" #  <-- will come from WGCNA
tfsUsed = "LambertAndJaspar"
