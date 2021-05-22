# :) Please note that this R file by Saniya has all the parameters we need

bodyRegion = "Hippocampus"

inputGeneExpressionDataFilePath = "F://organizedAlzheimers//Setup//InputDataFiles//originalGeneExpressionDataSets//originalHippocampalCa1RegionGeneExpressionData.csv"
phenotypesFilePath = "F://organizedAlzheimers//Setup//InputDataFiles//phenotypes//AlzheimersHippocampusCa1Phenotypes.csv"
parentPath = "F://organizedAlzheimers//GeneExpressionPreparation//"
outputPath = "F://organizedAlzheimers//GeneExpressionPreparation//outputs"
log2transformInputData = "TRUE" # (should we apply a log2(x+1) transform on data
scaleInputData = "TRUE" # should we apply a scale on data

print(paste0(":) please note body region = ", bodyRegion))
wgcnaAndKMeansOutputPath = "F://organizedAlzheimers//GeneExpressionPreparation//outputs"
geneToEntrezIDMappingPath = "F://organizedAlzheimers//entrezMappingInfo//infoDFMappingsGeneSymbolAndID.csv"
transcriptionFactorsFilePath = "F:\\organizedAlzheimers\\TranscriptionalRegulator\\inputs\\LambertAndJasparTFs.csv"
percentilesForCoexpressionNetworkVec = c(0, 0.25, 0.50, 0.75, 0.85, 0.90, 0.95, 1)
minTOMThresholdValue = 0.0001
print(paste0(":) Please note that the body region ye have selected is: ", bodyRegion))
tissueName = "Brain"
disease = "Alzheimers" #Alzheimers"
minModuleSize = 30;
maxNumPowers = 50
min.genes.for.grey=0
nIterations = 201 
meg = 0
netType="signed"
numberOfRoundingDigits = 3
pValueCutOffForSignificance = 0.05
includeTimeStamp = "TRUE"
useRecommendedPower = "TRUE"
ourPower = NULL 
saveOutputs = "TRUE"
positiveCorrelationsOnly = FALSE #TRUE

tfsUsed = "LambertAndJaspar"
powerVal = "Unknown"
#fullDiseaseName = paste0(disease, "Disease")
# 
# if (bodyRegion == "HippocampusCa1"){
#   inputGeneExpressionDataFilePath = "F://organizedAlzheimers//Setup//InputDataFiles//originalGeneExpressionDataSets//originalHippocampalCa1RegionGeneExpressionData.csv"
#   #inputGeneExpressionDataFilePath = "F://organizedAlzheimers//DataSets//GeneExpressionDataSet//finalizedHippocampalCa1RegionData.csv"
#   
#   
#   phenotypesFilePath = "F://organizedAlzheimers//Setup//InputDataFiles//phenotypes//AlzheimersHippocampusCa1Phenotypes.csv"
#   #phenotypesFilePath = "F://organizedAlzheimers//DataSets//PhenotypesAndTraits//AlzheimersTraits.csv"
#   
#   parentPath = "F://organizedAlzheimers//GeneExpressionPreparation//"
#   outputPath = "F://organizedAlzheimers//GeneExpressionPreparation//outputs"
#   
#   log2transformInputData = "TRUE" # (should we apply a log2(x+1) transform on data
#   scaleInputData = "TRUE" # should we apply a scale on data
# 
# 
# } else if (bodyRegion == "LateralTemporalLobe") {
#   inputGeneExpressionDataFilePath = "F://organizedAlzheimers//Setup//InputDataFiles//originalGeneExpressionDataSets//originalLateralTemporalLobeRegionGeneExpressionData.csv"
#   
#   #inputGeneExpressionDataFilePath = "F://organizedAlzheimers//GeneExpressionPreparation//mostUpdatedLateralTemporalLobeGroupedMeanNewest.csv"
#   phenotypesFilePath = "F://organizedAlzheimers//Setup//InputDataFiles//phenotypes//AlzheimersLateralTemporalLobePhenotypesUpdated.csv"
#   log2transformInputData = "TRUE" # (should we apply a log2(x+1) transform on data
#   scaleInputData = "FALSE" # should we apply a scale on data
#   
# }else if (bodyRegion == "DorsoLateralPrefrontalCortex") {
#  # inputGeneExpressionDataFilePath = "F://organizedAlzheimers//Setup//InputDataFiles//originalGeneExpressionDataSets//originalDLPFCRegionGeneExpressionData.csv"
#   #inputGeneExpressionDataFilePath = "F://organizedAlzheimers//DLPFC//updatedAveragedDLPFCDataSet.csv"
#   #inputGeneExpressionDataFilePath = "F://organizedAlzheimers//DLPFC//updatedAveragedDLPFCDataSetNewer.csv"
#   inputGeneExpressionDataFilePath = "F://organizedAlzheimers//DLPFC//finalAlzheimersDLPFCDataSetForSaniyaToUse.csv"
#   
#    # mostUpdatedDLPFCRegionGeneExpressionData
#   #inputGeneExpressionDataFilePath = "F://organizedAlzheimers//GeneExpressionPreparation//mostUpdatedLateralTemporalLobeGroupedMeanNewest.csv"
#   #phenotypesFilePath = "F://organizedAlzheimers//DLPFC//finalAlzheimersDLPFCDataSetForSaniyaToUse.csv"#pleaseUse//updatedAveragedDLPFCDataSet.csv"#AlzheimersLateralTemporalLobePhenotypesUpdated.csv"
#   
#   phenotypesFilePath = "F://organizedAlzheimers//Setup//InputDataFiles//phenotypes//AlzheimersDLPFCPhenotypesUpdatedForSaniya.csv"
#   log2transformInputData = "TRUE" # (should we apply a log2(x+1) transform on data
#   scaleInputData = "FALSE" # should we apply a scale on data
# } else if (bodyRegion == "myeloid") {
#    #inputGeneExpressionDataFilePath = "F://organizedAlzheimers//DLPFC//finalAlzheimersDLPFCDataSetForSaniyaToUse.csv"
#    inputGeneExpressionDataFilePath = "F://organizedAlzheimers//Setup//InputDataFiles//originalGeneExpressionDataSets//cellSampleMyeloidNEW.csv"
#    #F:\organizedAlzheimers\Setup\InputDataFiles\originalGeneExpressionDataSets
#   # mostUpdatedDLPFCRegionGeneExpressionData
#   #inputGeneExpressionDataFilePath = "F://organizedAlzheimers//GeneExpressionPreparation//mostUpdatedLateralTemporalLobeGroupedMeanNewest.csv"
#   #phenotypesFilePath = "F://organizedAlzheimers//DLPFC//finalAlzheimersDLPFCDataSetForSaniyaToUse.csv"#pleaseUse//updatedAveragedDLPFCDataSet.csv"#AlzheimersLateralTemporalLobePhenotypesUpdated.csv"
#   
#   phenotypesFilePath = "F://organizedAlzheimers//Setup//InputDataFiles//phenotypes//cellTypeTraits.csv" #"F://organizedAlzheimers//Setup//InputDataFiles//phenotypes//AlzheimersDLPFCPhenotypesUpdatedForSaniya.csv"
#   log2transformInputData = "TRUE" # (should we apply a log2(x+1) transform on data
#   scaleInputData = "TRUE" # should we apply a scale on data
# } else if (bodyRegion == "astrocyte") {
#   #inputGeneExpressionDataFilePath = "F://organizedAlzheimers//DLPFC//finalAlzheimersDLPFCDataSetForSaniyaToUse.csv"
#   inputGeneExpressionDataFilePath = "F://organizedAlzheimers//Setup//InputDataFiles//originalGeneExpressionDataSets//cellSampleAstrocyteNEW.csv"
#   #F:\organizedAlzheimers\Setup\InputDataFiles\originalGeneExpressionDataSets
#   # mostUpdatedDLPFCRegionGeneExpressionData
#   #inputGeneExpressionDataFilePath = "F://organizedAlzheimers//GeneExpressionPreparation//mostUpdatedLateralTemporalLobeGroupedMeanNewest.csv"
#   #phenotypesFilePath = "F://organizedAlzheimers//DLPFC//finalAlzheimersDLPFCDataSetForSaniyaToUse.csv"#pleaseUse//updatedAveragedDLPFCDataSet.csv"#AlzheimersLateralTemporalLobePhenotypesUpdated.csv"
#   
#   phenotypesFilePath = "F://organizedAlzheimers//Setup//InputDataFiles//phenotypes//cellTypeTraits.csv" #"F://organizedAlzheimers//Setup//InputDataFiles//phenotypes//AlzheimersDLPFCPhenotypesUpdatedForSaniya.csv"
#   log2transformInputData = "TRUE" # (should we apply a log2(x+1) transform on data
#   scaleInputData = "TRUE" # should we apply a scale on data
# } else if (bodyRegion == "neuron") {
#   #inputGeneExpressionDataFilePath = "F://organizedAlzheimers//DLPFC//finalAlzheimersDLPFCDataSetForSaniyaToUse.csv"
#   inputGeneExpressionDataFilePath = "F://organizedAlzheimers//Setup//InputDataFiles//originalGeneExpressionDataSets//cellSampleNeuronNEW.csv"
#   #F:\organizedAlzheimers\Setup\InputDataFiles\originalGeneExpressionDataSets
#   # mostUpdatedDLPFCRegionGeneExpressionData
#   #inputGeneExpressionDataFilePath = "F://organizedAlzheimers//GeneExpressionPreparation//mostUpdatedLateralTemporalLobeGroupedMeanNewest.csv"
#   #phenotypesFilePath = "F://organizedAlzheimers//DLPFC//finalAlzheimersDLPFCDataSetForSaniyaToUse.csv"#pleaseUse//updatedAveragedDLPFCDataSet.csv"#AlzheimersLateralTemporalLobePhenotypesUpdated.csv"
#   
#   phenotypesFilePath = "F://organizedAlzheimers//Setup//InputDataFiles//phenotypes//cellTypeTraits.csv" #"F://organizedAlzheimers//Setup//InputDataFiles//phenotypes//AlzheimersDLPFCPhenotypesUpdatedForSaniya.csv"
#   log2transformInputData = "TRUE" # (should we apply a log2(x+1) transform on data
#   scaleInputData = "TRUE" # should we apply a scale on data
# } else if (bodyRegion == "endothelial") {
#   #inputGeneExpressionDataFilePath = "F://organizedAlzheimers//DLPFC//finalAlzheimersDLPFCDataSetForSaniyaToUse.csv"
#   inputGeneExpressionDataFilePath = "F://organizedAlzheimers//Setup//InputDataFiles//originalGeneExpressionDataSets//cellSampleEndothelialNEW.csv"
#   #F:\organizedAlzheimers\Setup\InputDataFiles\originalGeneExpressionDataSets
#   # mostUpdatedDLPFCRegionGeneExpressionData
#   #inputGeneExpressionDataFilePath = "F://organizedAlzheimers//GeneExpressionPreparation//mostUpdatedLateralTemporalLobeGroupedMeanNewest.csv"
#   #phenotypesFilePath = "F://organizedAlzheimers//DLPFC//finalAlzheimersDLPFCDataSetForSaniyaToUse.csv"#pleaseUse//updatedAveragedDLPFCDataSet.csv"#AlzheimersLateralTemporalLobePhenotypesUpdated.csv"
#   
#   phenotypesFilePath = "F://organizedAlzheimers//Setup//InputDataFiles//phenotypes//cellTypeTraits.csv" #"F://organizedAlzheimers//Setup//InputDataFiles//phenotypes//AlzheimersDLPFCPhenotypesUpdatedForSaniya.csv"
#   log2transformInputData = "TRUE" # (should we apply a log2(x+1) transform on data
#   scaleInputData = "TRUE" # should we apply a scale on data
# } else if (bodyRegion == "neuroscience") {
#   #inputGeneExpressionDataFilePath = "F://organizedAlzheimers//DLPFC//finalAlzheimersDLPFCDataSetForSaniyaToUse.csv"
#   inputGeneExpressionDataFilePath = "F://organizedAlzheimers//Setup//InputDataFiles//originalGeneExpressionDataSets//neuroscienceInfo.csv"
#   #F:\organizedAlzheimers\Setup\InputDataFiles\originalGeneExpressionDataSets
#   # mostUpdatedDLPFCRegionGeneExpressionData
#   #inputGeneExpressionDataFilePath = "F://organizedAlzheimers//GeneExpressionPreparation//mostUpdatedLateralTemporalLobeGroupedMeanNewest.csv"
#   #phenotypesFilePath = "F://organizedAlzheimers//DLPFC//finalAlzheimersDLPFCDataSetForSaniyaToUse.csv"#pleaseUse//updatedAveragedDLPFCDataSet.csv"#AlzheimersLateralTemporalLobePhenotypesUpdated.csv"
#   
#   phenotypesFilePath = "F://organizedAlzheimers//Setup//InputDataFiles//phenotypes//neurosciencePhenos.csv" #"F://organizedAlzheimers//Setup//InputDataFiles//phenotypes//AlzheimersDLPFCPhenotypesUpdatedForSaniya.csv"
#   log2transformInputData = "FALSE" # (should we apply a log2(x+1) transform on data
#   scaleInputData = "FALSE" # should we apply a scale on data
# } else if (bodyRegion == "covid") {
#   #inputGeneExpressionDataFilePath = "F://organizedAlzheimers//DLPFC//finalAlzheimersDLPFCDataSetForSaniyaToUse.csv"
#   inputGeneExpressionDataFilePath = "F://organizedAlzheimers//Setup//InputDataFiles//originalGeneExpressionDataSets//originalCovid19Severity.csv"
#   #F:\organizedAlzheimers\Setup\InputDataFiles\originalGeneExpressionDataSets
#   # mostUpdatedDLPFCRegionGeneExpressionData
#   #inputGeneExpressionDataFilePath = "F://organizedAlzheimers//GeneExpressionPreparation//mostUpdatedLateralTemporalLobeGroupedMeanNewest.csv"
#   #phenotypesFilePath = "F://organizedAlzheimers//DLPFC//finalAlzheimersDLPFCDataSetForSaniyaToUse.csv"#pleaseUse//updatedAveragedDLPFCDataSet.csv"#AlzheimersLateralTemporalLobePhenotypesUpdated.csv"
#   
#   phenotypesFilePath = "F://organizedAlzheimers//Setup//InputDataFiles//phenotypes//covidSeverity.csv" #"F://organizedAlzheimers//Setup//InputDataFiles//phenotypes//AlzheimersDLPFCPhenotypesUpdatedForSaniya.csv"
#   log2transformInputData = "FALSE" # (should we apply a log2(x+1) transform on data
#   scaleInputData = "FALSE" # should we apply a scale on data
# }
# 
# 




tfsUsed = "LambertAndJaspar"




powerVal = "Unknown"
#dataScaling = "log2TransformedAndScaled" #  <-- will come from WGCNA
fullDiseaseName = paste0(disease, "Disease")




# WGCNAwithSpecificPower
# WGCNA with recommended power

#transcriptionalFilePath = 
# 
# newDate = str_replace_all(Sys.Date(), "-", "_")
# csvPart = paste0(tfsUsed, "TFs_", disease, "_", bodyRegion, "_", dataScaling, "_power", powerVal, "_", newDate, "Final.csv")
# rDataPart = paste0(tfsUsed, "TFs_", disease, "_", bodyRegion, "_", dataScaling, "_power", powerVal, "_", newDate, "Final.RData")
# oldRDataPart = paste0(filePathOfRTN, tfsUsed, "TFs_", disease, "_", bodyRegion, "_", dataScaling, "_power", powerVal, "_", newDate, "_Original.RData")
# RTN_rDataPart = paste0(tfsUsed, "TFs_", disease, "_", bodyRegion, "_", dataScaling, "_power", powerVal, "_FinalDataObjects.RData")
# 
# 
# 
# genie3CsvPart = paste0(filePathOfGenie3, tfsUsed, "TFs_", disease, "_", bodyRegion, "_", dataScaling, "_power", powerVal, "_", newDate, "Final.csv")
# genie3RDataPart = paste0(filePathOfGenie3, tfsUsed, "TFs_", disease, "_", bodyRegion, "_", dataScaling, "_power", powerVal, "_", newDate, "Final.RData")
# oldGenie3RDataPart = paste0(filePathOfGenie3, tfsUsed, "TFs_", disease, "_", bodyRegion, "_", dataScaling, "_power", powerVal, "_", newDate, "_Original.RData")
# genie3_rDataPart = paste0(filePathOfGenie3, tfsUsed, "TFs_", disease, "_", bodyRegion, "_", dataScaling, "_power", powerVal, "_FinalDataObjects.RData")
# 
# # printing out a text file or csv with information on location of each element we have,
# # along with what it means
# 
# # SNP threshold, etc.
# # key KEGG Pathways to visualize
# # WGCNA userListEnrichment

