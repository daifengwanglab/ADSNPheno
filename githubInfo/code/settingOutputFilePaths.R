# please note that the user will then get a csv of an explanation of the output
# files that will be generated

source("parametersNeeded.R")
if (bodyRegion == "Hippocampus"){
  tissueName = "Brain"
  disease = "Alzheimers"
  powerVal = "Unknown"
  dataScaling = "log2TransformedAndScaled" #  <-- will come from WGCNA
  tfsUsed = "LambertAndJaspar"
  
  
}


source("parametersNeeded.R")
newDate = str_replace_all(Sys.Date(), "-", "_")
csvPart = paste0(tfsUsed, "TFs_", disease, "_", bodyRegion, "_", dataScaling, "_power", powerVal, "_", newDate, "Final.csv")
rDataPart = paste0(tfsUsed, "TFs_", disease, "_", bodyRegion, "_", dataScaling, "_power", powerVal, "_", newDate, "Final.RData")
oldRDataPart = paste0(filePathOfRTN, tfsUsed, "TFs_", disease, "_", bodyRegion, "_", dataScaling, "_power", powerVal, "_", newDate, "_Original.RData")
RTN_rDataPart = paste0(tfsUsed, "TFs_", disease, "_", bodyRegion, "_", dataScaling, "_power", powerVal, "_FinalDataObjects.RData")

genie3CsvPart = paste0(filePathOfGenie3, tfsUsed, "TFs_", disease, "_", bodyRegion, "_", dataScaling, "_power", powerVal, "_", newDate, "Final.csv")
genie3RDataPart = paste0(filePathOfGenie3, tfsUsed, "TFs_", disease, "_", bodyRegion, "_", dataScaling, "_power", powerVal, "_", newDate, "Final.RData")
oldGenie3RDataPart = paste0(filePathOfGenie3, tfsUsed, "TFs_", disease, "_", bodyRegion, "_", dataScaling, "_power", powerVal, "_", newDate, "_Original.RData")
genie3_rDataPart = paste0(filePathOfGenie3, tfsUsed, "TFs_", disease, "_", bodyRegion, "_", dataScaling, "_power", powerVal, "_FinalDataObjects.RData")
