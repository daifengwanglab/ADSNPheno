geneIDMappingFunction <- function(x){
  entrezIDMappingVec[x]
}
returnCurrentDateToFilePath <- function(){
  library("stringr")
  currentDate = Sys.Date()
  outputAddOn = paste0("_", str_replace_all(currentDate, "-", "_"))
  return(outputAddOn)
}
# https://github.com/juanbot/CoExpNets/blob/master/R/main.R

# corDistance = function(a,b,signed=TRUE){
#   # if(signed)
#   #   return(0.5 * (1 + WGCNA::corFast(a,b)))
#   # return(abs(WGCNA::corFast(a,b)))
#   if(signed)
#     return(0.5 * (1 + cor(a,b)))
#   return(abs(cor(a,b)))
# }

corDistance = function(a,b,signed=TRUE,cor.type="pearson"){
  if(cor.type=="pearson"){
    if(signed)
      #return(0.5 + 0.5*WGCNA::corFast(x=a,y=b)) #(Note they are equivalent)
      return(0.5 * (1 + stats::cor(a,b,use="pairwise.complete.obs")))
    return(abs(stats::cor(a,b)))
  }else{
    if(signed)
      #return(0.5 + 0.5*WGCNA::corFast(a,b)) #(Note they are equivalent)
      return(0.5 * (1 + stats::cor(a,b,method=cor.type)))
    return(abs(stats::cor(a,b,method=cor.type,use="pairwise.complete.obs")))
  }
}

createCentroidMatrix <- function(eigengenes){
  my.matrix <- NULL
  for(eigengene in eigengenes){
    my.matrix <- cbind(my.matrix,eigengene)
  }
  return(my.matrix)
}

getNewCentroids <- function(expr.data,partition.in.colors,centroid.labels,mgg){
  
  # if (sum(partition.in.colors == "grey") == 0) {
  #   eg.vectors = moduleEigengenes(expr.data,partition.in.colors)$eigengenes
  # }
  if(sum(partition.in.colors == "grey") < mgg)
    eg.vectors = moduleEigengenes(expr.data,partition.in.colors, excludeGrey=TRUE)$eigengenes
  else
    eg.vectors = moduleEigengenes(expr.data,partition.in.colors, excludeGrey=F)$eigengenes
  
  names(eg.vectors) <- substring(names(eg.vectors),3)
  eg.vectors <- eg.vectors[,centroid.labels]
  return(eg.vectors)
}

getExchangedGenes <- function(old.partition,new.partition){
  stopifnot(length(old.partition) == length(new.partition))
  return(old.partition[old.partition != new.partition])
}


getBestModuleCor <- function(gene,centroids,signed=TRUE){
  
  return(which.max(corDistance(centroids,gene,signed)))
}


geneNameToIDMappingFunction <- function(geneSymbolsVec, infoDF){
  infoDF = buildingOutBestMappingInfoFunction(infoDF, geneSymbolsVec)
  
  entrezIDMappingVec = as.vector(infoDF[,2])
  names(entrezIDMappingVec) = infoDF[,1]
  entrezIDMappingVec
  
  print(paste0(":) welcome to the geneNameToIDMappingFunction :)"))
  print(":) Here, we map gene symbols to the id using the mapping function from entrezIDMappingVec.")
  geneIDsVec = sapply(geneSymbolsVec,geneIDMappingFunction)
  print(paste0(":) Please note our updated ", length(geneIDsVec) , " gene ids: "))
  print(geneIDsVec)
  geneDataFrame = data.frame(geneSymbolsVec)
  geneDataFrame$entrezID  = as.vector(geneIDsVec)
  colnames(geneDataFrame)[0] = "geneSymbolName"
  print(":) please note the first few entries of geneDataFrame:  geneSymbolName and entrezID respectively:")
  print(head(geneDataFrame))
  return(geneDataFrame)
  
  
}

buildingOutBestMappingInfoFunction <- function(infoDF, geneSymbolsVec){
  print(":) please note that this function builds out the updated infoDF by getting any missing Entrez terms for us")
  library("biomaRt")
  library("org.Hs.eg.db")
  
  newEntrezIDs = mapIds(org.Hs.eg.db, geneSymbolsVec, 'ENTREZID', 'SYMBOL')
  mappedDF = data.frame(as.vector(newEntrezIDs))
  mappedDF$geneName = as.vector(names(newEntrezIDs))
  mappedDF
  updatedDF = mappedDF[,c(2, 1)]
  colnames(updatedDF) = colnames(infoDF)[1:2]
  head(updatedDF)
  rowsToExclude = which(is.na(updatedDF[,2])) # na values
  print(paste0(":( Please note that EntrezIDs are missing for ", length(rowsToExclude),
               " genes, including: "))
  print(updatedDF[rowsToExclude, 1])
  updatedDF = updatedDF[-rowsToExclude,]
  updatedDF$status = rep("known", nrow(updatedDF))
  print(paste0(":) original dim(infoDF) is: ", dim(infoDF)))
  
  infoDF = rbind(infoDF, updatedDF)
  dim(infoDF)
  infoDF = unique(infoDF)
  print(paste0(":) updated dim(infoDF) is: ", dim(infoDF)))
  return(infoDF)
}


identifyingTFsInGeneExpressionDataSet <- function(allsamples, tfsDF, geneToEntrezIDPathForGeneExpressionDataCSV){
  geneExpressionSamplesVec = row.names(allsamples)
  print(":) Please note that this function returns the list of Transcription Factors ")
  print("that are found in the gene expression data set by mapping on Common EntrezID")
  print("This avoids problems with different gene symbol names  :)")
  #geneDataFrame = geneNameToIDMappingFunction(geneExpressionSamplesVec, infoDF)
  geneDataFrame = read.csv(geneToEntrezIDPathForGeneExpressionDataCSV, header = TRUE)
  tfIds = as.vector(tfsDF[,2])
  geneExpressionIds = as.vector(geneDataFrame[,2])
  tfsFound = which(tfIds %in% geneExpressionIds)
  #tfsFound = match(as.vector(geneDataFrame[,2]), as.vector(tfsDF[,2]))
  print(paste0(":) Please note that ", length(tfsFound), " TFs are found in the gene expression data set"))
  geneSymbolsOfTFs = unique(as.vector(geneDataFrame[tfsFound,1]))
  return(geneSymbolsOfTFs)
}


# def geneExpressionDataSetToEntrezIDMapping(inputGeneExpressionDataFilePath,
# geneToEntrezIDPathForGeneExpressionDataCSV, 
# geneToEntrezIDPathForGeneExpressionDataRData):
geneExpressionDataSetToEntrezIDMapping <- function(inputGeneExpressionDataFilePath,
                                                   geneToEntrezIDPathForGeneExpressionDataCSV, 
                                                   geneToEntrezIDPathForGeneExpressionDataRData){
  alzheimers_meanRMAData_raw = read.csv(inputGeneExpressionDataFilePath, #ROSMAP.normdata2 - Copy.csv", 
                                        header = TRUE, row.names = 1) #read.csv(GeneExpressionFileName, 
  dim(alzheimers_meanRMAData_raw)
  head(alzheimers_meanRMAData_raw)
  entrezIDsVec = c()   
  # if our 2nd column has the entrezIDs, please collect them but remove from the dataframe :)
  if (grepl(colnames(alzheimers_meanRMAData_raw)[1], "EntrezID")){
    entrezIDsVec = alzheimers_meanRMAData_raw[,1]
    alzheimers_meanRMAData_raw = alzheimers_meanRMAData_raw[,-1]
  }  else if (grepl(colnames(alzheimers_meanRMAData_raw)[1], "entrezID")){
    entrezIDsVec = alzheimers_meanRMAData_raw[,1]
    alzheimers_meanRMAData_raw = alzheimers_meanRMAData_raw[,-1]
  }
  dim(alzheimers_meanRMAData_raw) # 13076    31
  
  geneAndEntrezIDMappingDF = data.frame(entrezIDsVec)
  colnames(geneAndEntrezIDMappingDF) = "entrezID"
  row.names(geneAndEntrezIDMappingDF) = row.names(alzheimers_meanRMAData_raw)
  geneAndEntrezIDMappingDF = unique(geneAndEntrezIDMappingDF)
  write.csv(geneAndEntrezIDMappingDF, geneToEntrezIDPathForGeneExpressionDataCSV)
  save(geneAndEntrezIDMappingDF, alzheimers_meanRMAData_raw, file = geneToEntrezIDPathForGeneExpressionDataRData)
  print(geneAndEntrezIDMappingDF)
  return(geneAndEntrezIDMappingDF)
}
#    header = TRUE,
#   check.names = FALSE, row.names = 1)
# "one could also start with normalized counts (or RPKM/FPKM data) and log-transform them using log2(x+1)"


selectGenes <- function(counts, min.count=10, N=0.90){
  
  lib.size <- colSums(counts)
  MedianLibSize <- median(lib.size)
  CPM.Cutoff <- min.count / MedianLibSize*1e6
  CPM <- edgeR::cpm(counts,lib.size=lib.size)
  
  min.samples <- round(N * ncol(counts))
  
  f1 <- genefilter::kOverA(min.samples, CPM.Cutoff)
  flist <- genefilter::filterfun(f1)
  keep <- genefilter::genefilter(CPM, flist)
  
  ## the same as:
  #keep <- apply(CPM, 1, function(x, n = min.samples){
  #  t = sum(x >= CPM.Cutoff) >= n
  #  t
  #})
  
  return(keep)
}

creatingDatExprFromGeneExpressionDataset <- function(inputGeneExpressionDataFilePath,
                                                     geneToEntrezIDPathForGeneExpressionDataCSV, 
                                                     geneToEntrezIDPathForGeneExpressionDataRData,
                                                     geneExpressionOutputNameCSV_withDate,
                                                     geneExpressionOutputNameCSV,
                                                     originalGeneFilteredCSV,
                                                     originalGeneFilteredRData,
                                                     geneExpressionOutputNameRData,
                                                     geneExpressionOutputNameRData_withDate,
                                                     updatedOutputAddOnAndDatExpr_withDate,
                                                     updatedOutputAddOnAndDatExpr
){
  print(":) Please note this is function: creatingDatExprFromGeneExpressionDataset by Saniya")
  alzheimers_meanRMAData_raw = read.csv(inputGeneExpressionDataFilePath, #ROSMAP.normdata2 - Copy.csv", 
                                        header = TRUE, row.names = 1) 
  if (grepl(colnames(alzheimers_meanRMAData_raw)[1], "EntrezID")){
    alzheimers_meanRMAData_raw = alzheimers_meanRMAData_raw[,-1]
  } else if (grepl(colnames(alzheimers_meanRMAData_raw)[1], "entrezID")){
    alzheimers_meanRMAData_raw = alzheimers_meanRMAData_raw[,-1]
  }
  geneAndEntrezIDMappingDF = geneExpressionDataSetToEntrezIDMapping(inputGeneExpressionDataFilePath,
                                                                    geneToEntrezIDPathForGeneExpressionDataCSV, 
                                                                    geneToEntrezIDPathForGeneExpressionDataRData)
  currentDate = as.character(now()) # current date and time str_replace_all(Sys.Date(), "-", "_")
  
  alzheimersDFold = as.data.frame(t(alzheimers_meanRMAData_raw))
  # We first check for genes and samples with too many missing values:
  gsg = goodSamplesGenes(alzheimersDFold, verbose = 4)#3);
  gsg$allOK # this is TRUE, so all genes have passed the test :)
  library("edgeR")
  #keep.exprs <- filterByExpr(alzheimersDFold, min.count=10)
  #filt1 <- x[keep.exprs,]
  #keep.exprs <- selectGenes(alzheimersDFold, min.count=10, N=0.90)
  
  
  if (!gsg$allOK)
  {
    # Optionally, print the gene and sample names that were removed:
    if (sum(!gsg$goodGenes)>0)
      printFlush(paste("Removing genes:", paste(names(alzheimersDFold)[!gsg$goodGenes], collapse = ", ")));
    if (sum(!gsg$goodSamples)>0)
      printFlush(paste("Removing samples:", paste(rownames(alzheimersDFold)[!gsg$goodSamples], collapse = ", ")));
    # Remove the offending genes and samples from the data:
    alzheimersDFold = alzheimersDFold[gsg$goodSamples, gsg$goodGenes]
    write.csv(alzheimersDFold, originalGeneFilteredCSV)
    save(alzheimersDFold, currentDate, file = originalGeneFilteredRData)
  }
  
  
  alzheimersDFnew = alzheimers_meanRMAData_raw[gsg$goodGenes, gsg$goodSamples]
  
  
  
  if (log2transformInputData == "TRUE"){
    sdataAlzh2 = log2(alzheimersDFnew + 1)  # log2 transform of the data
  } else {
    sdataAlzh2 = alzheimersDFnew  # no log2 transform applied
    
  }
  if (scaleInputData == "TRUE"){
    ndataAlzh2 = scale(sdataAlzh2) # finding the z-score values for this data
  } else {
    ndataAlzh2 = sdataAlzh2 # no scale transform applied
    
  }
  
  #sdataAlzh2 = log2(alzheimers_meanRMAData_raw + 1)  # log2 transform of the data
  #ndataAlzh2 = scale(sdataAlzh2) # finding the z-score values for this data
  colnames(ndataAlzh2) = str_replace_all(colnames(ndataAlzh2), "X", "")
  
  write.csv(ndataAlzh2, geneExpressionOutputNameCSV_withDate)
  write.csv(ndataAlzh2, geneExpressionOutputNameCSV)
  
  
  
  #ndataAlzh2 = alzheimers_meanRMAData_raw
  alzheimersDF = as.data.frame(ndataAlzh2)
  head(alzheimersDF)
  dim(alzheimersDF)
  colnames(alzheimersDF) = str_replace_all(colnames(alzheimersDF), "X", "")

  colnames(alzheimersDFnew)
  
  datExpr = as.data.frame(t(alzheimersDF))
  nGenes = ncol(datExpr)
  nSamples = nrow(datExpr)
  print(paste0(":) Please note that there are: ", nGenes, " genes and ",
               nSamples, " total patient samples."))
  rm(alzheimersDF)
  
  print(paste0(":) Please note the # of genes in the final gene expression dataset is: ", nGenes))
  print(paste0(":) Please note the # of samples in the final gene expression dataset is: ", nSamples))
  
  save(ndataAlzh2, datExpr, geneAndEntrezIDMappingDF, nGenes, nSamples, currentDate, file = geneExpressionOutputNameRData)
  save(ndataAlzh2, datExpr, geneAndEntrezIDMappingDF, nGenes, nSamples, currentDate, file = geneExpressionOutputNameRData_withDate)
  write.csv(datExpr, updatedOutputAddOnAndDatExpr_withDate)
  write.csv(datExpr, updatedOutputAddOnAndDatExpr)
  print(dim(datExpr))
  print(head(datExpr))
  return(datExpr)
}



returningSoftThresholdPowerResults <- function(datExpr = datExpr, maxNumPowers = maxNumPowers,
                                               netType = netType, ourPower = ourPower, sftThresholdingOutputNameCSV = sftThresholdingOutputNameCSV) {
  # Choose a set of soft-thresholding powers
  powers = c(1:maxNumPowers)
  powers
  # [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21
  # [22] 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42
  # [43] 43 44 45 46 47 48 49 50
  # Call the network topology analysis function
  print(paste0("please note :) we will check soft-thresholding power values for powers from 1 to maxNumPowers = ", maxNumPowers))
  #library("WGCNA")
  
  # 
  # print(paste0(":) Please note that there are: ", nGenes, " genes and ",
  #              nSamples, " total patient samples."))
  # 
  geneNames = colnames(datExpr) # list of genes :)
  
  sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, networkType = netType)
  print("Please note the results from the Soft-Threshold Selection:")
  print(sft$fitIndices)
  sftDF = data.frame(sft$fitIndices)
  powerEstimate = sft$powerEstimate
  sftDF$recommendedPower = rep(powerEstimate, nrow(sftDF))
  
  sftThresholdingOutputNameCSV = paste0(wgcnaAndKMeansPathsForOutputs, "//softThresholdPowerDF_",netType, "_network_recommendedPower_",fullDiseaseName, "_", bodyRegion, outputAddOn, ".csv")
  sftThresholdingOutputNameRData = paste0(wgcnaAndKMeansPathsForOutputs, "//softThresholdPowerDF_",netType, "_network_recommendedPower_",fullDiseaseName, "_", bodyRegion, outputAddOn, ".RData")
  
  print(paste0("Please note that the Soft Threshold Power Data Frame has been written here: ", sftThresholdingOutputNameCSV))
  write.csv(sft$fitIndices, sftThresholdingOutputNameCSV)
  save(sftDF, powerEstimate, file = sftThresholdingOutputNameRData)
  ## methodology: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-02-networkConstr-man.pdf
  
  print(paste0("Please note that the Soft-Threshold Power Selected by WGCNA is: ", powerEstimate))
  if(is.na(powerEstimate)){
    print(paste0(":( Please note that this is NULL as an estimate for the power.  Please re-run analysis."))  
    if(log2transformInputData && scaleInputData){ # (should we apply a log2(x+1) transform on data
      print(paste0("Please try setting scaleInputData = FALSE instead of TRUE and try this function again with a log2(x+1) transform only (on the gene expression data) instead."))
    }
    break # return
  }
  
  
  return(sft)
}


selectingWGCNAPower <- function(datExpr = datExpr, maxNumPowers = maxNumPowers,
                                netType = netType, ourPower = ourPower, useRecommendedPower = useRecommendedPower,
                                log2transformInputData = log2transformInputData, scaleInputData = scaleInputData){ 
  # strategy if it is null for sft$powerEstimate
  
  sft = returningSoftThresholdPowerResults(datExpr, maxNumPowers,
                                           netType, ourPower, sftThresholdingOutputNameCSV) 


if (useRecommendedPower == "TRUE"){
  powerEstimate = sft$powerEstimate
  if(is.na(powerEstimate)){
    print(paste0(":( Please note that this is NULL as an estimate for the power.  Please re-run analysis."))  
    if(log2transformInputData && scaleInputData){ # (should we apply a log2(x+1) transform on data
      print(paste0("Please try setting scaleInputData = FALSE instead of TRUE and try this function again with a log2(x+1) transform only (on the gene expression data) instead."))
    }
    break # return
  }
} else {
  if (is.null(ourPower)){
    print("Please note ye did not specify a power to use.  Please note we are using WGCNA default instead.  Otherwise, please specify a non-null power")
  } else if (ourPower <= 0 || is.integer(ourPower) == FALSE) {
    print(paste0(":( Please note that ye specified ", ourPower, " as your selected power.  However, the power should be an integer that is 1 and above.  Please note we are using WGCNA default instead."))
    
  } else {
    powerEstimate = ourPower
    print(paste0("Please note we will instead use the power ye gave us: ", ourPower))
    
  }
}
print(paste0(":) ************* \n Please note Power Being Used: ", powerEstimate))
print(paste0(":) ************* \n Please note Network Type: ", netType))
return(powerEstimate)
}




# please note that this function by Saniya returns information on a numeric vector
# in terms of the key percentiles
percentileAssignmentsForAVector <- function(numericVectorToAnalyze, percentilesToUseVec, numRoundingDigits){
  #numericVectorToAnalyze = moduleInteractionsDF$TOMvalue
  percentilesVec = percentilesToUseVec
  #percentilesVec = percentilesForCoexpressionNetworkVec
  print("Please note that this function by Saniya returns information on a numeric vector in terms of the key percentiles.")
  numericVectorToAnalyze = as.numeric(numericVectorToAnalyze)
  percentileResults = quantile(numericVectorToAnalyze, probs = percentilesVec)
  
  # please note this percentileGroupingVec holds the groupings for the values in 
  # the numeric vector based on percentile groupings :)
  percentileGroupingVec = rep(NA, length(numericVectorToAnalyze))
  groupBoundariesVec =rep(NA, length(numericVectorToAnalyze))
  for (i in 1:(length(percentilesVec) - 1)){
    j = i + 1
    # getting the respective percentile values from the percentileResults
    minPercentileResultVal = as.numeric(percentileResults[i])
    maxPercentileResultVal = as.numeric(percentileResults[j])
    
    intervalResults = paste0(round(minPercentileResultVal, numRoundingDigits),
                             " to ", round(maxPercentileResultVal, numRoundingDigits))
    
    # getting the associated percentile #s
    minPercentileInInterval = percentilesVec[i] * 100
    maxPercentileInInterval = percentilesVec[j] * 100
    # please note this is the top grouping
    if (minPercentileInInterval == 0){ # ex. please note: bottom 25%
      groupName = paste0("Bottom ", maxPercentileInInterval, "%")
    } else if (maxPercentileInInterval == 100){ # ex. please note: top 5%
      differenceFrom100 = 100 - minPercentileInInterval
      groupName = paste0("Top ", differenceFrom100, "%")
      
    } else { # all other situations like 25% to 50% etc.
      groupName = paste0(minPercentileInInterval, "% to ", 
                         maxPercentileInInterval, "%")
      
    }
    if (j == length(percentilesVec)){ # if this is last interval, lets get all
      indexesInGroupVec = intersect(which(numericVectorToAnalyze >= minPercentileResultVal),
                                    which(numericVectorToAnalyze <= maxPercentileResultVal))
    } else {
      indexesInGroupVec = intersect(which(numericVectorToAnalyze >= minPercentileResultVal),
                                    which(numericVectorToAnalyze < maxPercentileResultVal))
    }
    percentileGroupingVec[indexesInGroupVec] = groupName
    groupBoundariesVec[indexesInGroupVec]= intervalResults
  }
  print(paste0(":) Please note the breakdown for these ", length(percentileGroupingVec), "entries: "))
  print(table(percentileGroupingVec))
  print(table(groupBoundariesVec))
  summaryStatisticsOnPercentilesDF = unique(cbind(percentileGroupingVec, groupBoundariesVec))
  percentileResultsToReturnList = list()
  percentileResultsToReturnList$summaryStatisticsOnPercentilesDF = summaryStatisticsOnPercentilesDF
  percentileResultsToReturnList$groupBoundariesVec = groupBoundariesVec
  percentileResultsToReturnList$percentileGroupingVec = percentileGroupingVec
  print(":) ******************************************** :)")
  print(paste0(":) Please note we return percentileResultsToReturnList, which is a list with 3 elements:"))
  print(paste0("\n\n:) percentileResultsToReturnList$summaryStatisticsOnPercentilesDF: has summary statistics for  percentile intervals, and results that define them"))
  print(paste0("\n :) percentileResultsToReturnList$groupBoundariesVec: has the group boundaries (percentile results) assignment for each entry in the numericVectorToAnalyze (ex. 0.004 to 0.025)"))
  print(paste0("\n :) percentileResultsToReturnList$percentileGroupingVec: has the group names (percentile results) assignment for each entry in the numericVectorToAnalyze (ex. Bottom 25%)"))
  print(str(percentileResultsToReturnList))
  return(percentileResultsToReturnList)
  
}


pleaseGetCoexpressionNetworkForAWGCNAModuleAfterKmeans <- function(module, 
                                                                   percentilesVec = percentilesForCoexpressionNetworkVec,
                                                                   datExprMat = datExpr,
                                                                   moduleColors = dynamicColorsKmeans,
                                                                   TOM = TOMPower, powerUsed = powerEstimate,
                                                                   networkType = netType,
                                                                   dataScale = dataScaling,
                                                                   bodyRegionName = bodyRegion,
                                                                   tissue = tissueName,
                                                                   diseaseName = fullDiseaseName,
                                                                   numRoundingDigits = numberOfRoundingDigits,
                                                                   minTOMThreshold = minTOMThresholdValue
){
  
  print(":) pleaseGetCoexpressionNetworkForAWGCNAModuleAfterKmeans")
  print(paste0(":) Please note that this function by Saniya tries to find the co-expression module ",
               " relationships between genes in a module to build a module co-expression network"))
  # please load TOM matrix
  # Select module
  #module = uniqueModules[1]
  # Select module probes
  probes = names(datExprMat)
  inModule = (moduleColors==module);
  modProbes = probes[inModule];
  # Select the corresponding Topological Overlap
  modTOM = TOM[inModule, inModule];
  modTOM[upper.tri(modTOM)] <- NA
  dimnames(modTOM) = list(modProbes, modProbes)
  
  moduleInteractionsDF = reshape2::melt(modTOM, varnames = c('row', 'col'), na.rm = TRUE) 
  
  print(dim(moduleInteractionsDF))
  print(head(moduleInteractionsDF))
  mainDiagIndices = which(moduleInteractionsDF[,1] == moduleInteractionsDF[,2])
  # removing main diagonal (since that is the same gene again)
  moduleInteractionsDF = moduleInteractionsDF[-mainDiagIndices,]
  print(paste0("Please note that this original dataframe for gene co-expression module ", 
               module, " has: ", nrow(moduleInteractionsDF), " total unique edges and ", 
               ncol(moduleInteractionsDF), " total columns"))
  str(moduleInteractionsDF)
  print(paste0(":) please note we are filtering out (and removing) the edges with a TOM value that is below: minTOMThreshold = ", minTOMThreshold))
  colnames(moduleInteractionsDF) = c("gene1", "gene2", "TOMvalue")
  
  # dataFrame[rows, columns]
  edgesToFilterRemove = which(moduleInteractionsDF$TOMValue < minTOMThreshold)
  numEdgesFilteredOut = length(edgesToFilterRemove)
  
  if (numEdgesFilteredOut == 0){
    print(paste0(":) Please note that all ", nrow(moduleInteractionsDF), 
                 "TOM edge values in the ", module, " were greater than minTOMThreshold = ",
                 minTOMThreshold, ".  --> NO Edges were removed."))
    print(paste0("Please set minTOMThreshold higher than ", min(as.numeric(moduleInteractionsDF[,3])), 
                 " to remove edges. Thank ye :)"))
  } else {
    print(paste0(":) Please note we are removing ", numEdgesFilteredOut, " that have edges with TOM values than ", minTOMThreshold))
    moduleInteractionsDF = moduleInteractionsDF[-edgesToFilterRemove,]
    print(paste0(":) Please note we will calculate metrics for the module and calculations based on these remaining edges"))
  }
  moduleInteractionsDF = unique(moduleInteractionsDF)
  
  numModuleGenes = length(which(inModule))
  
  moduleInteractionsDF$edgesConsidered = rep(paste0("only: all edges between any 2 genes in ", module, " module"),
                                             nrow(moduleInteractionsDF))
  
  rankOrderingVec = order(moduleInteractionsDF$TOMvalue, decreasing = TRUE)
  if (numEdgesFilteredOut == 0){
    rankOrder = paste0("edge rank: ", rankOrderingVec, " of ", nrow(moduleInteractionsDF), 
                       " edges within ", module, " gene module.")
  } else {
    rankOrder = paste0("edge rank: ", rankOrderingVec, " of ", nrow(moduleInteractionsDF), 
                       " remaining (unfiltered) edges within ", module, " gene module.")
  }
  # bottom 25%
  # 25% to 50%
  # 50% to 75%
  # 75% to 85%
  # 85 to 90%
  # 90 to 95%
  # 95 and above 
  
  
  percentileResultsToReturnList = percentileAssignmentsForAVector(moduleInteractionsDF$TOMvalue, percentilesVec, numRoundingDigits)
  
  
  moduleInteractionsDF$moduleForGenes =  rep(module, nrow(moduleInteractionsDF))
  moduleInteractionsDF$numOfGenesInThatModule =  rep(numModuleGenes, nrow(moduleInteractionsDF))
  moduleInteractionsDF$totalNumOfGenesInGeneExpressionDataset =  rep(ncol(datExpr), nrow(moduleInteractionsDF))
  
  if (numEdgesFilteredOut == 0){
    moduleInteractionsDF$analysis =  rep(paste0("WGCNA with KMeans (Results Calculated on ALL Edges for ", module, " Module)"), nrow(moduleInteractionsDF))
    
  } else {
    moduleInteractionsDF$analysis =  rep(paste0("WGCNA with KMeans (Results Calculated on Filtered List of Edges for ", module, " Module)"), nrow(moduleInteractionsDF))
    
  }
  #dataScale = dataScaling,
  #bodyRegionName = bodyRegion,
  moduleInteractionsDF$disease =  rep(diseaseName, nrow(moduleInteractionsDF))
  moduleInteractionsDF$tissue =  rep(tissueName, nrow(moduleInteractionsDF))
  moduleInteractionsDF$bodyRegion = rep(bodyRegionName, nrow(moduleInteractionsDF))
  moduleInteractionsDF$dataScaling = rep(dataScale, nrow(moduleInteractionsDF))
  moduleInteractionsDF$powerUsed = rep(powerUsed, nrow(moduleInteractionsDF))
  moduleInteractionsDF$networkType = rep(networkType, nrow(moduleInteractionsDF))
  
  moduleInteractionsDF$modulePercentileGroupForTOM =  percentileResultsToReturnList$percentileGroupingVec
  moduleInteractionsDF$moduleTomIntervalForPercentileGroup =  percentileResultsToReturnList$groupBoundariesVec
  moduleInteractionsDF$rankOrder = rankOrderingVec
  moduleInteractionsDF$rankOrderExplanation = rankOrder
  moduleInteractionsDF$minTOMThreshold = rep(minTOMThreshold, nrow(moduleInteractionsDF))
  
  if (numEdgesFilteredOut == 0){
    moduleInteractionsDF$numEdgesFilteredOut = rep("no edges filtered out", nrow(moduleInteractionsDF))
    
  } else {
    moduleInteractionsDF$numEdgesFilteredOut = rep(edgesFilteredOut, nrow(moduleInteractionsDF))
    
  }
  moduleInteractionsDF$totalNumOfEdgesAfterFiltering = rep(nrow(moduleInteractionsDF), nrow(moduleInteractionsDF))
  
  
  head(moduleInteractionsDF) 
  dim(moduleInteractionsDF)
  
  wgcnaAndKmeansModuleCoexpressionNetworkForModuleALLInteractionsCSV = paste0(wgcnaAndKMeansPathsForOutputs, "//CoexpressionNetwork_WGCNAAndKmeans//", bodyRegion, "_wgcna_with_kmeans_Power", powerEstimate, "_CoexpressionNetworkForModule_", module,outputAddOn,".csv")
  wgcnaAndKmeansModuleCoexpressionNetworkForModuleALLInteractionsRData = paste0(wgcnaAndKMeansPathsForOutputs, "//CoexpressionNetwork_WGCNAAndKmeans//", bodyRegion, "_wgcna_with_kmeans_Power", powerEstimate, "_CoexpressionNetworkForModule_", module,outputAddOn,".RData")
  
  
  wgcnaAndKmeansModuleCoexpressionNetworkForModuleTOMPercentileSummaryStatisticsCSV = paste0(wgcnaAndKMeansPathsForOutputs, "//CoexpressionNetwork_WGCNAAndKmeans//", bodyRegion, "_wgcna_with_kmeans_Power", powerEstimate, "_TOMPercentileSummaryStatisticsForModule_", module,outputAddOn,".csv")
  currentDate = as.character(now())
  write.csv(moduleInteractionsDF, wgcnaAndKmeansModuleCoexpressionNetworkForModuleALLInteractionsCSV)
  save(moduleInteractionsDF, currentDate, file = wgcnaAndKmeansModuleCoexpressionNetworkForModuleALLInteractionsRData)
  
  # generating the summary statistics based on the TOM edge values
  summaryStatisticsOnPercentilesDF = data.frame(percentileResultsToReturnList$summaryStatisticsOnPercentilesDF)
  colnames(summaryStatisticsOnPercentilesDF) = c("modulePercentileGroupForTOM", "moduleTomIntervalForPercentileGroup")
  summaryStatisticsOnPercentilesDF$module =  rep(module, nrow(summaryStatisticsOnPercentilesDF))
  summaryStatisticsOnPercentilesDF$numOfGenesInThatModule =  rep(numModuleGenes, nrow(summaryStatisticsOnPercentilesDF))
  summaryStatisticsOnPercentilesDF$minTOMThreshold = rep(minTOMThreshold, nrow(summaryStatisticsOnPercentilesDF))
  
  if (numEdgesFilteredOut == 0){
    summaryStatisticsOnPercentilesDF$numEdgesFilteredOut = rep("no edges filtered out", nrow(summaryStatisticsOnPercentilesDF))
    
  } else {
    summaryStatisticsOnPercentilesDF$numEdgesFilteredOut = rep(edgesFilteredOut, nrow(summaryStatisticsOnPercentilesDF))
    
  }
  summaryStatisticsOnPercentilesDF$totalNumOfEdgesAfterFiltering = rep(nrow(moduleInteractionsDF), nrow(summaryStatisticsOnPercentilesDF))
  summaryStatisticsOnPercentilesDF$edgesConsidered = rep(paste0("only: all edges between any 2 genes in ", module, " module"),
                                                         nrow(summaryStatisticsOnPercentilesDF))
  summaryStatisticsOnPercentilesDF$dateTheFileWasCreated =  rep(currentDate, nrow(summaryStatisticsOnPercentilesDF))
  
  write.csv(summaryStatisticsOnPercentilesDF, wgcnaAndKmeansModuleCoexpressionNetworkForModuleTOMPercentileSummaryStatisticsCSV)
  
  moduleCoExpressionResultsList = list()
  
  moduleCoExpressionResultsList$module = module 
  moduleCoExpressionResultsList$moduleInteractionsDF = moduleInteractionsDF
  moduleCoExpressionResultsList$summaryStatisticsOnPercentilesDF = summaryStatisticsOnPercentilesDF
  
  print(":) ******************************************** :)")
  print(":) Please note we return moduleCoExpressionResultsList, which is a list with 3 elements:")
  print("\n\n :) moduleCoExpressionResultsList$module: has the summary statistics in a dataframe for the percentile intervals, and results that define them")
  print("\n :) moduleCoExpressionResultsList$moduleInteractionsDF: has the interactions for the module genes based on the Topological Overlap Matrix (TOM)")
  print(paste0("\n :) moduleCoExpressionResultsList$summaryStatisticsOnPercentilesDF:  has the updated summary statistics for the ", module, 
               " gene coexpression module Topological Overlap Matrix (TOM) percentile intervals, and results that define them"))
  
  print(paste0("\n\n Please note we are done with gene module: ", module))
  print(str(moduleCoExpressionResultsList))
  return(moduleCoExpressionResultsList)
  
}

# calculate significant positive phenotypes for each gene please

pleaseGetTheGeneTraitSignificanceValuesForTheTraits <- function(datTraits,
                                                                datExpr,
                                                                wgcnaWithKmeansGeneTraitSignificancePath,
                                                                wgcnaWithKmeansGeneTraitSigPValuePath,
                                                                bodyRegion, tissueName, fullDiseaseName, dataScaling, geneAndEntrezIDMappingDF) {
  gene_names = colnames(datExpr)
  
  GSPvalueDF = data.frame()
  geneTraitSignificanceDF = data.frame()
  for (featureIndex in 1:ncol(datTraits)){
    print(paste0("featureIndex: ", featureIndex))
    weight = as.data.frame(datTraits[,featureIndex])
    weightName = colnames(datTraits)[featureIndex]
    geneTraitSignificance = as.data.frame(cor(datExpr, weight, use = "p"));
    colnames(geneTraitSignificance) = paste("GeneSignificance_", weightName, sep="");
    nSamples = nrow(datExpr)
    GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
    colnames(GSPvalue) = paste("GeneSignificancePValue_", weightName, sep="");
    
    #names(geneTraitSignificance) = paste("GeneSignificance", weightName, sep="");
    #names(GSPvalue) = paste("GeneSignificancePValue", weightName, sep="");
    if (nrow(GSPvalueDF) == 0){
      GSPvalueDF = GSPvalue
      
    } else {
      GSPvalueDF = cbind(GSPvalueDF, GSPvalue)
      
    }
    if (nrow(geneTraitSignificanceDF) == 0){
      geneTraitSignificanceDF = geneTraitSignificance
      
    } else {
      geneTraitSignificanceDF = cbind(geneTraitSignificanceDF, geneTraitSignificance)
    }
  }
  geneTraitSignificanceDF$explanation = paste0("Correlation between gene ", row.names(geneTraitSignificanceDF), 
                                               " and each ", disease, " trait observed in the ", bodyRegion, " of the ",
                                               tissueName, " for the ", dataScaling, " gene expression data.")
  
  geneTraitSignificanceDF$region = rep(bodyRegion, nrow(geneTraitSignificanceDF))
  geneTraitSignificanceDF$disease = rep(fullDiseaseName, nrow(geneTraitSignificanceDF))
  geneTraitSignificanceDF$dataScaling = rep(dataScaling, nrow(geneTraitSignificanceDF))
  
  
  GSPvalueDF$explanation = paste0("Associated p-value for the correlation between gene ", row.names(geneTraitSignificanceDF), 
                                  " and each ", disease, " trait observed in the ", bodyRegion, " of the ",
                                  tissueName, " for the ", dataScaling, " gene expression data.")
  
  GSPvalueDF$region = rep(bodyRegion, nrow(geneTraitSignificanceDF))
  GSPvalueDF$disease = rep(fullDiseaseName, nrow(geneTraitSignificanceDF))
  GSPvalueDF$dataScaling = rep(dataScaling, nrow(geneTraitSignificanceDF))
  
  geneTraitSignificanceDF = merge(geneAndEntrezIDMappingDF, geneTraitSignificanceDF, by = 0)
  rownames(geneTraitSignificanceDF)= geneTraitSignificanceDF[,1] #gene_names
  geneTraitSignificanceDF = geneTraitSignificanceDF[,-1]
  colnames(geneTraitSignificanceDF)[1] = "entrezID"
  GSPvalueDF = merge(geneAndEntrezIDMappingDF, GSPvalueDF, by = 0)
  colnames(GSPvalueDF)[1] = "entrezID"
  rownames(GSPvalueDF)= GSPvalueDF[,1] #gene_names
  GSPvalueDF = GSPvalueDF[,-1]
  head(GSPvalueDF)
  #colnames(clusterInfo) = c("module", "entrezID")
  
  
  
  write.csv(geneTraitSignificanceDF, wgcnaWithKmeansGeneTraitSignificancePath)
  write.csv(GSPvalueDF, wgcnaWithKmeansGeneTraitSigPValuePath)
  
  geneTraitSignificanceResultsList = list()
  geneTraitSignificanceResultsList$geneTraitSignificanceDF = geneTraitSignificanceDF
  geneTraitSignificanceResultsList$GSPvalueDF = GSPvalueDF
  print(paste0(":) Please note we return the results for gene trait significance correlations and associated correlation p-values here: ",
               "geneTraitSignificanceResultsList"))
  
  print(paste0("Please note :) geneTraitSignificanceResultsList$geneTraitSignificanceDF has the geneTraitSignificance correlations"))
  print(paste0("Please note :) geneTraitSignificanceResultsList$GSPvalueDF has the p-values for the geneTraitSignificance correlations"))
  
  print(str(geneTraitSignificanceResultsList))
  return(geneTraitSignificanceResultsList)
}
