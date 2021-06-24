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

   colnames(ndataAlzh2) = str_replace_all(colnames(ndataAlzh2), "X", "")
if (includeTimeStamp == "TRUE"){
  write.csv(ndataAlzh2, geneExpressionOutputNameCSV_withDate)
  }
  write.csv(ndataAlzh2, geneExpressionOutputNameCSV)

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
  geneNames = colnames(datExpr) # list of genes :)

  sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, networkType = netType)
  print("Please note the results from the Soft-Threshold Selection:")
  print(sft$fitIndices)
  sftDF = data.frame(sft$fitIndices)
  powerEstimate = sft$powerEstimate
  sftDF$recommendedPower = rep(powerEstimate, nrow(sftDF))

 
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
    break 
  }


  return(sft)
}


selectingWGCNAPower <- function(datExpr = datExpr, maxNumPowers = maxNumPowers,
                                netType = netType, ourPower = ourPower, useRecommendedPower = useRecommendedPower,
                                log2transformInputData = log2transformInputData, scaleInputData = scaleInputData){
  # strategy if it is null for sft$powerEstimate




if (useRecommendedPower == "TRUE"){
    print(paste0(":( Please note that since useRecommendedPower == TRUE, we are getting the best soft-thresholding beta power-estimate for ye from WGCNA :)"))#  is NULL as an estimate for the power.  Please re-run analysis."))

	sft = returningSoftThresholdPowerResults(datExpr, maxNumPowers,
                                           netType, ourPower, sftThresholdingOutputNameCSV)
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
  } else if (ourPower < 1 || is.integer(as.integer(ourPower))== FALSE) {
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
                                                                   minTOMThreshold = minTOMThresholdValue,
                                                                   wgcnaAndKMeansPathsForOutputs = pathForWGCNA
){

  dir.create(paste0(wgcnaAndKMeansPathsForOutputs, "//Module_CoexpressionNetwork//"))

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

  wgcnaAndKmeansModuleCoexpressionNetworkForModuleALLInteractionsCSV = paste0(wgcnaAndKMeansPathsForOutputs, "//Module_CoexpressionNetwork//", bodyRegion, "_wgcnaWithKmeans_Pow", powerEstimate, "_CoexpressionNetworkForModule_", module,outputAddOn,".csv")
  wgcnaAndKmeansModuleCoexpressionNetworkForModuleALLInteractionsRData = paste0(wgcnaAndKMeansPathsForOutputs, "//Module_CoexpressionNetwork//", bodyRegion, "_wgcnaWithKmeans_Pow", powerEstimate, "_CoexpressionNetworkForModule_", module,outputAddOn,".RData")

  wgcnaAndKmeansModuleCoexpressionNetworkForModuleTOMPercentileSummaryStatisticsCSV = paste0(wgcnaAndKMeansPathsForOutputs, "//Module_CoexpressionNetwork//", bodyRegion, "_wgcnaWithKmeans_Pow", powerEstimate, "_TOMPercentModSumStats_", module,outputAddOn,".RData")

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



##############################################################
### ENRICHMENTS:

# MEDICAL SUBJECT HEADINGS (MESH)

runAllMeSHEnrichments <- function(categoryVec,# = c("A", "B", "C", "D", "E", "F", "G", "H"),
                                  powerVal,# = 10,
                                  region,# = "Lateral Temporal Lobe",
                                  geneIdCol = entrezIDCol,
                                  geneClusterDFToUse,
                                  parentName = "D:\\organizedAlzheimers\\enrichments\\outputs\\MeSH Enrichment Analysis\\",
								  colNum = modColNum, kMeans = performAdditionalKMeansStep,
								  meshDatabaseNames = c("gendoo", "gene2pubmed", "RBBH")) {
  listOfMeshEnrichments <- list()

  for (i in 1:length(categoryVec)){
    categoryName = categoryVec[i]
    print("\n*****************************************")
    print(paste0(" i = ", i, " of ", length(categoryVec),
                 " :) Please note that the MeSH category is ", categoryName))
    meshDF = meshesFunctionsEnrichment(categoryName = categoryName,
                                       powerVal = powerVal,
                                       geneClusterDF = geneClusterDFToUse,
									   colNum = colNum,
									   kMeans = performAdditionalKMeansStep,
									   databaseNamesVec = meshDatabaseNames)
    listOfMeshEnrichments[categoryName] <- meshDF

  }
  print(str(listOfMeshEnrichments))
  print(":) Please note that we have returned a list of the different MeSH Enrichments :)")
  return(listOfMeshEnrichments)
}

#####################################

meshesFunctionsEnrichment <- function(categoryName, # = "A",
                                      powerVal, # = 10,
                                      region = regionName,# = "Lateral Temporal Lobe",
                                      geneIdCol = entrezIDCol,
                                      geneClusterDF,
                                      parentName = "D:\\organizedAlzheimers\\enrichments\\outputs\\MeSH Enrichment Analysis\\",
									  colNum = modColNum, kMeans = performAdditionalKMeansStep,   databaseNamesVec = c("gendoo", "gene2pubmed", "RBBH")){
  #categoryName = "H"
  if (kMeans){
  outputName = paste0(parentName, "mesh_", categoryName, "_", region,
                      "_WGCNA_with_KmeansAfter_pow_", powerVal, ".csv")
					  }
	else {
	  outputName = paste0(parentName, "mesh_", categoryName, "_", region,
                      "_WGCNA_Only_pow_", powerVal, ".csv")
	}
  print("******************************************************* \n")
  if (categoryName == "A"){
    fullCategoryName = "anatomyMedSubHeadCategA"
    longName = "A: Anatomy"
    print(paste0(":) Please note that the MeSH (Medical Subject Heading category selected is: ", longName, ", which will be encoded as: ",
                 fullCategoryName, " :)"))
  } else if (categoryName == "B"){
    fullCategoryName = "organismsMedSubHeadCategB"
    longName = "B: Organisms"
    print(paste0(":) Please note that the MeSH (Medical Subject Heading category selected is: ", longName, ", which will be encoded as: ",
                 fullCategoryName, " :)"))
  } else if (categoryName == "C"){
    fullCategoryName = "diseasesMedSubHeadCategC"
    longName = "C: Diseases"
    print(paste0(":) Please note that the MeSH (Medical Subject Heading category selected is: ", longName, ", which will be encoded as: ",
                 fullCategoryName, " :)"))
  } else if (categoryName == "D"){
    fullCategoryName = "chemicalsAndDrugsMedSubHeadCategD"
    longName = "D: Chemicals and Drugs"
    print(paste0(":) Please note that the MeSH (Medical Subject Heading category selected is: ", longName, ", which will be encoded as: ",
                 fullCategoryName, " :)"))
  } else if (categoryName == "E"){
    fullCategoryName = "diagnosticMedSubHeadCategE"
    longName = "E: Analytical, Diagnostic, and Therapeutic Techniques and Equipment"
    print(paste0(":) Please note that the MeSH (Medical Subject Heading category selected is: ", longName, ", which will be encoded as: ",
                 fullCategoryName, " :)"))
  } else if (categoryName == "F"){
    fullCategoryName = "psychiatryAndPsychologyMedSubHeadCategF"
    longName = "F: Psychiatry and Psychology"
    print(paste0(":) Please note that the MeSH (Medical Subject Heading category selected is: ", longName, ", which will be encoded as: ",
                 fullCategoryName, " :)"))
  } else if (categoryName == "G"){
    fullCategoryName = "phenomenaStateMedSubHeadCategG"
    longName = "G: Biological Sciences"
    print(paste0(":) Please note that the MeSH (Medical Subject Heading category selected is: ", longName, ", which will be encoded as: ",
                 fullCategoryName, " :)"))
  } else if (categoryName == "H"){
    fullCategoryName = "physicalScienceMedSubHeadCategH"
    longName = "H: Physical Sciences"
    print(paste0(":) Please note that the MeSH (Medical Subject Heading category selected is: ", longName, ", which will be encoded as: ",
                 fullCategoryName, " :)"))
  }

  outputName = paste0(parentName, "mesh_", categoryName, "_", fullCategoryName, "_", region,
                      "_WGCNA_with_KmeansAfter_power_", powerVal, "_ALL.csv")

  print(paste0(":) Please note this is MeSH (Medical Subject Headings) Enrichment Analysis, for category:", categoryName))

  uniqueMods = as.vector(unique(geneClusterDF[,colNum]))
  firstMod = as.vector(as.character(geneClusterDF[which(geneClusterDF$module == uniqueMods[1]),entrezIDCol]))

  #databaseNamesVec = c("gendoo", "gene2pubmed", "RBBH")


  diseaseStateMeSH_C_Diseases = data.frame()
  namesPart = ""
  for (j in 1:length(databaseNamesVec)){
    databaseName = databaseNamesVec[j]
    if (j == 1){ # for database name
      namesPart = databaseName
    } else if (j == 2) {
      namesPart = "gendoo_gene2pubmed"
    } else {
      namesPart = "gendoo_gene2pubmed_RBBH"
    }

    outputNameMini = paste0(parentName, "mesh_", categoryName, "_", fullCategoryName, "_", region,
                        "_WGCNA_with_KmeansAfter_power_", powerVal, "_", namesPart, "Only.csv")

    print(outputNameMini)
    for (i in 1:length(uniqueMods)){
      tryCatch(
        {
          mod = uniqueMods[i]
          print(paste0(":) i = ", i, " of ", length(uniqueMods),
                       " ; please note MeSH database used is ", databaseName,
                       " and power ", powerVal, " and module : " , mod))
          nextMod = as.vector(as.character(geneClusterDF[which(geneClusterDF[,colNum] == mod),entrezIDCol]))
          nextMod_MeSH_C_Diseases = enrichMeSH(nextMod, MeSHDb = "MeSH.Hsa.eg.db", database=databaseName, category = categoryName)

          minusLog10AdjustedPval = -1*log10(nextMod_MeSH_C_Diseases$p.adjust)
          database = rep(databaseName, nrow(nextMod_MeSH_C_Diseases))
          module = rep(mod, nrow(nextMod_MeSH_C_Diseases))
          power = rep(powerVal, nrow(nextMod_MeSH_C_Diseases))
          regionName = rep(region, nrow(nextMod_MeSH_C_Diseases))
          meshHeading = rep(fullCategoryName, nrow(nextMod_MeSH_C_Diseases))
          explanation = rep(longName, nrow(nextMod_MeSH_C_Diseases))
          nextMod_MeSH_C_Diseases_DF = data.frame(cbind(data.frame(nextMod_MeSH_C_Diseases),
                                                        minusLog10AdjustedPval,
                                                        meshHeading, explanation,
                                                        database,
                                                        module, power, regionName))

          if (nrow(diseaseStateMeSH_C_Diseases) == 0){  # new dataframe
            diseaseStateMeSH_C_Diseases = nextMod_MeSH_C_Diseases_DF
          } else {
            diseaseStateMeSH_C_Diseases = rbind(diseaseStateMeSH_C_Diseases, nextMod_MeSH_C_Diseases_DF)
          }
        }, error= function(cond) {
          message(cond)
          print(":( skipped")

        })

    }
    write.csv(diseaseStateMeSH_C_Diseases, outputNameMini)
  }
  print(paste0(":)! YAY!  Please note we have finished the enrichments here and wrote the output here: ",
               outputName))
  print(head(diseaseStateMeSH_C_Diseases))
  write.csv(diseaseStateMeSH_C_Diseases, outputName)
  return(diseaseStateMeSH_C_Diseases)
}


# Molecular Signatures DB:

molSigDBFunctionsEnrichment <- function(categoryName = "C1", powerVal = 10,
                                        geneIdCol = 2, colNum,
                                        region = "Lateral Temporal Lobe",
                                        geneClusterDF,
                                        parentName = "D:\\organizedAlzheimers\\enrichments\\outputs\\MolecularSignaturesDatabase\\",
										kMeans = performAdditionalKMeansStep){


  m_t2g_C5 <- msigdbr(species = "Homo sapiens", category = categoryName) %>%
    dplyr::select(gs_name, entrez_gene)
  head(m_t2g_C5)
  data.frame(m_t2g_C5)
  if (kMeans){
  outputName = paste0(parentName, "molSigDB_", categoryName, "_", region,
                      "_WGCNA_then_Kmeans_power_", powerVal, ".csv")
					  }
	else {
	  outputName = paste0(parentName, "molSigDB_", categoryName, "_", region,
                      "_WGCNA_Only_power_", powerVal, ".csv")
	}
  if (categoryName == "C1"){
    molSigDBHeading = "C1 Positional Gene Sets"

  } else if (categoryName == "C2"){
    molSigDBHeading = "C2 Curated Gene Sets"
  } else if (categoryName == "C3"){
    molSigDBHeading = "C3 Motif Gene Sets"
  } else if (categoryName == "C4"){
    molSigDBHeading = "C4 Computational Gene Sets"
  } else if (categoryName == "C5"){
    molSigDBHeading = "C5 GO Gene Sets"
  } else if (categoryName == "C6"){
    molSigDBHeading = "C6 Oncogenic Signatures"
  } else if (categoryName == "C7"){
    molSigDBHeading = "C7 Immunologic Signatures"
  } else if (categoryName == "H"){
    molSigDBHeading = "H Hallmark Gene Sets"
  }


  uniqueMods = as.vector(unique(geneClusterDF[,colNum]))


  C5_df = data.frame() #emDF

  for (i in 1:length(uniqueMods)){
    tryCatch(
      {
        mod = uniqueMods[i]
        print(paste0(i, ":) please note module ", powerVal, ": ", mod))
        nextMod = as.vector(as.character(geneClusterDF[which(geneClusterDF[,colNum] == mod),geneIdCol]))

        em <- enricher(nextMod, TERM2GENE=m_t2g_C5)
        emDF = data.frame(em)
        molecularSignaturesDatabase = rep(molSigDBHeading, nrow(emDF))
        module = rep(mod, nrow(emDF))
        power = rep(powerVal, nrow(emDF))
        emDF$molecularSignaturesDatabase = molecularSignaturesDatabase
        emDF$module = module
        emDF$power = power
        emDF$region = region

        if (nrow(C5_df) == 0){
          C5_df = emDF
        } else {
          C5_df = rbind(C5_df, emDF)
        }


      }, error= function(cond) {
        message(cond)
        print(":( skipped")
        #next

      })

  }
  print(paste0(":) Please note that the output path is: ", outputName))

  write.csv(C5_df, outputName)
  print(C5_df)
  return(C5_df)
}



##############################################
### STEP 5: CHROMATIN INTERACTION WORK:


# please note that this function processes the set of enhancers data for scgrnom pipeline
processEnhancerFile <- function(enhancerChromatinFilePath){
  enhancer = read.csv(enhancerChromatinFilePath, header = TRUE)#[,2:4]
  head(enhancer)
  # please find respective columns from the Enhancer Chromatin data
  # please note that we need the columns to be in this form:  chr, start, end
  chrCol = grep("chr", enhancer[1,])[1] # please take 1st entry with chr entry in it
  enhancerColNamesVec = c(chrCol, chrCol + 1, chrCol + 2) # please note the input should be chr, start, end
  enhancerColNamesVec
  enhancer = data.frame(enhancer)[,enhancerColNamesVec] # please select those respective columns
  head(enhancer)

  colnames(enhancer) = c("chr", "start", "end")
  enhancer = unique(enhancer)
  head(enhancer)
  return(enhancer)
}

# please note that this function processes the set of HIC Interaction (Promoter and Enhancer) data for scgrnom pipeline
processHICInteractionFile <- function(chromatinInteractionDataFilePath){
  hic_interaction = read.csv(chromatinInteractionDataFilePath, header = TRUE)
  dim(hic_interaction) # 2,451,345      7 -->  2,805,134       7
  hicColNamesVec = colnames(hic_interaction)
  hicColNamesVec
  # please find respective columns from the HI-C Chromatin interaction data
  # please note that we need the columns to be in this form:  chr1, start1, end1, chr2, start2, end2
  colsNeededHIC = c(grep("chr1", hicColNamesVec),
                    grep("start1", hicColNamesVec), grep("end1", hicColNamesVec),
                    grep("chr2", hicColNamesVec),
                    grep("start2", hicColNamesVec), grep("end2", hicColNamesVec))
  colsNeededHIC
  hic_interaction = data.frame(hic_interaction)[,colsNeededHIC]

  #hic_interaction = data.frame(hic_interaction)[,2:7]

  hic_interaction = unique(hic_interaction)
  head(hic_interaction)
  return(hic_interaction)
}


############ PLEASE NOTE THAT THIS CODE BELOW HAS BEEN TAKEN FROM SCGRNOM:


scGRN_interactionOLD = function(hic_interaction, enhancers, ref_promoters = 'all',up_stream = 2500,
                             down_stream = 2500, link_type = 'within' ,target_genes='all',
                             gene_id_option = 'hgnc_symbol',
                             mart = biomaRt::useMart(biomart="ENSEMBL_MART_ENSEMBL",
                                                     dataset="hsapiens_gene_ensembl",
                                                     host="uswest.ensembl.org")){


  ehs <- GenomicRanges::GRanges(seqnames = enhancers$chr, ranges = IRanges::IRanges(start = enhancers$start,
                                                                                    end = enhancers$end))
  names(ehs) <- paste("ENH", as.character(ehs), sep = "_")

  # get hi-tad data
  # I assume there are strictly 5 columns or 6 columns
  if(ncol(hic_interaction)==6){
    anchor.one <- GenomicRanges::GRanges(hic_interaction$chr1,
                                         IRanges::IRanges(start = hic_interaction$start1, end = hic_interaction$end1))
    anchor.two <- GenomicRanges::GRanges(hic_interaction$chr2,
                                         IRanges::IRanges(start = hic_interaction$start2, end = hic_interaction$end2))
    hic <- GenomicInteractions::GenomicInteractions(anchor.one, anchor.two)

  }else{
    anchor.one <- GenomicRanges::GRanges(hic_interaction$chr,
                                         IRanges::IRanges(start = hic_interaction$start1, end = hic_interaction$end1))
    anchor.two <- GenomicRanges::GRanges(hic_interaction$chr,
                                         IRanges::IRanges(start = hic_interaction$start2, end = hic_interaction$end2))
    hic <- GenomicInteractions::GenomicInteractions(anchor.one, anchor.two)
  }


  # annotation
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
  hg19.genes <- GenomicFeatures::genes(txdb)
  hg19.transcripts <- GenomicFeatures::transcriptsBy(txdb, by="gene")
  hg19.transcripts <- hg19.transcripts[ names(hg19.transcripts)  %in% unlist(hg19.genes$gene_id) ]
  hg19_refseq_promoters <- GenomicFeatures::promoters(hg19.transcripts,upstream = up_stream,
                                                      downstream = down_stream)

  seqnames_record <- GenomeInfoDb::seqnames(hg19_refseq_promoters)
  hg19_refseq_promoters <- unlist(hg19_refseq_promoters[S4Vectors::`%in%`(seqnames_record,
                                                                          c('chrX','chrY',paste('chr',1:22,sep='')))])
  hg19_refseq_promoters <- unique(hg19_refseq_promoters)
  gene_names <- biomaRt::getBM(attributes = c("hgnc_symbol",'entrezgene_id','ensembl_gene_id'),
                               filters = "entrezgene_id",
                               values = names(hg19_refseq_promoters), mart = mart)
  hg19_refseq_promoters$geneSymbol <- gene_names$hgnc_symbol[match(names(hg19_refseq_promoters),
                                                                   gene_names$entrezgene_id)]
  hg19_refseq_promoters$ensembl_id <- gene_names$ensembl_gene_id[match(names(hg19_refseq_promoters),
                                                                       gene_names$entrezgene_id)]

  if(gene_id_option=="hgnc_symbol"){
    names(hg19_refseq_promoters) <- hg19_refseq_promoters$geneSymbol
    hg19_refseq_promoters <- hg19_refseq_promoters[(!is.na(names(hg19_refseq_promoters)) )&
                                                     (names(hg19_refseq_promoters)!='')]
  }else if(gene_id_option=="ensembl_gene_id"){
    names(hg19_refseq_promoters) <- hg19_refseq_promoters$ensembl_id
    hg19_refseq_promoters <-hg19_refseq_promoters[!is.na(names(hg19_refseq_promoters))]
  }

  # if target gene is not string 'all' then
  if(target_genes != 'all'){
    hg19_refseq_promoters <- hg19_refseq_promoters[ names(hg19_refseq_promoters)
                                                    %in% target_genes]
  }


  if(typeof(ref_promoters) == typeof('all')){
    annotation_promoters <- hg19_refseq_promoters
  }else{
    ref_promoters_granges <-GenomicRanges::GRanges(
      seqnames = ref_promoters$chr,
      ranges = IRanges::IRanges(start = ref_promoters$start, end = ref_promoters$end))
    overlap_granges <- IRanges::findOverlapPairs(ref_promoters_granges,hg19_refseq_promoters,
                                                 type=link_type, ignore.strand=T)

    annotation_promoters <- overlap_granges@first
    names(annotation_promoters) <- names(overlap_granges@second)
    annotation_promoters <- unique(annotation_promoters)
  }


  annotation.features <- list(promoter = annotation_promoters, enhancer = ehs)
  GenomicInteractions::annotateInteractions(hic, annotation.features)
  interaction_index <- GenomicInteractions::isInteractionType(hic, "promoter", "enhancer")
  hic_subset <- hic[interaction_index]

  if(ncol(hic_interaction) == 6){
    df1 <- data.table::data.table(promoters = S4Vectors::mcols(GenomicInteractions::anchorOne(hic_subset))$promoter.id,
                                  promoter_chr = hic_interaction$chr1[interaction_index],
                                  enhs = S4Vectors::mcols(GenomicInteractions::anchorTwo(hic_subset))$enhancer.id,
                                  enh_chr = hic_interaction$chr2[interaction_index])
    df1 <- df1[ (!is.na(df1$promoters)) & (!is.na(df1$enhs)) ,]

    df2 <- data.table::data.table(promoters = S4Vectors::mcols(GenomicInteractions::anchorTwo(hic_subset))$promoter.id,
                                  promoter_chr = hic_interaction$chr2[interaction_index],
                                  enhs = S4Vectors::mcols(GenomicInteractions::anchorOne(hic_subset))$enhancer.id,
                                  enh_chr = hic_interaction$chr1[interaction_index]
    )
    df2 <- df2[ (!is.na(df2$promoters)) & (!is.na(df2$enhs)) ,]
  }else if( ncol(hic_interaction) == 5){
    df1 <- data.table::data.table(promoters = S4Vectors::mcols(GenomicInteractions::anchorOne(hic_subset))$promoter.id,
                                  promoter_chr = hic_interaction$chr[interaction_index],
                                  enhs = S4Vectors::mcols(GenomicInteractions::anchorTwo(hic_subset))$enhancer.id,
                                  enh_chr = hic_interaction$chr[interaction_index]
    )
    df1 <- df1[ (!is.na(df1$promoters)) & (!is.na(df1$enhs)) ,]

    df2 <- data.table::data.table(promoters = S4Vectors::mcols(GenomicInteractions::anchorTwo(hic_subset))$promoter.id,
                                  promoter_chr = hic_interaction$chr[interaction_index],
                                  enhs = S4Vectors::mcols(GenomicInteractions::anchorOne(hic_subset))$enhancer.id,
                                  enh_chr = hic_interaction$chr[interaction_index]
    )
    df2 <- df2[ (!is.na(df2$promoters)) & (!is.na(df2$enhs)) ,]
  }


  df <- data.table::rbindlist(list(df1,df2))
  df$id <- seq.int(nrow(df))
  df_enhancers <- df[,c('enhs','id')]

  df_enhancers <- df_enhancers[, .(enhancer = unlist(enhs)),by = id]
  df <- dplyr::left_join(df,df_enhancers,by = 'id')[,c('promoters','promoter_chr','enhancer','enh_chr')]
  df <- data.table::data.table(df)
  df$id <- seq.int(nrow(df))
  df_promoters <- df[,c('promoters','id')]
  df_promoters<- df_promoters[, .(promoter = unlist(promoters)),by = id]
  df <- dplyr::left_join(df,df_promoters,by = 'id')[,c('promoter','promoter_chr','enhancer','enh_chr')]
  df <- dplyr::distinct(df)
  sp_list  <- data.table::tstrsplit(df$enhancer,"_|:|-", keep = c(2,3,4))
  df$enh_start <- as.numeric(sp_list[[2]])
  df$enh_end <- as.numeric(sp_list[[3]])


  ref_df <- data.frame(promoter = names(annotation_promoters),
                       promoter_start = annotation_promoters@ranges@start,
                       promoter_end = annotation_promoters@ranges@start +
                         annotation_promoters@ranges@width, stringsAsFactors = F)


  final_df <- dplyr::full_join(ref_df,df,by = 'promoter')
  final_df <- na.omit(final_df)
  final_df <- final_df[,c('promoter','promoter_chr',
                          'promoter_start',
                          'promoter_end','enh_chr',
                          'enh_start','enh_end')]
  colnames(final_df) <- c('gene','gene_chr',
                          'promoter_start',
                          'promoter_end','enh_chr',
                          'enh_start','enh_end')

  return(final_df)
}









library(withr)


with_makevars(c(PKG_CFLAGS = "-std=c11"),
              install.packages("plyr", repos = NULL, type = "source"),
              assignment = "+=")


library(readxl)


library("DirichletMultinomial")
library("GenomicRanges")
library("GenomicFeatures")
library("GenomeInfoDb")
library("IRanges")
library("S4Vectors")
library("biomaRt")

library("gsl")
library("GenomicInteractions")

library("JASPAR2018")
library("BSgenome.Hsapiens.UCSC.hg19")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")

library("TFBSTools")
library("motifmatchr")
library("BSgenome.Hsapiens.UCSC.hg19")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")



print(paste0(":) Please note that the code for scGRN comes from: https://github.com/daifengwanglab/scGRN/tree/master/R"))
print("Please note that scGRN was written in Professor Daifeng Wang's Lab by Ting Jin, Peter Rehani, Mufang Ying, Jiawei Huang, Shuang Liu, Panagiotis Roussos & Professor Daifeng Wang")
print("loading in functions...")

scGRN_interaction = function(hic_interaction, enhancers, ref_promoters = 'all',up_stream = 2500,
                             down_stream = 2500, link_type = 'within' ,target_genes='all',
                             gene_id_option = 'hgnc_symbol',
                             mart = biomaRt::useMart(biomart="ENSEMBL_MART_ENSEMBL",
                                                     dataset="hsapiens_gene_ensembl",
                                                     host="uswest.ensembl.org")){


  ehs <- GenomicRanges::GRanges(seqnames = enhancers$chr, ranges = IRanges::IRanges(start = enhancers$start,
                                                                                    end = enhancers$end))
  names(ehs) <- paste("ENH", as.character(ehs), sep = "_")

  # get hi-tad data
  # I assume there are strictly 5 columns or 6 columns
  if(ncol(hic_interaction)==6){
    anchor.one <- GenomicRanges::GRanges(hic_interaction$chr1,
                                         IRanges::IRanges(start = hic_interaction$start1, end = hic_interaction$end1))
    anchor.two <- GenomicRanges::GRanges(hic_interaction$chr2,
                                         IRanges::IRanges(start = hic_interaction$start2, end = hic_interaction$end2))
    hic <- GenomicInteractions::GenomicInteractions(anchor.one, anchor.two)

  }else{
    anchor.one <- GenomicRanges::GRanges(hic_interaction$chr,
                                         IRanges::IRanges(start = hic_interaction$start1, end = hic_interaction$end1))
    anchor.two <- GenomicRanges::GRanges(hic_interaction$chr,
                                         IRanges::IRanges(start = hic_interaction$start2, end = hic_interaction$end2))
    hic <- GenomicInteractions::GenomicInteractions(anchor.one, anchor.two)
  }


  # annotation
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
  hg19.genes <- GenomicFeatures::genes(txdb)
  hg19.transcripts <- GenomicFeatures::transcriptsBy(txdb, by="gene")
  hg19.transcripts <- hg19.transcripts[ names(hg19.transcripts)  %in% unlist(hg19.genes$gene_id) ]
  hg19_refseq_promoters <- GenomicFeatures::promoters(hg19.transcripts,upstream = up_stream,
                                                      downstream = down_stream)

  seqnames_record <- GenomeInfoDb::seqnames(hg19_refseq_promoters)
  hg19_refseq_promoters <- unlist(hg19_refseq_promoters[S4Vectors::`%in%`(seqnames_record,
                                                                          c('chrX','chrY',paste('chr',1:22,sep='')))])
  hg19_refseq_promoters <- unique(hg19_refseq_promoters)
  gene_names <- biomaRt::getBM(attributes = c("hgnc_symbol",'entrezgene_id','ensembl_gene_id'),
                               filters = "entrezgene_id",
                               values = names(hg19_refseq_promoters), mart = mart)
  hg19_refseq_promoters$geneSymbol <- gene_names$hgnc_symbol[match(names(hg19_refseq_promoters),
                                                                   gene_names$entrezgene_id)]
  hg19_refseq_promoters$ensembl_id <- gene_names$ensembl_gene_id[match(names(hg19_refseq_promoters),
                                                                       gene_names$entrezgene_id)]

  if(gene_id_option=="hgnc_symbol"){
    names(hg19_refseq_promoters) <- hg19_refseq_promoters$geneSymbol
    hg19_refseq_promoters <- hg19_refseq_promoters[(!is.na(names(hg19_refseq_promoters)) )&
                                                     (names(hg19_refseq_promoters)!='')]
  }else if(gene_id_option=="ensembl_gene_id"){
    names(hg19_refseq_promoters) <- hg19_refseq_promoters$ensembl_id
    hg19_refseq_promoters <-hg19_refseq_promoters[!is.na(names(hg19_refseq_promoters))]
  }

  # if target gene is not string 'all' then
  if(target_genes != 'all'){
    hg19_refseq_promoters <- hg19_refseq_promoters[ names(hg19_refseq_promoters)
                                                    %in% target_genes]
  }


  if(typeof(ref_promoters) == typeof('all')){
    annotation_promoters <- hg19_refseq_promoters
  }else{
    ref_promoters_granges <-GenomicRanges::GRanges(
      seqnames = ref_promoters$chr,
      ranges = IRanges::IRanges(start = ref_promoters$start, end = ref_promoters$end))
    overlap_granges <- IRanges::findOverlapPairs(ref_promoters_granges,hg19_refseq_promoters,
                                                 type=link_type, ignore.strand=T)

    annotation_promoters <- overlap_granges@first
    names(annotation_promoters) <- names(overlap_granges@second)
    annotation_promoters <- unique(annotation_promoters)
  }


  annotation.features <- list(promoter = annotation_promoters, enhancer = ehs)
  GenomicInteractions::annotateInteractions(hic, annotation.features)
  interaction_index <- GenomicInteractions::isInteractionType(hic, "promoter", "enhancer")
  hic_subset <- hic[interaction_index]

  if(ncol(hic_interaction) == 6){
    df1 <- data.table::data.table(promoters = S4Vectors::mcols(GenomicInteractions::anchorOne(hic_subset))$promoter.id,
                                  promoter_chr = hic_interaction$chr1[interaction_index],
                                  enhs = S4Vectors::mcols(GenomicInteractions::anchorTwo(hic_subset))$enhancer.id,
                                  enh_chr = hic_interaction$chr2[interaction_index])
    df1 <- df1[ (!is.na(df1$promoters)) & (!is.na(df1$enhs)) ,]

    df2 <- data.table::data.table(promoters = S4Vectors::mcols(GenomicInteractions::anchorTwo(hic_subset))$promoter.id,
                                  promoter_chr = hic_interaction$chr2[interaction_index],
                                  enhs = S4Vectors::mcols(GenomicInteractions::anchorOne(hic_subset))$enhancer.id,
                                  enh_chr = hic_interaction$chr1[interaction_index]
    )
    df2 <- df2[ (!is.na(df2$promoters)) & (!is.na(df2$enhs)) ,]
  }else if( ncol(hic_interaction) == 5){
    df1 <- data.table::data.table(promoters = S4Vectors::mcols(GenomicInteractions::anchorOne(hic_subset))$promoter.id,
                                  promoter_chr = hic_interaction$chr[interaction_index],
                                  enhs = S4Vectors::mcols(GenomicInteractions::anchorTwo(hic_subset))$enhancer.id,
                                  enh_chr = hic_interaction$chr[interaction_index]
    )
    df1 <- df1[ (!is.na(df1$promoters)) & (!is.na(df1$enhs)) ,]

    df2 <- data.table::data.table(promoters = S4Vectors::mcols(GenomicInteractions::anchorTwo(hic_subset))$promoter.id,
                                  promoter_chr = hic_interaction$chr[interaction_index],
                                  enhs = S4Vectors::mcols(GenomicInteractions::anchorOne(hic_subset))$enhancer.id,
                                  enh_chr = hic_interaction$chr[interaction_index]
    )
    df2 <- df2[ (!is.na(df2$promoters)) & (!is.na(df2$enhs)) ,]
  }


  df <- data.table::rbindlist(list(df1,df2))
  df$id <- seq.int(nrow(df))
  df_enhancers <- df[,c('enhs','id')]

  df_enhancers <- df_enhancers[, .(enhancer = unlist(enhs)),by = id]
  df <- dplyr::left_join(df,df_enhancers,by = 'id')[,c('promoters','promoter_chr','enhancer','enh_chr')]
  df <- data.table::data.table(df)
  df$id <- seq.int(nrow(df))
  df_promoters <- df[,c('promoters','id')]
  df_promoters<- df_promoters[, .(promoter = unlist(promoters)),by = id]
  df <- dplyr::left_join(df,df_promoters,by = 'id')[,c('promoter','promoter_chr','enhancer','enh_chr')]
  df <- dplyr::distinct(df)
  sp_list  <- data.table::tstrsplit(df$enhancer,"_|:|-", keep = c(2,3,4))
  df$enh_start <- as.numeric(sp_list[[2]])
  df$enh_end <- as.numeric(sp_list[[3]])


  ref_df <- data.frame(promoter = names(annotation_promoters),
                       promoter_start = annotation_promoters@ranges@start,
                       promoter_end = annotation_promoters@ranges@start +
                         annotation_promoters@ranges@width, stringsAsFactors = F)


  final_df <- dplyr::full_join(ref_df,df,by = 'promoter')
  final_df <- na.omit(final_df)
  final_df <- final_df[,c('promoter','promoter_chr',
                          'promoter_start',
                          'promoter_end','enh_chr',
                          'enh_start','enh_end')]
  colnames(final_df) <- c('gene','gene_chr',
                          'promoter_start',
                          'promoter_end','enh_chr',
                          'enh_start','enh_end')

  return(final_df)
}

library("stringr")
scGRN_getTF <- function(df, database = JASPAR2018::JASPAR2018, species_type = 9606, min_score = 0.9,
                        pwm_type = 'prob',num_cores = 2){

  opts <- list()
  opts[["species"]] <- species_type
  opts[["all_versions"]] <- TRUE
  PFMatrixList <- TFBSTools::getMatrixSet(database, opts)
  pwmlist <- TFBSTools::toPWM(PFMatrixList, type = pwm_type)
  TF_names <- TFBSTools::name(pwmlist)
  names(TF_names) = NULL

  TF_names_splited = sapply(TF_names,data.table::tstrsplit,'::|\\(var.2\\)|\\(var.3\\)')

  df$promoter_id <- paste(df$gene_chr,':',df$promoter_start,'-',df$promoter_end,sep = '')
  df$enhancer_id <- paste(df$enh_chr,':',df$enh_start,'-',df$enh_end,sep = '')
  df = data.table::data.table(df)

  df_p <- dplyr::distinct(df[,c('gene_chr','promoter_start','promoter_end','promoter_id')])
  df_e <- dplyr::distinct(df[,c('enh_chr','enh_start','enh_end','enhancer_id')])

  suppressWarnings( G1 <- GenomicRanges::GRanges(seqnames = df_p$gene_chr,
                                                 IRanges::IRanges(start=df_p$promoter_start,
                                                                  end=df_p$promoter_end),
                                                 seqlengths = GenomeInfoDb::seqlengths(BSgenome.Hsapiens.UCSC.hg19::Hsapiens)[1:24]))
  G1 <- GenomicRanges::trim(G1)
  suppressWarnings( G2 <- GenomicRanges::GRanges(seqnames = df_e$enh_chr,
                                                 IRanges::IRanges(start=df_e$enh_start,
                                                                  end=df_e$enh_end),
                                                 seqlengths = GenomeInfoDb::seqlengths(BSgenome.Hsapiens.UCSC.hg19::Hsapiens)[1:24]))
  G2 <- GenomicRanges::trim(G2)


  cl <- parallel::makeCluster(num_cores) # not overload your computer
  doParallel::registerDoParallel(cl)

  `%dopar%` <- foreach::`%dopar%`
  df_p$promoter_TF <- foreach::foreach(i = 1:nrow(df_p), .combine = rbind,
                                       .packages = c('data.table','motifmatchr')) %dopar% {
                                         peak <- G1[i]
                                         motif_ix <- matchMotifs(pwmlist, peak,
                                                                 genome = "hg19",
                                                                 out = "scores"
                                         )
                                         result <- motifScores(motif_ix)[1,]
                                         curr_TF <- unique(unlist(TF_names_splited[result > quantile(result,min_score)]))
                                         if(length(curr_TF) == 0){
                                           curr_TF <- NA
                                         }
                                         data.table(promoter_TF = list(curr_TF))

                                       }
  parallel::stopCluster(cl)

  cl <- parallel::makeCluster(num_cores)
  doParallel::registerDoParallel(cl)
  df_e$enhancer_TF <- foreach::foreach(i = 1:nrow(df_e), .combine = rbind,
                                       .packages = c('data.table','motifmatchr')) %dopar% {
                                         peak <- G2[i]
                                         motif_ix <- matchMotifs(pwmlist, peak,
                                                                 genome = "hg19",
                                                                 out = "scores")
                                         result <- motifScores(motif_ix)[1,]
                                         curr_TF <- unique(unlist(TF_names_splited[result > quantile(result,min_score)]))
                                         if(length(curr_TF) == 0){
                                           curr_TF <- NA
                                         }
                                         data.table(promoter_TF = list(curr_TF))
                                       }
  parallel::stopCluster(cl)


  df$promoter_TF <- df_p$promoter_TF[match(df$promoter_id, df_p$promoter_id)]
  df$enhancer_TF <- df_e$enhancer_TF[match(df$enhancer_id, df_e$enhancer_id)]

  df <- df[, c('gene','promoter_id','enhancer_id',
               'promoter_TF','enhancer_TF')]
  colnames(df) <- c('gene','promoter','enhancer',
                    'promoter_TF','enhancer_TF')
  return(df)

}


scGRN_getNt <- function(df, gexpr, df_gene_id = 'hgnc_symbol', gexpr_gene_id = 'hgnc_symbol',
                        cutoff_by = 'quantile', cutoff_percentage = 0.9, cutoff_absolute = 0.1,scaleby = 'no',
                        train_ratio = 0.7, num_cores = 2,
                        mart = biomaRt::useMart(biomart="ENSEMBL_MART_ENSEMBL",
                                                dataset="hsapiens_gene_ensembl",
                                                host="uswest.ensembl.org")){

  # here df is actually a data.table and is previously cleaned(no NA value)
  # typically we can generate df by calling function gremen_getTF.R
  # expression_data is a data.frame : rownames are gene names and each column
  #                                   is a observation

  # data preprocessing and unlist the nested TFs
  # scaleby can be one of c('gene', 'sample', 'no')
  # one can use mart = useMart(biomart="ENSEMBL_MART_ENSEMBL",
  # dataset="hsapiens_gene_ensembl",
  # host="uswest.ensembl.org")
  # cutoff_by can be either absolute or quantile

  # Get all the TFs
  cl <- parallel::makeCluster(num_cores) # not overload your computer
  doParallel::registerDoParallel(cl)
  `%dopar%` <- foreach::`%dopar%`

  df$TFs <- foreach::foreach(i = 1:nrow(df), .combine = rbind,
                             .packages = c('data.table')) %dopar% {
                               curr_TF <- unique(c(df$enhancer_TF[[i]], df$promoter_TF[[i]]))
                               curr_TF <- curr_TF[!is.na(curr_TF)]
                               if(length(curr_TF) == 0){
                                 curr_TF <- NA
                               }
                               data.table::data.table(TFs = list(curr_TF))
                             }
  parallel::stopCluster(cl)

  df <- df[, c('gene','enhancer','promoter','enhancer_TF','promoter_TF','TFs')]
  df <- df[!is.na(df$TFs),]
  df$id <- seq.int(nrow(df))
  df_TF <- df[,c('TFs','id')][,.(TF = unlist(TFs)), by = id]
  df <- dplyr::left_join(df,df_TF,by = 'id')

  df$id <- NULL
  df$TFs <- NULL


  cl <- parallel::makeCluster(num_cores) # not overload your computer
  doParallel::registerDoParallel(cl)
  df$TFbs <- foreach::foreach(i = 1:nrow(df), .combine = rbind
  ) %dopar% {

    if(is.na(df$TF[i])){
      print("error")
    }

    if((df$TF[i] %in% df$enhancer_TF[[i]]) & (df$TF[i] %in% df$promoter_TF[[i]])){
      binding_site <- 'both'
    }else if(df$TF[i] %in% df$enhancer_TF[[i]]){
      binding_site <- 'enhancer'
    }else if(df$TF[i] %in% df$promoter_TF[[i]]){
      binding_site <- 'promoter'
    }
    binding_site
  }
  parallel::stopCluster(cl)

  df <- df[,c('gene','promoter','enhancer','TFbs','TF')]
  df <- df[df$TF != df$gene, ]


  ###### change TF names from hgnc_symbol to emsembl_id
  ###### delete the TF whose ensembl_id is NA

  if(gexpr_gene_id == 'hgnc_symbol'){
    if(df_gene_id != 'hgnc_symbol'){
      gene_names <- biomaRt::getBM(attributes = c("hgnc_symbol","ensembl_gene_id"), filters = "ensembl_gene_id",
                                   values = unique(df$gene), mart = mart)
      df$gene <- gene_names$hgnc_symbol[match(df$gene, gene_names$ensembl_gene_id)]
      df <- na.omit(df)
      df <- df[df$gene != '',]
    }
  }

  if(gexpr_gene_id == 'ensembl_gene_id'){
    if(df_gene_id == 'ensembl_gene_id'){
      gene_names <- biomaRt::getBM(attributes = c("hgnc_symbol","ensembl_gene_id"), filters = "hgnc_symbol",
                                   values = unique(df$TF), mart = mart)
      df$TF <- gene_names$ensembl_gene_id[match(df$TF, gene_names$hgnc_symbol)]
      df <- na.omit(df)
    }else{
      gene_names <- biomaRt::getBM(attributes = c("hgnc_symbol","ensembl_gene_id"), filters = "hgnc_symbol",
                                   values = c(unique(df$gene),unique(df$TF)), mart = mart)
      df$gene <- gene_names$ensembl_gene_id[match(df$gene, gene_names$hgnc_symbol)]
      df$TF <- gene_names$ensembl_gene_id[match(df$TF, gene_names$hgnc_symbol)]
      df <- na.omit(df)
    }
  }

  tgs <- unique(df$gene)
  cl <- parallel::makeCluster(num_cores)
  # not overload your computer
  doParallel::registerDoParallel(cl)
  output_df <- foreach::foreach(i = 1:length(tgs), .combine = rbind, .packages='glmnet') %dopar% {

    selgene <- tgs[i]

    if(selgene %in% rownames(gexpr)){
      selTFs <- unique(df[df$gene == selgene, 'TF'])
      train_cols <- sample(1:ncol(gexpr), round(train_ratio*ncol(gexpr)))
      if(sum(rownames(gexpr) %in% selTFs) > 1 ){

        x.train <- t(gexpr[rownames(gexpr) %in% selTFs,train_cols])
        x.test <- t(gexpr[rownames(gexpr) %in% selTFs,-train_cols])
        y.train <- t(gexpr[selgene,train_cols])
        y.test <- t(gexpr[selgene,-train_cols])

        yfit <- vector('list',11)
        yhat <- vector('list',length(yfit))
        mse <- rep(Inf,length(yfit))

        if(sd(y.train) == 0){
          curr_df <- data.frame(TG = NULL, TF = NULL,
                                coef = NULL)
        }else{

          for (j in 0:10) {
            # assign(paste("fit", i, sep=""), cv.glmnet(x.train, y.train,
            # type.measure="mse", alpha=i/10,family="gaussian"))
            yfit[[j+1]] <- glmnet::cv.glmnet(x.train, y.train, type.measure="mse",
                                             alpha=j/10,family="gaussian",standardize= F,intercept=T)
            yhat[[j+1]] <- as.numeric(predict(yfit[[j+1]], s=yfit[[j+1]]$lambda.1se, newx=x.test))
            mse[j+1] <- mean((y.test - yhat[[j+1]])^2)
          }
          fitcoef <- coef(yfit[[which.min(mse)]], s = "lambda.min")


          TF_coef <- as.matrix(fitcoef)
          TF_coef <- TF_coef[2:nrow(TF_coef),]
          if(cutoff_by == 'quantile'){
            TF_coef <- TF_coef[abs(TF_coef) > quantile(abs(TF_coef),1 - cutoff_percentage)]
          }else if(cutoff_by == 'absolute'){
            TF_coef <- TF_coef[abs(TF_coef) > cutoff_absolute]
          }else{
            print('cutoff_by can only take absolute or quantile')
          }

          if(length(TF_coef) > 0){
            curr_df <- data.frame(TG = rep(selgene,length(TF_coef)), TF = names(TF_coef),
                                  coef = unname(TF_coef),stringsAsFactors = F)
          }else{
            curr_df <- data.frame(TG = NULL, TF = NULL,
                                  coef = NULL)
          }

        }
      }else{
        curr_df = data.frame(TG = NULL,TF = NULL, coef = NULL)
      }

    }else{
      curr_df = data.frame(TG = NULL,TF = NULL, coef = NULL)
    }

    curr_df
  }

  parallel::stopCluster(cl)

  # df and output_df
  # df has gene, enhancer, TF
  # output_df has TG, TF, coef

  # data cleaning
  if(nrow(output_df) > 0){
    df$id <- paste(df$gene,'-',df$TF,sep = '')
    output_df$id <- paste(output_df$TG,'-',output_df$TF,sep = '')
    output_df$TF <- NULL
    df <- dplyr::full_join(df,output_df,by = 'id')

    df <- na.omit(df)
    df <- df[,c('TG','TF','enhancer','promoter','TFbs','coef')]
  }else{
    df <- data.frame(TG = NULL, TF = NULL, enhancer = NULL, promoter = NULL, TFbs = NULL,
                     coef = NULL)
  }

  return(df)

}

print("... :) Please note we are done loading in scGRN :)")













#### Please note that here, we call these scgrnom functions
# please note that this gives us the interaction information from scGRNom
chromosome_scGRNom_interactionInfo <- function(chromName = "chr1", outputName = chromatinInteractionFilePath,
                                             hic_interaction, enhancers, bodyRegion){
  interactome_data = as_tibble(hic_interaction[which(hic_interaction[,1] == chromName),])
  enhancers = as_tibble(enhancer[which(enhancer[,1] == chromName),])
  df <- scGRN_interaction(interactome_data, enhancers)
  print(paste0(":) Please note the chromosome name is: ", chromName))
  fileOut = paste0(outputName, "scgrnom_",bodyRegion,"ChromatinInteractions_",chromName, ".csv")
  write.csv(df, fileOut)
  fileOutR = paste0(outputName, "scgrnom_",bodyRegion,"ChromatinInteractions_",chromName, ".RData")
  info = "ScGRNom Interations Between Promoters and Enhancers"
  save(df, info, chromName, bodyRegion, file = fileOutR)

  return(df)
}



pleaseGetDataFrameOfChromosomeRegionInfo <- function(regionVector, regulatoryElementInfo = ""){
  #regionVector
  # please create a function to return dataframe of the chromosome, regulatory region, start and end
  # from chr#:start-end
  # regionVector is a vector of  chr#:start-end values.  For example:
  # [1] "chr20:43277877-43282877" "chr20:43277877-43282877" "chr20:43277877-43282877"
  # [4] "chr20:43277877-43282877" "chr20:43277877-43282877" "chr20:43277877-43282877"
  #regulatoryElementInfo is the info on the regulatoryElement and if it is an enhancer, or promoter, etc.

  chromAndRegionsVec = unlist(strsplit(regionVector, ":"))

  chromVec = chromAndRegionsVec[seq(1, length(chromAndRegionsVec), by = 2)] # please take every other element as the chromosome info
  startAndEndVec = chromAndRegionsVec[seq(2, length(chromAndRegionsVec), by = 2)]
  startAndEndVec = unlist(strsplit(startAndEndVec, "-"))

  startVec = startAndEndVec[seq(1, length(startAndEndVec), by = 2)] # please take every other element as the start info (odd # indices)
  endVec = startAndEndVec[seq(2, length(startAndEndVec), by = 2)] # please take every other element as the end info (even # indices)

  dataframeToReturn = data.frame(cbind(chromVec, startVec, endVec))

  #colnames(dataframeToReturn) = c("chromosome", "start", "end")
  if (regulatoryElementInfo == ""){
  colnames(dataframeToReturn) = c("chromosome", "start", "end")
  } else {
     colnames(dataframeToReturn) = c("chromosome", paste0(regulatoryElementInfo, "_start"), paste0(regulatoryElementInfo, "_end"))

  }

  dataframeToReturn$regulatoryRegion = rep(regulatoryElementInfo, nrow(dataframeToReturn))
  return(dataframeToReturn)
  }



pleaseProcessPromoterRowAfterScgrnomGetTF <- function(promoterRow, includeIndividualTFsInGroup){
  # please note that this function takes in a promoter row vector returned by scgrnom getTF
  # it then processes that row vector so that it is more organized
  # each promoter TF is a separate element in a vector
  # TFs with "-" in them are split into multiple TFs as well
  # for example
  #promoterRow = tfInteractions_df$promoter_TF[row]
  promoterRow = str_replace_all(promoterRow, "c[(]", "")
  promoterRow = str_replace_all(promoterRow, "[)]", "")
  promoterRow = str_replace_all(promoterRow, " ", "")
  promoterRow = scan(what = "", text = promoterRow, sep = ",",  quiet = TRUE)
  promoterRow

  if (includeIndividualTFsInGroup){ # if this is true, we will also split upon the "-" and include multiple TFs there as indiviudal TFs as well
    listOfMultipleTFs_Prom = grep("-", promoterRow)
	for (j in 1:length(listOfMultipleTFs_Prom)){
      differentGenes = unlist(strsplit(promoterRow[listOfMultipleTFs_Prom[j]], "-"))
      for (k in 1:length(differentGenes)){
        promoterRow = c(promoterRow, differentGenes[k])
      }
    }
  }
  return(promoterRow)
}



pleaseProcessEnhancerRowAfterScgrnomGetTF_NEW <- function(enhancerRow, includeIndividualTFsInGroup){
  # please note that this function takes in a enhancer row vector returned by scgrnom getTF
  # it then processes that row vector so that it is more organized
  # each enhancer TF is a separate element in a vector
  # TFs with "-" in them are split into multiple TFs as well
  # for example
  #enhancerRow = tfInteractions_df$enhancer_TF[row]
  enhancerRow = str_replace_all(enhancerRow, "c[(]", "")
  enhancerRow = str_replace_all(enhancerRow, "[)]", "")
  enhancerRow = str_replace_all(enhancerRow, " ", "")
  enhancerRow = scan(what = "", text = enhancerRow, sep = ",", quiet = TRUE)
  enhancerRow

  if (includeIndividualTFsInGroup){ # if this is true, we will also split upon the "-" and include multiple TFs there as indiviudal TFs as well
    listOfMultipleTFs_Enh = grep("-", enhancerRow)

	for (j in 1:length(listOfMultipleTFs_Enh)){
      differentGenes = unlist(strsplit(enhancerRow[listOfMultipleTFs_Enh[j]], "-"))
      for (k in 1:length(differentGenes)){
        enhancerRow = c(enhancerRow, differentGenes[k])
      }
    }
  }
  return(enhancerRow)
}


###################
## Creating the interactions df


# Please note that this function returns the information on Promoters given a # of bases of interest for the promoter length.
# The default is 5,000 bases for the promoter length
pleaseGetPromoterInformation <- function(maxBasesUpstreamForPromoterLength = 5000){
  suppressMessages(library(TxDb.Hsapiens.UCSC.hg19.knownGene))
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  promoterDF = promoters(genes(txdb), upstream = maxBasesUpstreamForPromoterLength)
  promoterDF = data.frame(promoterDF)
  colnames(promoterDF) = c("chromosome", "promoter_start", "promoter_end", "width", "strand", "entrezID")
  promoterDF$regionOnDNA = rep("Promoter", nrow(promoterDF))
  # this is information on promoter interactions based on TxDb.Hsapiens.UCSC.hg19.knownGene
  return(promoterDF)
}



pleaseGetTFsFromPromoterAndInteractionDataStep2 <- function(interactionsDF, chromName, numOfCoresToUseForScgrnom,
                                                            bodyRegion, otherInfo, folderName, chromatinInteraction_Step2_FilePath,
                                                            chromatinInteraction_Step2_FilePath_InitialRDataObjects,
                                                            chromatinInteraction_Step2_FilePath_CleanedRDataObjects,
															chromatinInteraction_Step2_FilePath_CleanedCSVObjects,
															includeIndividualTFsInGroup = includeTheIndividualTFsInGroup, pleasePrintProgressAfterNumIterations_GetTF = 1000
                                                            #includeIndividualTFsInGroup = includeTheIndividualTFsInGroup
                                                            ){


  library(stringr)
  library(dplyr)
  library("tictoc")

  print(paste0("Please note that we are working on scGRN part 2 (getting the TFs) for: chromosome ",
               chromName, " (", bodyRegion, " ", otherInfo, ") :)"))

  # please convert to numeric to be safe:
interactionsDF[,grep("promoter_start", colnames(interactionsDF))] = as.numeric(interactionsDF[,grep("promoter_start", colnames(interactionsDF))])
interactionsDF[,grep("promoter_end", colnames(interactionsDF))] = as.numeric(interactionsDF[,grep("promoter_end", colnames(interactionsDF))])
interactionsDF[,grep("enh_start", colnames(interactionsDF))] = as.numeric(interactionsDF[,grep("enh_start", colnames(interactionsDF))])
interactionsDF[,grep("enh_end", colnames(interactionsDF))] = as.numeric(interactionsDF[,grep("enh_end", colnames(interactionsDF))])


  df_chrPart2 = scGRN_getTF(interactionsDF, num_cores = numOfCoresToUseForScgrnom)
  head(df_chrPart2)
  print(paste0("Please note that we have finished scGRN part 2 (getting the TFs) for: chromosome ",
               chromName, " (", bodyRegion, " ", otherInfo, ") :) and are saving our initial results: "))

  outPathFinalCSV_simple = paste0(chromatinInteraction_Step2_FilePath_CleanedCSVObjects, "UpdatedNetworkFor", chromName, ".csv")
  print(head(df_chrPart2))
  outPathInitialRData = paste0(chromatinInteraction_Step2_FilePath_InitialRDataObjects, "scgrnomGetTFs_", chromName, "_", folderName, "_initialResults.RData")
  outPathFinalRData = paste0(chromatinInteraction_Step2_FilePath_CleanedRDataObjects, "scgrnomGetTFs_", chromName, "_", folderName, "_finalResults.RData")
  outPathFinalCSV = paste0(chromatinInteraction_Step2_FilePath_CleanedCSVObjects, "scgrnomGetTFs_", chromName, "_", folderName, "_finalResults.csv")

  save(df_chrPart2, file = outPathInitialRData)

  tfInteractions_df <- apply(df_chrPart2, 2,as.character)
  #save(df_chrPart2,  tfInteractions_df, file = chromatinInteraction_Step2_FilePath_InitialRDataObjects)

  head(tfInteractions_df)
  tfInteractions_df = data.frame(tfInteractions_df)


  rateFor319561_seconds = 546
  if (nrow(tfInteractions_df) < 319561){
    estimatedTime_seconds = (nrow(tfInteractions_df)*rateFor319561_seconds)/(319561)
    print(paste0("Please note that the estimated time for this code to run will be approximately: ",
                 estimatedTime_seconds, " seconds! :("))
  } else {
  estimatedTime_mins = (nrow(tfInteractions_df)*rateFor319561_seconds)/(60.0*319561)
  print(paste0("Please note that the estimated time for this code to run will be approximately: ",
               estimatedTime_mins, " minutes! :("))
}
  #################################################
  targetGenesVec = as.vector(tfInteractions_df$gene)
  #organizedGetTFInfoForChromDF = data.frame()
  columnNamesGetTFOrganized = c("targetGene", "TF", "chromosome",
                                "start", "end", "regulatoryRegion")

  promoterRegionsVec = tfInteractions_df[,grep("\\bpromoter\\b", colnames(tfInteractions_df))] # please note that the promoter is column 2
  enhancerRegionsVec = tfInteractions_df[,grep("\\benhancer\\b", colnames(tfInteractions_df))]
  head(promoterRegionsVec)
  head(enhancerRegionsVec)

  promoterRegionDF_ALL = pleaseGetDataFrameOfChromosomeRegionInfo(promoterRegionsVec, "promoter")
  enhancerRegionDF_ALL = pleaseGetDataFrameOfChromosomeRegionInfo(enhancerRegionsVec, "enhancer")

  columnNamesGetTFOrganizedList <- list()
  counter = 1
  tic("running please :)")
  for (i in 1:nrow(tfInteractions_df)){
    if (i %% pleasePrintProgressAfterNumIterations_GetTF == 1){
      print(paste0(":) Please note i = ", i, " of ", nrow(tfInteractions_df)))
    }
    targetGene = targetGenesVec[i]
    row = i
    # please format promoterRow properly

    promoterRow = tfInteractions_df$promoter_TF[row]
    enhancerRow =  tfInteractions_df$enhancer_TF[row]
    promoterRow = pleaseProcessPromoterRowAfterScgrnomGetTF(promoterRow, includeIndividualTFsInGroup)
    enhancerRow = pleaseProcessEnhancerRowAfterScgrnomGetTF_NEW(enhancerRow, includeIndividualTFsInGroup)


    # target gene, TF, chromosome, start, end, regulatory region
    promoterRegionDF = promoterRegionDF_ALL[row,]
    enhancerRegionDF = enhancerRegionDF_ALL[row,]

    promTF_DF = cbind(targetGene, promoterRow, promoterRegionDF$chromosome,
                      promoterRegionDF$promoter_start, promoterRegionDF$promoter_end,
                      promoterRegionDF$regulatoryRegion)

    enhTF_DF = cbind(targetGene, enhancerRow, enhancerRegionDF$chromosome, 
                     enhancerRegionDF$enhancer_start,
                     enhancerRegionDF$enhancer_end,enhancerRegionDF$regulatoryRegion)



    columnNamesGetTFOrganizedList[[counter]] = promTF_DF
    counter = counter + 1
    columnNamesGetTFOrganizedList[[counter]] = enhTF_DF
    counter = counter + 1
  }
  print(paste0(":) Yay!  Please note we finished organizing the results for ", folderName, " and chromosome ", chromName))
organizedGetTFInfoForChromDF = do.call("rbind", columnNamesGetTFOrganizedList)
head(organizedGetTFInfoForChromDF)
organizedGetTFInfoForChromDF = unique(organizedGetTFInfoForChromDF)
colnames(organizedGetTFInfoForChromDF) = columnNamesGetTFOrganized
organizedGetTFInfoForChromDF = data.frame(organizedGetTFInfoForChromDF)
fileOut = paste0(chromatinInteraction_Step2_FilePath, bodyRegion,"_scgrnomChromatinGetTFs_",chromName, ".csv")
write.csv(organizedGetTFInfoForChromDF, fileOut)
fileOutR = paste0(chromatinInteraction_Step2_FilePath, bodyRegion,"_scgrnomChromatinGetTFs_",chromName, ".RData")
info = "ScGRNom Step 2 GetTFs Between Promoters and Enhancers"
save(organizedGetTFInfoForChromDF, info, chromName, bodyRegion, file = fileOutR)
save(organizedGetTFInfoForChromDF, info, chromName, bodyRegion, file = outPathFinalRData)
write.csv(data.frame(organizedGetTFInfoForChromDF), outPathFinalCSV)
write.csv(data.frame(organizedGetTFInfoForChromDF), outPathFinalCSV_simple)

print(head(organizedGetTFInfoForChromDF))
print(fileOutR)
print(toc())
return(organizedGetTFInfoForChromDF)
}


### OLD:gettingTheTFsForAChromosomeScGRNPart2

gettingTheTFsForAChromosomeScGRNPart2 <- function(chromosomeNum, dfToUse = miniDF,
                                                  bodyRegion = bodyRegion,
                                                  parentPath = part2outputName,
                                                  otherInfo = "",
                                                  numOfCoresToUse = 2
){
  print("******************************************************")
  print(paste0("Please note that we are working on scGRN part 2 (getting the TFs) for: chromosome ",
               chromosomeNum, " (", bodyRegion, " ", otherInfo, ") :)"))

  outputNamesVec = outputFileNameGetTF(parentPath, chromosomeNum, bodyRegion, otherInfo)
  interactDFForChromosome = getInteractionDataForAChromosome(dfToUse, chromosomeNum)
  df_chrPart2 = scGRN_getTF(interactDFForChromosome, num_cores = numOfCoresToUse)
  save(df_chrPart2, chromosomeNum, file = outputNamesVec[3])
  #save(df2_chr1, file = "df2_chr1.RData")
  df <- apply(df_chrPart2, 2,as.character)
  print(dim(df))
  if (nrow(df) > 0){
    rm(df_chrPart2)
    print(":) YAY! Please note we removed the huge object since we successfully converted this to a dataframe")
    print(paste0("Please note that there are: ", nrow(df), " in this dataframe for ",
                 chromosomeNum, " for part 2."))
  } else {
    print("Please note there was an error when we tried to turn the data object into a dataframe. :(")
    return(df_chrPart2) # we exit with error :(
  }

  View(df)

  head(df)
  write.csv(df, outputNamesVec[1])
  save(df, chromosomeNum, file = outputNamesVec[2])

  print("******************************************************")
  print(paste0("Please note that we have finished working on scGRN part 2 for: chromosome ", chromosomeNum, " :)"))
  print("******************************************************")
  print("******************************************************")
  print(head(df))
  return(df)
}

####


#### Please note that here, we call these scgrnom functions
# please note that this gives us the interaction information from Enhancers DF and Promoters Manually
chromosome_AdSNPhenoPromEnh_interactionInfo <- function(chromName = "chr1",
                                               outputName = chromatinInteraction_Step1ADSNPheno_FilePath,
                                               bodyRegion,
                                               maxBasesUpstreamForPromoterLength,
                                               enhancerChromatinFilePath,
                                               geneIDColumnNameEnhancerDF,
                                               enhancerRegionCol,
                                               pleaseOrganizeEnhancerRegions){

  interactionsDF_AllChrom = gettingFullPromAndEnhInteractionsFromEnhancerDataOnly(enhancerChromatinFilePath,
                                                                    geneIDColumnNameEnhancerDF,
                                                                    enhancerRegionCol,
                                                                    outputName,
                                                                    pleaseOrganizeEnhancerRegions)

  df = interactionsDF_AllChrom[which(interactionsDF_AllChrom$gene_chr == chromName),]
  # scGRN_interaction(interactome_data, enhancers)
  print(paste0(":) Please note the chromosome name is: ", chromName))
  fileOut = paste0(outputName, "adsnpheno_",bodyRegion,"ChromatInterac_",chromName, "_PromLen", maxBasesUpstreamForPromoterLength, "bp.csv")
  write.csv(df, fileOut)
  fileOutR = paste0(outputName, "adsnpheno_",bodyRegion,"ChromatInterac_",chromName,"_PromLen", maxBasesUpstreamForPromoterLength, "bp.RData")

  info = paste0("ADSNPheno Interations Between Promoters and Enhancers from Enhancers and Promoters ", maxBasesUpstreamForPromoterLength)
  save(df, info, chromName, bodyRegion, file = fileOutR)

  return(df)
}



gettingFullPromAndEnhInteractionsFromEnhancerDataOnly <- function(enhancerChromatinFilePath,
                                                                  geneIDColumnNameEnhancerDF,
                                                                  enhancerRegionCol,
                                                                  chromatinInteractionFilePath,
                                                                  pleaseOrganizeEnhancerRegions){
  print(paste0(":) please note that this function by Saniya gets the full promoters and enhancers dataframe of interactions only based on the "))
  print(paste0("enhancers dataframe and the idea that promoter lengths are a maximum of maxBasesUpstreamForPromoterLength = ", maxBasesUpstreamForPromoterLength, " base pairs upstream"))

  finalColNames = c("gene", "gene_chr", "promoter_start", "promoter_end", "enh_chr", "enh_start", "enh_end")
  print("Please note that the final interactions Dataframe is similar to what would be returned from HI-C Chromatin Interaction data in Step 1 of scgrnom and will be input for step 2: ")
  print(finalColNames)
  print("Above are the final column names for the interactions DF")
  enhancerInteractionsData <- read.csv(enhancerChromatinFilePath, header = TRUE)

  # please look for gene column there
  geneCol = unique(c(grep("Gene", colnames(enhancerInteractionsData)),
  grep("gene", colnames(enhancerInteractionsData))))[1]
  colnames(enhancerInteractionsData)[geneCol] = "gene"
  head(enhancerInteractionsData)
  # please update the gene_chr column:
  colnames(enhancerInteractionsData)[grep("chr", enhancerInteractionsData[1,])[1]] = "gene_chr"
  if (pleaseOrganizeEnhancerRegions){
    if (is.null(enhancerRegionCol)){
      enhancerRegionCol = intersect(grep(":",enhancerInteractionsData[1,]), grep("-",enhancerInteractionsData[1,]))
      print(paste0(":) Please note that since no enhancerRegionCol was given (it was null), then we thought that the enhancerRegionCol is in column: ",
                   enhancerRegionCol, " based on the location of : and - in the dataframe.  We thought that would be chr#:start-end for enhancer.  Please fix if this is wrong and specify non-null value for enhancerRegionCol otherwise!"))
    }

    # please note we look for : and - in the 1st row of the data to determine where the enhancer region column could be
    colnames(enhancerInteractionsData)[enhancerRegionCol] = "EnhancerRegion"
    #chr_col = setdiff(chr_col, enhancerRegionCol)[1] # please note we only take the column that does not have the

    print(paste0(":) Please note that we assume that the EnhancerRegion is in column: ", enhancerRegionCol))
    # please note that these are the first 5 rows of the enhancer data
    head(enhancerInteractionsData)

    enhancerRegion = enhancerInteractionsData[,enhancerRegionCol]
    enhancerDF = pleaseGetDataFrameOfChromosomeRegionInfo(enhancerRegion, "enhancer")
    colnames(enhancerDF)[1] = "gene_chr"
    enhancerInteractionsData = cbind(enhancerInteractionsData, enhancerDF)

    head(enhancerInteractionsData)
  }

  promoterDF = pleaseGetPromoterInformation(maxBasesUpstreamForPromoterLength = maxBasesUpstreamForPromoterLength)

  #geneIDColumnNameEnhancerDF = "entrezID" # what is the exact column name for the column with gene IDs in enhancer data? Please doublecheck if it is in the 1st row that the column name is correct.
  geneIDColumnEnhancers = grep(geneIDColumnNameEnhancerDF, colnames(enhancerInteractionsData))
  enhStartCol = c(grep("start", colnames(enhancerInteractionsData)),
                  grep("Start", colnames(enhancerInteractionsData)))[1]
  enhEndCol = c(grep("end", colnames(enhancerInteractionsData)),
                  grep("End", colnames(enhancerInteractionsData)))[1]
  colnames(enhancerInteractionsData)[enhStartCol] = "enh_start"
  colnames(enhancerInteractionsData)[enhEndCol] = "enh_end"

  promAndEnhancerInteractionsDF = merge(enhancerInteractionsData, promoterDF, by.x = geneIDColumnNameEnhancerDF, by.y = "entrezID")

  # please note that this is information on the interactions between promoters and enhancers in the Lateral Temporal Lobe
  head(promAndEnhancerInteractionsDF)

  colNumsVecEnh = c(grep("\\bgene\\b", colnames(promAndEnhancerInteractionsDF))[1],
  grep("gene_chr", colnames(promAndEnhancerInteractionsDF))[1],
  grep("promoter_start", colnames(promAndEnhancerInteractionsDF))[1],
  grep("promoter_end", colnames(promAndEnhancerInteractionsDF))[1],
  grep("gene_chr", colnames(promAndEnhancerInteractionsDF))[1],
  grep("\\benh_start\\b", colnames(promAndEnhancerInteractionsDF))[1],
  grep("\\benh_end\\b", colnames(promAndEnhancerInteractionsDF))[1])

  # dataFrame[rows, cols]
  interactionsDF_AllChrom = promAndEnhancerInteractionsDF[,colNumsVecEnh]

  colnames(interactionsDF_AllChrom) = finalColNames
  head(interactionsDF_AllChrom)

  interactionsDF_AllChrom$enh_start = as.numeric(interactionsDF_AllChrom$enh_start)
  interactionsDF_AllChrom$enh_end = as.numeric(interactionsDF_AllChrom$enh_end)

  fileOut = paste0(chromatinInteractionFilePath, "final_",bodyRegion,"PromEnhChromatInteract_PromLength_", maxBasesUpstreamForPromoterLength, "bp_ALL_Chrom", ".csv")
  print(paste0(":) Please note we printed out the promoters and enhancers interactions dataframe here: ",
               fileOut))
  write.csv(interactionsDF_AllChrom, fileOut)
  fileOutR = paste0(chromatinInteractionFilePath, "final_",bodyRegion,"PromEnhChromatInteract_PromLength_", maxBasesUpstreamForPromoterLength, "bp_ALL_Chrom", ".RData")
  save(interactionsDF_AllChrom, file = fileOutR)

  return(interactionsDF_AllChrom)
}

# please note that this method returns the interactionsDF we use
pleaseGetInteractionsDFForChromosome_haveChromatinRegulatoryNetwork_Done <- function(chromName, caseInfo){
  # haveChromatinInteractionRegulatoryNetwork = TRUE
  chromRegulatoryNetwork = read.csv(chromatinInteractionRegulatoryNetworkFilePath, header = TRUE)
  print(paste0("Please note the chromatin regulatory network: ", caseInfo, " and the chromosome interaction info is for: ", chromName))
  return(chromRegulatoryNetwork)
}


pleaseGetInteractionsDFForChromosome_haveEnhancerAndPromoterInteractionsData_Step1_Done <- function(chromName, caseInfo, enhAndPromoterInteractionsFilePath){
  # haveEnhancerAndPromoterInteractionsData = TRUE
  interactionsDF_ALL = read.csv(enhAndPromoterInteractionsFilePath, header = TRUE)

  interactionsDF = interactionsDF_ALL[which(interactionsDF_ALL[,chrCol] == chromName),]
  print(paste0("Please note the caseInfo: ", caseInfo, " and the chromosome interaction info is for: ", chromName))

  return(interactionsDF)
}



pleaseGetInteractionsDFForChromosome_haveHIC_PromoterAndEnhancerData <- function(chromName, chromatinInteraction_Step1_FilePath,
                                                                                 hic_interaction, enhancers, bodyRegion, caseInfo){
  interactionsDF = chromosome_scGRNom_interactionInfo(chromName = chromName,
                                                      outputName = chromatinInteraction_Step1_FilePath,
                                                      hic_interaction, enhancers, bodyRegion)
  print(paste0("Please note the caseInfo: ", caseInfo, " and the chromosome interaction info is for: ", chromName))

  return(interactionsDF)
}




pleaseGetInteractionsDFForChromosome_haveOnlyEnhancerData <- function(chromName, chromatinInteraction_Step1ADSNPheno_FilePath,
                                                                      bodyRegion,
                                                                      maxBasesUpstreamForPromoterLength,
                                                                      enhancerChromatinFilePath,
                                                                      geneIDColumnNameEnhancerDF,
                                                                      enhancerRegionCol,
                                                                      pleaseOrganizeEnhancerRegions, caseInfo){
  interactionsDF = chromosome_AdSNPhenoPromEnh_interactionInfo(chromName = chromName,
                                                               chromatinInteraction_Step1ADSNPheno_FilePath,
                                                               bodyRegion,
                                                               maxBasesUpstreamForPromoterLength,
                                                               enhancerChromatinFilePath,
                                                               geneIDColumnNameEnhancerDF,
                                                               enhancerRegionCol,
                                                               pleaseOrganizeEnhancerRegions)


  print(paste0("Please note the caseInfo: ", caseInfo, " and the chromosome interaction info is for: ", chromName))

  return(interactionsDF)
}


########################################
## GENE REGULATORY INFORMATION:

# please identify which geneExpressionSampleIds are TFs:
identifyingTFsInGeneExpressionDataSetOlder <- function(geneAndEntrezIDMappingDF,
                                                  allsamples, tfsDF){
  geneExpressionSamplesVec = row.names(allsamples)
  print(":) Please note that this function returns the list of Transcription Factors ")
  print("that are found in the gene expression data set by mapping on Common EntrezID")
  print("This avoids problems with different gene symbol names  :)")
  geneExpressionSamplesDF = data.frame(geneExpressionSamplesVec)

  # please note we can use our geneAndEntrezIDMappingDF to get the gene mappings from before
  geneDataFrame = data.frame(row.names(geneAndEntrezIDMappingDF), geneAndEntrezIDMappingDF[,1])
  colnames(geneDataFrame) = c("geneName", "entrezID")
  geneDataFrame$entrezID = as.character(geneDataFrame$entrezID)
  head(geneDataFrame)

  head(tfsDF)

  tfsInGeneExpressionData = inner_join(geneDataFrame,tfsDF,
                                       by="entrezID")#c("geneExpressionName" = "geneSymbolName"))#Customer.ID")


  geneExpressionIds = as.vector(geneDataFrame[,2])
  tfsInGeneExpressionData_Names = as.vector(tfsInGeneExpressionData$geneName)

  tfsFound = intersect(tfsInGeneExpressionData_Names, geneExpressionSamplesVec) #which(tfsInGeneExpressionData_Names %in% geneExpressionSamplesVec)
  print(paste0(":) Please note that ", length(tfsFound), " TFs are found in the gene expression data set"))

  geneSymbolsOfTFs = tfsFound
  return(geneSymbolsOfTFs)
}


identifyingTFsInGeneExpressionDataSetOlder2 <- function(geneAndEntrezIDMappingDF,
                                                  allsamples, tfsDF){
  geneExpressionSamplesVec = row.names(allsamples)
  print(":) Please note that this function returns the list of Transcription Factors ")
  print("that are found in the gene expression data set by mapping on Common EntrezID")
  print("This avoids problems with different gene symbol names  :)")
  geneExpressionSamplesDF = data.frame(geneExpressionSamplesVec)

  # please note we can use our geneAndEntrezIDMappingDF to get the gene mappings from before
  geneDataFrame = data.frame(row.names(geneAndEntrezIDMappingDF), geneAndEntrezIDMappingDF[,1])
  colnames(geneDataFrame) = c("geneName", "entrezID")
  geneDataFrame$entrezID = as.character(geneDataFrame$entrezID)
  head(geneDataFrame)

  head(tfsDF)

  tfsInGeneExpressionData_Part1 = inner_join(geneDataFrame,tfsDF,
                                       by="entrezID")#c("geneExpressionName" = "geneSymbolName"))#Customer.ID")
  colnames(tfsDF)
  otherTFsDF = tfsDF
  colnames(otherTFsDF) = colnames(geneDataFrame)
  head(otherTFsDF)
  tfsInGeneExpressionData_Part2 = inner_join(geneDataFrame,otherTFsDF,
                                             by="geneName")#c("geneExpressionName" = "geneSymbolName"))#Customer.ID")

  part1 = tfsInGeneExpressionData_Part1[,c(1,2)]
  part2 = tfsInGeneExpressionData_Part2[,c(1,2)]
  colnames(part1) = c("geneName", "entrezID")
  colnames(part2) = c("geneName", "entrezID")
  tfsInGeneExpressionData = unique(rbind(part1, part2))
  head(part1)
  dim(part1)
  dim(part2)
  dim(tfsInGeneExpressionData)
  head(tfsInGeneExpressionData_Part2)
  head(tfsInGeneExpressionData_Part1)
  geneExpressionIds = as.vector(geneDataFrame[,2])
  tfsInGeneExpressionData
  tfsInGeneExpressionData_Names = as.vector(tfsInGeneExpressionData$geneName)

  tfsFound = intersect(tfsInGeneExpressionData_Names, geneExpressionSamplesVec) #which(tfsInGeneExpressionData_Names %in% geneExpressionSamplesVec)
  print(paste0(":) Please note that ", length(tfsFound), " TFs are found in the gene expression data set"))

  geneSymbolsOfTFs = tfsFound
  return(geneSymbolsOfTFs)
}


identifyingTFsInGeneExpressionDataSetolderr <- function(geneAndEntrezIDMappingDF,
                                                  allsamples, tfsDF){
  geneExpressionSamplesVec = row.names(allsamples)
  print(":) Please note that this function returns the list of Transcription Factors ")
  print("that are found in the gene expression data set by mapping on Common EntrezID")
  print("This avoids problems with different gene symbol names  :)")
  geneExpressionSamplesDF = data.frame(geneExpressionSamplesVec)

  # please note we can use our geneAndEntrezIDMappingDF to get the gene mappings from before
  geneDataFrame = data.frame(row.names(geneAndEntrezIDMappingDF), geneAndEntrezIDMappingDF[,1])
  colnames(geneDataFrame) = c("geneName", "entrezID")
  geneDataFrame$entrezID = as.character(geneDataFrame$entrezID)
  head(geneDataFrame)

  head(tfsDF)

  tfsInGeneExpressionData_Part1 = inner_join(geneDataFrame,tfsDF,
                                             by="entrezID")#c("geneExpressionName" = "geneSymbolName"))#Customer.ID")
  head(tfsInGeneExpressionData_Part1)
  colnames(tfsDF)
  otherTFsDF = tfsDF
  colnames(otherTFsDF) = colnames(geneDataFrame)
  head(otherTFsDF)
  tfsInGeneExpressionData_Part2 = inner_join(geneDataFrame,otherTFsDF,
                                             by="geneName")#c("geneExpressionName" = "geneSymbolName"))#Customer.ID")

  part1 = tfsInGeneExpressionData_Part1[,c(1,2)]
  part2 = tfsInGeneExpressionData_Part2[,c(1,2)]

  tfsToUseVec = unique(c(part1[,1], part2[,1]))
   print(paste0(":) Please note that ", length(tfsToUseVec), " TFs are found in the gene expression data set"))

  geneSymbolsOfTFs = tfsToUseVec
  
  tfsToUseVec = str_replace_all(tfsToUseVec, "-", ".")

  geneSymbolsOfTFs = intersect(geneExpressionSamplesVec, tfsToUseVec) #tfsFound#tfsInGeneExpressionData_Names[tfsFound,] # unique(as.vector(geneExpressionSamplesVec[tfsFound]))
  print(paste0(":) Please note that now ", length(geneSymbolsOfTFs), " TFs are found in the gene expression data set"))



  return(geneSymbolsOfTFs)
}



identifyingTFsInGeneExpressionDataSet <- function(geneAndEntrezIDMappingDF,
                                                  allsamples, tfsDF){
  geneExpressionSamplesVec = row.names(allsamples)
  print(":) Please note that this function returns the list of Transcription Factors ")
  print("that are found in the gene expression data set by mapping on Common EntrezID")
  print("This avoids problems with different gene symbol names  :)")
  geneExpressionSamplesDF = data.frame(geneExpressionSamplesVec)

  # please note we can use our geneAndEntrezIDMappingDF to get the gene mappings from before
  geneDataFrame = data.frame(row.names(geneAndEntrezIDMappingDF), geneAndEntrezIDMappingDF[,1])
  colnames(geneDataFrame) = c("geneName", "entrezID")
  geneDataFrame$entrezID = as.character(geneDataFrame$entrezID)
  head(geneDataFrame)

  head(tfsDF)


  tfsInGeneExpressionData_Part1 = inner_join(geneDataFrame,tfsDF,
                                             by="entrezID")#c("geneExpressionName" = "geneSymbolName"))#Customer.ID")
  head(tfsInGeneExpressionData_Part1)
  colnames(tfsDF)
  otherTFsDF = tfsDF
  colnames(otherTFsDF) = colnames(geneDataFrame)
  head(otherTFsDF)
  tfsInGeneExpressionData_Part2 = inner_join(geneDataFrame,otherTFsDF,
                                             by="geneName")#c("geneExpressionName" = "geneSymbolName"))#Customer.ID")

  part1 = tfsInGeneExpressionData_Part1[,c(1,2)]
  part2 = tfsInGeneExpressionData_Part2[,c(1,2)]

  tfsToUseVec = unique(c(part1[,1], part2[,1]))
  print(paste0(":) Please note that ", length(tfsToUseVec), " TFs are found in the gene expression data set"))

  length(tfsToUseVec)
  geneSymbolsOfTFs = intersect(geneExpressionSamplesVec, tfsToUseVec) #tfsFound#tfsInGeneExpressionData_Names[tfsFound,] # unique(as.vector(geneExpressionSamplesVec[tfsFound]))
  if (length(geneSymbolsOfTFs) < length(tfsToUseVec)){
    differencesVec = setdiff(tfsToUseVec, geneSymbolsOfTFs)
    if(length(grep("-", differencesVec)) > 0){
      tfsToUseVec = str_replace_all(tfsToUseVec, "-", ".")
      geneSymbolsOfTFs = intersect(geneExpressionSamplesVec, tfsToUseVec) #tfsFound#tfsInGeneExpressionData_Names[tfsFound,] # unique(as.vector(geneExpressionSamplesVec[tfsFound]))

    } else if (length(grep("[.]", differencesVec)) > 0){
      tfsToUseVec = str_replace_all(tfsToUseVec, ".", "-")
      geneSymbolsOfTFs = intersect(geneExpressionSamplesVec, tfsToUseVec) #tfsFound#tfsInGeneExpressionData_Names[tfsFound,] # unique(as.vector(geneExpressionSamplesVec[tfsFound]))
    }
  }

  print(paste0(":) Please note that now ", length(geneSymbolsOfTFs), " TFs are found in the gene expression data set"))
  setdiff(tfsToUseVec, geneSymbolsOfTFs)


  return(geneSymbolsOfTFs)
}





findingTFsThatRegulateGeneModules <- function(clusterInfo, folderName, geneModMembershipOutputNameCSV,
                                              filePathOfTFtoModule, RTN_rDataPart,
                                              geneAndEntrezIDMappingDF, tfsDF, wgcnaWithKmeansAllObjectsExceptTOM,
											  enrichmentsForKeyModulesOnly = FALSE){#newRDataPart){
  load(RTN_rDataPart) #newRDataPart) # RTN info
  print(paste0(":) please note that next we can find the TFs that significantly regulate the gene modules:"))
  if (enrichmentsForKeyModulesOnly){
  # we only load data for the key modules then
  load(keyGeneAssignmentsWithKmeansOutputNameRData)
  #clusterInfo = clusterInfo[which(clusterInfo[,colNum] %in% keyModules),]

} else {
  # we load info on genes in all modules
  load(geneAssignmentsWithKmeansOutputNameRData)
}
  load(wgcnaWithKmeansAllObjectsExceptTOM)

  ############
  MEList = moduleEigengenes(datExpr, colors = dynamicColorsKmeans)
MEs = MEList$eigengenes
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
geneModuleMembership = merge(geneAndEntrezIDMappingDF, geneModuleMembership, by = 0)
gene_names = colnames(datExpr)
rownames(geneModuleMembership)=gene_names
geneModuleMembership = geneModuleMembership[,-1]
colnames(geneModuleMembership)[1] = "entrezID"
write.csv(geneModuleMembership, geneModMembershipOutputNameCSV)

###########################
  modMembershipCorrelationPhenotype = read.csv(geneModMembershipOutputNameCSV, header = TRUE)
  head(modMembershipCorrelationPhenotype)
  moduleNamesVec = str_replace_all(colnames(modMembershipCorrelationPhenotype), "ME", "")
  colnames(modMembershipCorrelationPhenotype) = moduleNamesVec
  mraAllModulesDF = data.frame()
  colnamesMRA = c("Regulon (TF)", paste0("Universe.Size (# of Genes in ", bodyRegion, " Gene Expression Data Set"),
                  "Regulon.Size",                        "Total.Hits (# of Genes in Module)",
                  "Expected.Hits",                       "Observed.Hits",
                  "Pvalue",                              "Adjusted.Pvalue",
                  "currentTF_entrezID",
                  "geneModule",
                  "MRA_forModule", "region", "info", "phenotype", "hits", "numGenesInModule")

  moduleNamesVec = unique(clusterInfo$module)

  tni.regulon.summary(rtni)
  MMP_Alzh = modMembershipCorrelationPhenotype
  for (i in 1:length(moduleNamesVec)){
    moduleName = moduleNamesVec[i]
    print(paste0(":) i = ", i, " and please note we are on: ", moduleName))

    head(MMP_Alzh)

    geneNames = as.vector(MMP_Alzh[,1])
    moduleColumn = grep(paste0("\\b", moduleName, "\\b"), colnames(MMP_Alzh))[1]
    MMP_moduleBlack = MMP_Alzh[,moduleColumn]#grep(paste0("\\b", moduleName, "\\b"), colnames(MMP_Alzh))[1]# MMP_Alzh$black
    names(MMP_moduleBlack) = geneNames
    head(MMP_moduleBlack)
    # dataFrame[rows, columns]  :)
    blackGenes = as.vector(as.character(MMP_Alzh[which(clusterInfo$module == moduleName), 1]))

    rtna_black =  tni2tna.preprocess(object = rtni,
                                     phenotype = MMP_moduleBlack,
                                     hits = blackGenes)#,

    rtna_black <- tna.mra(rtna_black)

    mra_black <- tna.get(rtna_black, what="mra")#, ntop = -1)
    head(mra_black)
    mra_blackDF = data.frame(mra_black)
    # please note we can use our geneAndEntrezIDMappingDF to get the gene mappings from before
    geneDataFrame = data.frame(row.names(geneAndEntrezIDMappingDF), geneAndEntrezIDMappingDF[,1])
    colnames(geneDataFrame) = c("Regulon", "entrezID") # please name it so we can match with the mra output from RTN

    mra_blackDF = inner_join(mra_black, geneDataFrame, by = "Regulon")
    head(mra_blackDF)
    mra_blackDF = cbind(mra_blackDF, rep(moduleName, nrow(mra_black)))

    mra_blackDF = cbind(mra_blackDF, rep(paste0("mra_", moduleName), nrow(mra_black)))
    mra_blackDF = cbind(mra_blackDF, rep(bodyRegion, nrow(mra_black)))

    mra_blackDF = cbind(mra_blackDF, rep("RTN_MasterRegulatoryAnalysis", nrow(mra_black)))
    # phenotype
    phenotypeInfo = paste0("Module Membership of all genes in ", bodyRegion, " for the ",moduleName, " gene module")
    mra_blackDF = cbind(mra_blackDF, rep(phenotypeInfo, nrow(mra_black)))
    # hits
    hitsInfo = paste0("The ", length(blackGenes), " genes that are actually assigned to the ", moduleName, " gene module in the ", bodyRegion) # Membership of all genes in ", bodyRegion, " for the ",moduleName, " gene module")
    mra_blackDF = cbind(mra_blackDF, rep(hitsInfo, nrow(mra_black)))
    # numGenesInModule
    mra_blackDF = cbind(mra_blackDF, rep(length(blackGenes), nrow(mra_black)))
    colnames(mra_blackDF) = colnamesMRA

    head(mra_blackDF)

    blackGenes

    bodyRegion

    outputPath = paste0(filePathOfTFtoModule, "rtn_masterRegulatoryAnalysis_TFs_and_", moduleName, "_module_", folderName, ".csv")

    write.csv(mra_blackDF, outputPath)
    if (i == 1){
      mraAllModulesDF = mra_blackDF
    } else {
      mraAllModulesDF = rbind(mraAllModulesDF, mra_blackDF)
    }

  }
  mraAllModulesDF = unique(mraAllModulesDF)
  outputPath = paste0(filePathOfTFtoModule, "rtn_masterRegulatoryAnalysis_TFs_and_ALL_", length(moduleNamesVec), "modules_", folderName, ".csv")
  write.csv(mraAllModulesDF, outputPath)
  print(head(mraAllModulesDF))
  head(mraAllModulesDF)
  return(mraAllModulesDF)
}



# please remove:
findingTFsThatRegulateGeneModulesOld <- function(clusterInfo, folderName, geneModMembershipOutputNameCSV, filePathOfTFtoModule, RTN_rDataPart){#newRDataPart){
  load(RTN_rDataPart) #newRDataPart) # RTN info
  print(paste0(":) please note that next we can find the TFs that significantly regulate the gene modules:"))


  modMembershipCorrelationPhenotype = read.csv(geneModMembershipOutputNameCSV, header = TRUE)
  head(modMembershipCorrelationPhenotype)
  moduleNamesVec = str_replace_all(colnames(modMembershipCorrelationPhenotype), "ME", "")
  colnames(modMembershipCorrelationPhenotype) = moduleNamesVec
  mraAllModulesDF = data.frame()
  colnamesMRA = c("Regulon (TF)", paste0("Universe.Size (# of Genes in ", bodyRegion, " Gene Expression Data Set"),
                  "Regulon.Size",                        "Total.Hits (# of Genes in Module)",
                  "Expected.Hits",                       "Observed.Hits",
                  "Pvalue",                              "Adjusted.Pvalue",
                  "geneModule",
                  "MRA_forModule", "region", "info", "phenotype", "hits", "numGenesInModule")

  moduleNamesVec = unique(clusterInfo$module)

  tni.regulon.summary(rtni)
  MMP_Alzh = modMembershipCorrelationPhenotype
  for (i in 1:length(moduleNamesVec)){
    moduleName = moduleNamesVec[i]
    print(paste0(":) i = ", i, " and please note we are on: ", moduleName))

    head(MMP_Alzh)

    #condition = "Alzheimers"
    View(MMP_Alzh)
    geneNames = as.vector(MMP_Alzh[,1])
    moduleColumn = grep(paste0("\\b", moduleName, "\\b"), colnames(MMP_Alzh))[1]
    MMP_moduleBlack = MMP_Alzh[,moduleColumn]#grep(paste0("\\b", moduleName, "\\b"), colnames(MMP_Alzh))[1]# MMP_Alzh$black
    names(MMP_moduleBlack) = geneNames
    head(MMP_moduleBlack)
    # dataFrame[rows, columns]  :)
    blackGenes = as.vector(as.character(MMP_Alzh[which(clusterInfo$module == moduleName), 1]))

    rtna_black =  tni2tna.preprocess(object = rtni,
                                     phenotype = MMP_moduleBlack,
                                     hits = blackGenes)#,

    rtna_black <- tna.mra(rtna_black)

    mra_black <- tna.get(rtna_black, what="mra")#, ntop = -1)
    head(mra_black)
    mra_blackDF = data.frame(mra_black)
    mra_blackDF = cbind(mra_blackDF, rep(moduleName, nrow(mra_black)))

    mra_blackDF = cbind(mra_blackDF, rep(paste0("mra_", moduleName), nrow(mra_black)))
    mra_blackDF = cbind(mra_blackDF, rep(bodyRegion, nrow(mra_black)))

    mra_blackDF = cbind(mra_blackDF, rep("RTN_MasterRegulatoryAnalysis", nrow(mra_black)))
    # phenotype
    phenotypeInfo = paste0("Module Membership of all genes in ", bodyRegion, " for the ",moduleName, " gene module")
    mra_blackDF = cbind(mra_blackDF, rep(phenotypeInfo, nrow(mra_black)))
    # hits
    hitsInfo = paste0("The ", length(blackGenes), " genes that are actually assigned to the ", moduleName, " gene module in the ", bodyRegion) # Membership of all genes in ", bodyRegion, " for the ",moduleName, " gene module")
    mra_blackDF = cbind(mra_blackDF, rep(hitsInfo, nrow(mra_black)))
    # numGenesInModule
    mra_blackDF = cbind(mra_blackDF, rep(length(blackGenes), nrow(mra_black)))
    colnames(mra_blackDF) = colnamesMRA

    head(mra_blackDF)

    blackGenes
    #"MRA_forModule", "region", "info", "phenotype", "hits"

    bodyRegion

    outputPath = paste0(filePathOfTFtoModule, "rtn_masterRegulatoryAnalysis_TFs_and_", moduleName, "_module_", folderName, ".csv")

    write.csv(mra_blackDF, outputPath)
    if (i == 1){
      mraAllModulesDF = mra_blackDF
    } else {
      mraAllModulesDF = rbind(mraAllModulesDF, mra_blackDF)
    }

  }
  mraAllModulesDF = unique(mraAllModulesDF)
  outputPath = paste0(filePathOfTFtoModule, "rtn_masterRegulatoryAnalysis_TFs_and_ALL_", length(moduleNamesVec), "modules_", folderName, ".csv")
  write.csv(mraAllModulesDF, outputPath)
  print(head(mraAllModulesDF))
  head(mraAllModulesDF)
  return(mraAllModulesDF)
}










# dynamically figure out what ye need to load in :)

pleaseGetDatTraitsPhenotypesFile <- function(geneExpressionOutputNameRData, phenotypesFilePath, datTraitsRDataFilePath){
	load(geneExpressionOutputNameRData) # please do this to get information on the datExpr (and geneEntrezIDMapping)
	diseaseTraits = read.csv(phenotypesFilePath,
							 header = TRUE)

	print(paste0(":) Please note that we have laoded in the phenotypesfile (diseaseTraits) from the phenotypesFilePath: ", phenotypesFilePath))
	# information on the phenotypes, matching your phenotypes file with the samples in the gene expression data! YAY!"))
	print(paste0(":) Please ensure that these results are consistent with what you need: ", dim(diseaseTraits)))

	print("The first 5 rows are:")
	print(head(diseaseTraits))
	head(diseaseTraits)
	colnames(diseaseTraits)

	str(diseaseTraits)

	# chr

	class(diseaseTraits[,2])

	typesToRemoveVec = c(1)
	traitsToKeepVec = c()
	print(":) We only keep the numeric traits here:")
	for (k in 3:ncol(diseaseTraits)){
	  if (class(diseaseTraits[,k]) == "character"){
		typesToRemoveVec = c(typesToRemoveVec, k)
	  } else {
		traitsToKeepVec = c(traitsToKeepVec, k)
	  }

	}
	typesToRemoveVec
	print(paste0(":) Please note that we remove these ", length(typesToRemoveVec), " phenotypes/traits as they were non-numeric:"))
	print(typesToRemoveVec)
	allTraits = diseaseTraits[, -typesToRemoveVec]

	#colnames(diseaseTraits)
	#allTraits = diseaseTraits[, -c(1)]#)], 3,  6, 12)]
	#Condition, clinical_braak
	#no gender clinical_braak
	dim(allTraits)
	head(allTraits)
	names(allTraits)
	# [1] "??.."              "Patient"          "Condition"        "HealthStage"
	# [5] "Age"              "Gender"           "GenderGroup"      "nft"
	# [9] "braak"            "mmse"             "pmi"              "clinical_braak"
	# [13] "HealthStageBraak" "difference"       "controlStage"     "initialStage"
	# [17] "moderateStage"    "severeStage"      "asymptomatic"     "mild.symptoms"
	# [21] "dementia"         "braak1"           "braak2"           "braak3"
	# [44] "braak4"           "braak5"           "braak6"
	alzhSamples = rownames(datExpr);
	alzhSamples = str_replace_all(alzhSamples, "X", "")
	alzhSamples = str_replace_all(alzhSamples, "[.]", "-")
	#traitRows = match(alzhSamples, allTraits[,1]);
	traitRows = match(alzhSamples, diseaseTraits[,2]);
	datTraits = allTraits[traitRows, -1];
	length(traitRows)
	traitRows = traitRows[which(!isNA(traitRows))]
	length(traitRows)

	rownames(datTraits) = diseaseTraits[traitRows, 1];


	print(head(datTraits))
	print(paste0(":) Please note that we have recovered information on the phenotypes, matching your phenotypes file with the samples in the gene expression data! YAY!"))
	print(paste0(":) Please note the results are: ", dim(datTraits)))
	print(paste0(":) Please note that we will analyze gene and/or gene co-expression module relationships with these ", ncol(datTraits), " numeric phenotypes"))
	print(colnames(datTraits))
	save(datTraits, file = datTraitsRDataFilePath)
	print(paste0(":) Please note that the RData for this datTraits object has been stored here: ", datTraitsRDataFilePath))

	return(datTraits)

}



saveWGCNAPowerFilePathsInformation <- function(powerRelatedFilePathsRData, powerEstimate,  folderName, performAdditionalKMeansStep,
tomFilePath, fullDiseaseName, bodyRegion, outputAddOn, allObjectsFilePath, genePhenoAssociationFilePath

){
# please note that this function saves the file path information of the file paths that are based on WGCNA powerEstimate Used


dissTOMPowerOutputNameCSV = paste0(tomFilePath, "//dissTOMPower_",netType, "_pow",powerEstimate, "_", fullDiseaseName, "_", bodyRegion, outputAddOn, ".csv")
print(paste0("dissTOMPowerOutputNameCSV: ", dissTOMPowerOutputNameCSV))
dissTOMPowerOutputNameRData = paste0(tomFilePath, "//dissTOMPower_",netType, "_pow",powerEstimate, "_", fullDiseaseName, "_", bodyRegion, outputAddOn, ".RData")
print(paste0("dissTOMPowerOutputNameRData: ", dissTOMPowerOutputNameRData))

tomPowerOutputNameCSV = paste0(tomFilePath, "//TOMPower_",netType, "_pow",powerEstimate, "_", fullDiseaseName, "_", bodyRegion, outputAddOn, ".csv")
print(paste0("tomPowerOutputNameCSV: ", tomPowerOutputNameCSV))
tomPowerOutputNameRData = paste0(tomFilePath, "//TOMPower_",netType, "_pow",powerEstimate, "_", fullDiseaseName, "_", bodyRegion, outputAddOn, ".RData")
print(paste0("tomPowerOutputNameRData: ", tomPowerOutputNameRData))

geneModMembershipOutputNameCSV = paste0(geneModuleMembershipFilePath, wgcnaToAdd,netType, "_pow",powerEstimate, "_geneModuleMembershipCorr", outputAddOn, ".csv")
print(paste0("geneModMembershipOutputNameCSV: ", geneModMembershipOutputNameCSV))
mmpOutputNameCSV = paste0(geneModuleMembershipFilePath, wgcnaToAdd,netType, "_pow",powerEstimate, "_MMPvalues", outputAddOn, ".csv")
print(paste0("mmpOutputNameCSV: ", mmpOutputNameCSV))

modTraitCorrelationOutputNameCSV = paste0(modulePhenoAssociationFilePath, wgcnaToAdd,netType, "_pow",powerEstimate, "_moduleTraitCorrs", outputAddOn, ".csv")
modTraitCorrPValOutputNameCSV = paste0(modulePhenoAssociationFilePath, wgcnaToAdd,netType, "_pow",powerEstimate, "_moduleTraitCorrPValue", outputAddOn, ".csv")

print(paste0("modTraitCorrelationOutputNameCSV: ", modTraitCorrelationOutputNameCSV))
print(paste0("modTraitCorrPValOutputNameCSV: ", modTraitCorrPValOutputNameCSV))

modTraitCorrAndCorrPValueOutputNameRData = paste0(modulePhenoAssociationFilePath, wgcnaToAdd,netType, "_pow",powerEstimate, "_moduleTraitCorrAndP", outputAddOn, ".RData")
signifModTraitCorrelationOutputNameCSV = paste0(modulePhenoAssociationFilePath, wgcnaToAdd,netType, "_pow",powerEstimate, "_StatSignifModulePhenoCorr_Thresh", pValueCutOffForSignificance, outputAddOn, ".csv")

posCorsOnlyModuleTraitsOutputNameCSV = paste0(modulePhenoAssociationFilePath, "wgcna_with_kmeans_Power", powerEstimate, "_StatisticallySignificantModuleTraitCorrelations_POSITIVEonly_.csv")



print(paste0("modTraitCorrAndCorrPValueOutputNameRData: ", modTraitCorrAndCorrPValueOutputNameRData))
print(paste0("signifModTraitCorrelationOutputNameCSV: ", signifModTraitCorrelationOutputNameCSV))

updatedModuleTraitPhenotypeDFFilePath = paste0(modulePhenoAssociationFilePath, wgcnaToAdd, bodyRegion, "_", powerEstimate, "_StatSignificantModulePhenoCorrs_CombinedInfo",outputAddOn,".csv")

print(paste0("updatedModuleTraitPhenotypeDFFilePath: ", updatedModuleTraitPhenotypeDFFilePath))


wgcnaWithKmeansAllObjectsExceptTOM = paste0(allObjectsFilePath, wgcnaToAdd, bodyRegion, "_", powerEstimate, "_ALLKeyObjectsExceptTOM",outputAddOn,".RData")
wgcnaWithKmeansAllObjectsIncludingTOM = paste0(allObjectsFilePath, wgcnaToAdd, bodyRegion, "_", powerEstimate, "_ALLKeyObjectsIncludingTOM",outputAddOn,".RData")

print(paste0("wgcnaWithKmeansAllObjectsExceptTOM: ", wgcnaWithKmeansAllObjectsExceptTOM))
print(paste0("wgcnaWithKmeansAllObjectsIncludingTOM: ", wgcnaWithKmeansAllObjectsIncludingTOM))

wgcnaWithKmeansAllObjectsExceptTOMWithDate = paste0(allObjectsFilePath, wgcnaToAdd, bodyRegion, "_", powerEstimate, "_ALLKeyObjectsExceptTOM",outputAddOnWithDate,".RData")
wgcnaWithKmeansAllObjectsIncludingTOMWithDate = paste0(allObjectsFilePath, wgcnaToAdd, bodyRegion, "_", powerEstimate, "_ALLKeyObjectsIncludingTOM",outputAddOnWithDate,".RData")

print(paste0("wgcnaWithKmeansAllObjectsExceptTOMWithDate: ", wgcnaWithKmeansAllObjectsExceptTOMWithDate))
print(paste0("wgcnaWithKmeansAllObjectsIncludingTOMWithDate: ", wgcnaWithKmeansAllObjectsIncludingTOMWithDate))


wgcnaWithKmeansGeneTraitOrganizedPath = paste0(genePhenoAssociationFilePath, wgcnaToAdd, bodyRegion, "_", powerEstimate, "_GenePhenoSigCorrsOrganized",outputAddOn,".csv")
print(paste0("wgcnaWithKmeansGeneTraitOrganizedPath: ", wgcnaWithKmeansGeneTraitOrganizedPath))


wgcnaWithKmeansGeneTraitSignificancePath = paste0(genePhenoAssociationFilePath, wgcnaToAdd, bodyRegion, "_", powerEstimate, "_GenePhenoSigCorrs",outputAddOn,".csv")
wgcnaWithKmeansGeneTraitSigPValuePath = paste0(genePhenoAssociationFilePath, wgcnaToAdd, bodyRegion, "_", powerEstimate, "_GenePhenoCorrPValue",outputAddOn,".csv")

print(paste0("wgcnaWithKmeansGeneTraitSignificancePath: ", wgcnaWithKmeansGeneTraitSignificancePath))
print(paste0("wgcnaWithKmeansGeneTraitSigPValuePath: ", wgcnaWithKmeansGeneTraitSigPValuePath))

wgcnaWithKmeansGeneModuleMembershipPath  = paste0(geneModuleMembershipFilePath, wgcnaToAdd, bodyRegion, "_", powerEstimate, "_GeneModuleMembershipCorrelation",outputAddOn,".csv")
wgcnaWithKmeansGeneModMemberPValuePath= paste0(geneModuleMembershipFilePath, wgcnaToAdd, bodyRegion, "_", powerEstimate, "_GeneModMembershipCorrPValue",outputAddOn,".csv")

print(paste0("wgcnaWithKmeansGeneModuleMembershipPath: ", wgcnaWithKmeansGeneModuleMembershipPath))
print(paste0("wgcnaWithKmeansGeneModMemberPValuePath: ", wgcnaWithKmeansGeneModMemberPValuePath))

wgcnaWithKmeansGeneTraitsEdgeTablePathCSV  = paste0(genePhenoAssociationFilePath, wgcnaToAdd, bodyRegion, "_", powerEstimate, "_GeneTraitCorrForP",outputAddOn,".csv")
wgcnaWithKmeansGeneTraitsEdgeTablePathRData  = paste0(genePhenoAssociationFilePath, wgcnaToAdd, bodyRegion, "_", powerEstimate, "_GeneTraitCorrForP",outputAddOn,".RData")
signifGeneTraitCorrelationOutputNameCSV = paste0(genePhenoAssociationFilePath, wgcnaToAdd, bodyRegion, "_", powerEstimate, "_SigGeneTraitCorrForP", pValueCutOffForSignificance,outputAddOn,".csv")
signifPositiveGeneTraitCorrelationOutputNameCSV = paste0(genePhenoAssociationFilePath, wgcnaToAdd, bodyRegion, "_", powerEstimate, "_SigPosGenePhenoCorrForP", pValueCutOffForSignificance,outputAddOn,".csv")
updatedGeneTraitPhenotypeDFFilePath = paste0(genePhenoAssociationFilePath, wgcnaToAdd, bodyRegion, "_", powerEstimate, "finalCombo_SigGenePhenoCorrForP", pValueCutOffForSignificance,outputAddOn,".csv")

print(paste0("wgcnaWithKmeansGeneTraitsEdgeTablePathCSV: ", wgcnaWithKmeansGeneTraitsEdgeTablePathCSV))
print(paste0("wgcnaWithKmeansGeneTraitsEdgeTablePathRData: ", wgcnaWithKmeansGeneTraitsEdgeTablePathRData))
print(paste0("signifGeneTraitCorrelationOutputNameCSV: ", signifGeneTraitCorrelationOutputNameCSV))
print(paste0("signifPositiveGeneTraitCorrelationOutputNameCSV: ", signifPositiveGeneTraitCorrelationOutputNameCSV))
print(paste0("updatedGeneTraitPhenotypeDFFilePath: ", updatedGeneTraitPhenotypeDFFilePath))


############################################################################################

wgcnaPowerOutputNameRData = paste0(wgcnaFilePath, "//wgcnaInitialResultsObjects_",netType, "_pow",powerEstimate, "_", fullDiseaseName, "_", bodyRegion, outputAddOn, ".RData")

moduleCountsOutputName = paste0(wgcnaFilePath, "//WGCNAInitialModuleCounts_",netType, "_pow",powerEstimate, outputAddOn, ".csv")

wgcnaOutputName = paste0(wgcnaFilePath, "//WGCNAInitialResultsObjects_",netType, "_pow",powerEstimate, outputAddOn, ".RData")

#############################################


MEsWithKmeansOutputName = paste0(MEsFilePath, wgcnaToAdd, netType, "_pow",powerEstimate, "_ModuleEigengenesMEs", outputAddOn, ".csv")
geneAssignmentsWithKmeansOutputNameCSV = paste0(moduleAssignmentsFilePath, wgcnaToAdd, netType, "_pow",powerEstimate, "_ModuleAssignments", outputAddOn, ".csv")
geneAssignmentsWithKmeansOutputNameRData = paste0(moduleAssignmentsFilePath, wgcnaToAdd,netType, "_pow",powerEstimate, "_ModuleAssignments", outputAddOn, ".RData")



geneAssignmentsModAssignmentsPhenosWithKmeansOutputNameCSV = paste0(moduleAssignmentsFilePath, wgcnaToAdd, netType, "_pow",powerEstimate, "_ModuleAndGeneAssignAndPhenosSUMMARY", outputAddOn, ".csv")


keyGeneAssignmentsWithKmeansOutputNameCSV = paste0(moduleAssignmentsFilePath, wgcnaToAdd, netType, "_pow",powerEstimate, "_KeyModuleAssigns", outputAddOn, ".csv")
keyGeneAssignmentsWithKmeansOutputNameRData = paste0(moduleAssignmentsFilePath, wgcnaToAdd,netType, "_pow",powerEstimate, "_KeyModuleAssigns", outputAddOn, ".RData")

moduleCountsWithKmeansOutputNameCSV = paste0(pathForWGCNA, wgcnaToAdd,netType, "_pow",powerEstimate, "_ModuleCounts", outputAddOn, ".csv")



hubsCSV = paste0(pathForWGCNA, wgcnaToAdd,netType, "_pow",powerEstimate, "_hubGenes", outputAddOn, ".csv")

if (performAdditionalKMeansStep){
 allNetworkModulesCombinedTomSummaryStatisticsFileName = paste0(wgcnaAndKMeansOutputPath, "//", bodyRegion, "_wgcnaWith_kmeans_Power", powerEstimate, "finalCombo_CoExpressNetSummaryStats_forAllModules",outputAddOn,".csv")
} else {
 allNetworkModulesCombinedTomSummaryStatisticsFileName = paste0(wgcnaAndKMeansOutputPath, "//", bodyRegion, "_wgcnaOnly_Power", powerEstimate, "finalCombo_CoExpressNetSummaryStats_forAllModules",outputAddOn,".csv")
}



moduleEnrichmentsFilePath = paste0(pathForWGCNA, "//ModuleEnrichments")
dir.create(moduleEnrichmentsFilePath, showWarnings = FALSE)


save(dissTOMPowerOutputNameCSV,
dissTOMPowerOutputNameRData,
tomPowerOutputNameCSV,
tomPowerOutputNameRData,
geneModMembershipOutputNameCSV,
mmpOutputNameCSV,
posCorsOnlyModuleTraitsOutputNameCSV,
modTraitCorrelationOutputNameCSV,
modTraitCorrPValOutputNameCSV,
modTraitCorrAndCorrPValueOutputNameRData,
signifModTraitCorrelationOutputNameCSV,
updatedModuleTraitPhenotypeDFFilePath,
wgcnaWithKmeansAllObjectsExceptTOM,
wgcnaWithKmeansAllObjectsIncludingTOM,
wgcnaWithKmeansAllObjectsExceptTOMWithDate,
wgcnaWithKmeansAllObjectsIncludingTOMWithDate,
wgcnaWithKmeansGeneTraitOrganizedPath,
wgcnaWithKmeansGeneTraitSignificancePath,
wgcnaWithKmeansGeneTraitSigPValuePath,
wgcnaWithKmeansGeneModuleMembershipPath,
wgcnaWithKmeansGeneModMemberPValuePath,
wgcnaWithKmeansGeneTraitsEdgeTablePathCSV,
wgcnaWithKmeansGeneTraitsEdgeTablePathRData,
signifGeneTraitCorrelationOutputNameCSV,
signifPositiveGeneTraitCorrelationOutputNameCSV,
updatedGeneTraitPhenotypeDFFilePath,
wgcnaPowerOutputNameRData,
moduleCountsOutputName,
wgcnaOutputName,
MEsWithKmeansOutputName,
geneAssignmentsWithKmeansOutputNameCSV,
geneAssignmentsWithKmeansOutputNameRData,
keyGeneAssignmentsWithKmeansOutputNameCSV,
keyGeneAssignmentsWithKmeansOutputNameRData,
moduleCountsWithKmeansOutputNameCSV,
hubsCSV,
allNetworkModulesCombinedTomSummaryStatisticsFileName,
powerEstimate,
geneAssignmentsModAssignmentsPhenosWithKmeansOutputNameCSV,
file =  powerRelatedFilePathsRData)

print(paste0(":) please note that this stores the file path information related to the WGCNA power ", powerEstimate, " that will be used: ", powerRelatedFilePathsRData))

}






# dictionary or hashtable with gene symbol as key and entrez id as value :)
geneIDMappingFunction <- function(x){
  entrezIDMappingVec[x]
}

geneIDMappingHashTableFunction <- function(x){
  if (is.null(entrezIDMappingHashTable[[x]])){# == NULL
    return(x)
  } else {
    return(entrezIDMappingHashTable[[x]])
  }
}





pleaseCreateCombined_GeneRegNetwork_From_RTN_Genie3_TrenaEnsemble_and_Trrust2_Only <- function(trrust2DF,
                                                                                          RTN_rDataPart,
                                                                                          genie3RDataPartFinal,
                                                                                          trenaOutpathALLRData_Final,
                                                                                          grnResultsCombinedFileNameCSV,
                                                                                          grnResultsCombinedFileNameRData){

  # please note that this function by Saniya takes in the R Data files from RTN, genie3, and trena (along with the dataframe from TRRUST2.0)
  # it then combines them.
  # please note that each of the files should correspond to a dataframe that has these 5 same headings:
  # TF RegulatedGene Source                      Info    CombinedName
  # in general, as long as the file headings are the same for the 4 sources, this function can combine them and stitch them together
  # however, for subsequent functions, it may be important to ensure that the order is retained properly
  print(paste0(":) Please note that we load in the RData Files for the 3 Gene Regulatory Networks (GRN) from: "))
  load(RTN_rDataPart)
  print(paste0(":) RTN File Path (RData): ", RTN_rDataPart))

  print(paste0(":) RTN:  # of Transcription Factor - Target Gene Relationships in Final RTN GRN: ", nrow(rtnRegulonsDF)))


  load(genie3RDataPartFinal)
  print(paste0(":) Genie3 File Path (RData): ", genie3RDataPartFinal))
  print(paste0(":) Genie3:  # of Transcription Factor - Target Gene Relationships in Final Genie GRN: ", nrow(genie3df)))


  load(trenaOutpathALLRData_Final)
  print(paste0(":) Trena Ensemble Solver File Path (RData): ", trenaOutpathALLRData_Final))

  print(paste0(":) Trena Ensemble Solver:  # of Transcription Factor - Target Gene Relationships in Final Trena GRN: ", nrow(trenaEnsembleDF)))

  print(paste0(":) Please note that we also loaded in TRRUST2.0 gold standard regulatory information (found online) and found these # of Transcription Factor - Target Gene Relationships in the Final Trrust2.0 GRN: ",
               nrow(trrust2DF)))

  resultDF = rtnRegulonsDF
  resultDF = rbind(resultDF, genie3df)
  resultDF = rbind(resultDF, trenaEnsembleDF)
  resultDF = rbind(resultDF, trrust2DF)
  resultDF = unique(resultDF)

  print(paste0(":) Please note that altogether, we have these many rows in our final combined GRN from all 4 sources: ", nrow(resultDF)))
  print(paste0(":) Please note a summary breakdown of the # of unique Transcription Factor - Target Gene Edges from each source (please note some TF-TG relationships may fall in multiple sources (1 to 4 of them)"))
  print(table(resultDF$Source))
  print(":) Please note the first 5 rows:")
  print(head(resultDF))
  dim(resultDF)
  write.csv(resultDF, grnResultsCombinedFileNameCSV)
  save(resultDF, file = grnResultsCombinedFileNameRData)
  print(":) Please note that the results of this combined Gene Regulatory Network are stored here")
  print(paste0(":) in CSV form, please note: ", grnResultsCombinedFileNameCSV))
  print(paste0(":) in RData form, please note: ", grnResultsCombinedFileNameRData))

  return(resultDF)
}







pleaseCreateCombined_GeneRegNetwork_From_RTN_Genie3_TrenaEnsembleSolver_and_Trrust2 <- function(trrust2DF,
                                                                                          RTN_rDataPart,
                                                                                          genie3RDataPartFinal,
                                                                                          trenaOutpathALLRData_Final,
                                                                                          grnResultsCombinedFileNameCSV,
                                                                                          grnResultsCombinedFileNameRData,
																						  pathForPythonCode,
																						  parentName, disease, bodyRegion,
																						  dataScalingOutputMini, geneToEntrezIDMappingPath, minNumSourcesGeneGRN, progressBarPythonIterations){

  # please note that this function by Saniya takes in the R Data files from RTN, genie3, and trena (along with the dataframe from TRRUST2.0)
  # it then combines them.
  # please note that each of the files should correspond to a dataframe that has these 5 same headings:
  # TF RegulatedGene Source                      Info    CombinedName
  # in general, as long as the file headings are the same for the 4 sources, this function can combine them and stitch them together
  # however, for subsequent functions, it may be important to ensure that the order is retained properly
  print(paste0(":) minNumSourcesGeneGRN = ", minNumSourcesGeneGRN))

  print(paste0(":) Please note that we load in the RData Files for the 3 Gene Regulatory Networks (GRN) from: "))
  load(RTN_rDataPart)
  print(paste0(":) RTN File Path (RData): ", RTN_rDataPart))

  print(paste0(":) RTN:  # of Transcription Factor - Target Gene Relationships in Final RTN GRN: ", nrow(rtnRegulonsDF)))


  load(genie3RDataPartFinal)
  print(paste0(":) Genie3 File Path (RData): ", genie3RDataPartFinal))
  print(paste0(":) Genie3:  # of Transcription Factor - Target Gene Relationships in Final Genie GRN: ", nrow(genie3df)))


  load(trenaOutpathALLRData_Final)
  print(paste0(":) Trena Ensemble Solver File Path (RData): ", trenaOutpathALLRData_Final))

  print(paste0(":) Trena Ensemble Solver:  # of Transcription Factor - Target Gene Relationships in Final Trena GRN: ", nrow(trenaEnsembleDF)))

  print(paste0(":) Please note that we also loaded in TRRUST2.0 gold standard regulatory information (found online) and found these # of Transcription Factor - Target Gene Relationships in the Final Trrust2.0 GRN: ",
               nrow(trrust2DF)))

  resultDF = rtnRegulonsDF
  resultDF = rbind(resultDF, genie3df)
  resultDF = rbind(resultDF, trenaEnsembleDF)
  resultDF = rbind(resultDF, trrust2DF)
  resultDF = unique(resultDF)

  print(paste0(":) Please note that altogether, we have these many rows in our final combined GRN from all 4 sources: ", nrow(resultDF)))
  print(paste0(":) Please note a summary breakdown of the # of unique Transcription Factor - Target Gene Edges from each source (please note some TF-TG relationships may fall in multiple sources (1 to 4 of them)"))
  print(table(resultDF$Source))
  print(":) Please note the first 5 rows:")
  print(head(resultDF))
  dim(resultDF)
  write.csv(resultDF, grnResultsCombinedFileNameCSV)
  save(resultDF, file = grnResultsCombinedFileNameRData)
  print(":) Please note that the results of this combined Gene Regulatory Network are stored here")
  print(paste0(":) in CSV form, please note: ", grnResultsCombinedFileNameCSV))
  print(paste0(":) in RData form, please note: ", grnResultsCombinedFileNameRData))


	print(":) Next, please note we run this Python code to organize our GRN and put duplicate TF-TG edges in a more organized way:")
	source_python(pathForPythonCode)

	parentName = filePathOfAll4GRNs # "D:\\organizedAlzheimers\\TranscriptionalRegulator\\outputs\\"

	write.csv(resultDF, grnResultsComboForPythonInputFileNameCSV)

	 masterDF = organizingOutputFileForTFToGeneRelationship(grnResultsComboForPythonInputFileNameCSV,
	 													   parentName, disease,
	 													   bodyRegion, dataScalingOutputMini, geneToEntrezIDMappingPath, minNumSourcesGeneGRN, progressBarPythonIterations)

  return(masterDF)
}


# https://stackoverflow.com/questions/64519640/error-in-summary-connectionconnection-invalid-connection
unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}

buildGRNUsingGenie3 <- function(allsamples, tfsToUse, weightThresholdGenie3, genie3_rDataPart, genie3RDataPartFinal,
                                numOfCoresToUseForGenie3){
  print(paste0(":) please note that this function builds a gene regulatory network using genie3"))
  exprMatr = as.matrix(allsamples)
  row.names(exprMatr) = row.names(allsamples)
  colnames(exprMatr) = colnames(allsamples)
  
  unregister_dopar()
  weightMat <- GENIE3(exprMatr, regulators=tfsToUse, nCores=numOfCoresToUseForGenie3)
  linkList <- getLinkList(weightMat)

  nonZeroWeights = linkList[which(linkList$weight > 0),]
  save(weightMat, linkList, nonZeroWeights, file = genie3_rDataPart)
  write.csv(linkList, genie3CsvPart)
  quantile(nonZeroWeights$weight)
  weightsToUse = nonZeroWeights[which(nonZeroWeights$weight > weightThresholdGenie3),]
  print(paste0(":) Please note that we will use this weight threshold for Genie3 of ",
               weightThresholdGenie3, " and now will use these top: ", nrow(weightsToUse), " weighted edges for our final GRN :)")) # 146993 --> 692,987 # 5872078
  rankOrder = order(weightsToUse$weight, decreasing = TRUE)
  infoVector = paste0("GENIE3 (Rank: ", rankOrder, "; Weight: ", round(weightsToUse$weight, 4), ")")

  genie3df = weightsToUse[,c(1, 2)]
  genie3df$Source = rep("Genie3", nrow(genie3df))
  genie3df$Info = infoVector
  comboVec = paste(weightsToUse[,1], "||", weightsToUse[,2])
  genie3df$CombinedName = comboVec
  colnames(genie3df) = colnamesVec
  print(head(genie3df))
  write.csv(genie3df, genie3CsvPartFinal)
  save(genie3df, weightThreshold, folderName, file = genie3RDataPartFinal)

  print(paste0(":) Please note we have finished with Genie3 :)"))

  return(genie3df)
}



buildGRNUsingTrenaEnsembleSolver <- function(allsamples, tfsToUse, powerEstimate, tfsInGeneExpressionDataRData,
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
                                             numIterationsToPrintProgressForTrenaOutput){
  print(paste0(":) please note that this function builds a gene regulatory network using Trena Ensemble Solver"))
 load(tfsInGeneExpressionDataRData)
  allsamples <- read.csv(filePathOfGeneExpressionDataUsedByWGCNA,
                         header = TRUE, row.names = 1)

  dim(allsamples)
  geneExpressionSamplesVec = row.names(allsamples)
  exprMatr = as.matrix(allsamples)

  allsamplesMatrix = as.matrix(allsamples)
  row.names(allsamplesMatrix) = row.names(allsamples)
  colnames(allsamplesMatrix) = colnames(allsamples)
  zeroSamples = which(rowSums(allsamples) == 0) # please note these are indexes of genes with expression of 0
  if (length(zeroSamples) > 0){
    nonZeroSamples = allsamplesMatrix[-zeroSamples,]
  } else {
    nonZeroSamples = allsamplesMatrix
  }
  print(paste0(":) Please note that we remove these: ", length(zeroSamples),
               " genes that have a total expression of 0 for all samples!"))
  print(dim(allsamplesMatrix))
  print("to: --> ")
  print(dim(nonZeroSamples))
  allsamplesMatrix = nonZeroSamples
  print(dim(allsamplesMatrix))

  genesVec = row.names(allsamplesMatrix)

  print(paste0(":) Please note that there are: ", length(genesVec), " genes to Analyze"))
  print(paste0(":) ", length(tfsToUse), " of those genes are candidate regulators (Transcription Factors)"))

  trenaEnsembleDF = data.frame()
  errorVec = c()
  trenaOutputOrganizedList <- list()
  counter = 1
 # load(tfsInGeneExpressionDataRData) # loading info on TFsToUse
  if (speedUpTrenaEnsembleByEliminatingLassoRelatedModels){
    goalColNamesTrena = c("targetGene", "TF",
                          "pearsonCoeff",           "rfScore",
                          "betaRidge",              "spearmanCoeff",
                          "xgboost",                "region",
                          "tissue",                 "diseaseVec",
                          "power",                  "geneExpressionDataUsed",
                          "regNetworkDataSource",   "regNetworkMethodology")
    regNetworkMethodology = rep("Combo of: pearsonCoeff, rfScore, betaRidge, spearmanCoeff, xgboost")
    # please note that we remove lasso and lassopv from the list of solvers we use here, so we can speed up trena
    solverNamesVec = c("pearson", "randomForest", "ridge", "spearman",
                       "xgboost")
  } else {
    goalColNamesTrena = c("targetGene", "TF",
                          "betaLasso",              "lassoPValue",
                          "pearsonCoeff",           "rfScore",
                          "betaRidge",              "spearmanCoeff",
                          "xgboost",                "region",
                          "tissue",                 "diseaseVec",
                          "power",                  "geneExpressionDataUsed",
                          "regNetworkDataSource",   "regNetworkMethodology")
    regNetworkMethodology = rep("Combo of: betaLasso, lassoPValue, pearsonCoeff, rfScore, betaRidge, spearmanCoeff, xgboost")
    solverNamesVec = c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman",
                       "xgboost")
  }


  if (useVarianceFilterToReduceCandidateTFsForTrenaToSpeedUpProcess){
    if (is.null(varianceSizeForTFsAndTargetGeneForTrena)){
      varianceSizeFilterToUseForTrena = 0.25 # please note this is a default variance filter used just in case user forgot to provide one
      print(paste(":( Please note that since ye put TRUE for useVarianceFilterToReduceCandidateTFsForTrenaToSpeedUpProcess, but forgot to provide a non-null value for varianceSizeForTFsAndTargetGeneForTrena (detween 0 and 1)",
                  "we will be using our default value of: ", varianceSizeForTFsAndTargetGeneForTrena ))
      print(paste0("That is, based on Trena Documentation, this variance filter default we set of:", varianceSizeFilterToUseForTrena, " means that we will find all ",
                   "Transcription Factors (TFs) with variance that is within ", varianceSizeFilterToUseForTrena*100, "% of a given target gene's respective variance.  Thus, the list of candidate TFs will be based on each target gene.  This method will save some computational time :), hopefully!"))

    } else {
      varianceSizeFilterToUseForTrena = varianceSizeForTFsAndTargetGeneForTrena
      print(paste(":( Please note that we will be using your value of varianceSizeForTFsAndTargetGeneForTrena = ",
                  varianceSizeForTFsAndTargetGeneForTrena, " for the varianceSizeFilter for Trena."))
      print(paste0("That is, based on Trena Documentation, this variance filter ye chose of:", varianceSizeFilterToUseForTrena, " means that we will find all ",
                   "Transcription Factors (TFs) with variance that is within ", varianceSizeFilterToUseForTrena*100, "% of a given target gene's respective variance.  Thus, the list of candidate TFs will be based on each target gene.  This method will save some computational time :), hopefully!"))

    }
  }

  for (i in 1:length(genesVec)){
    tryCatch({
      if (i %% numIterationsAfterForSavingTrenaOutput == 0){
        outpath = paste0(filePathOfTrena, "//trena_", bodyRegion, "_endRegion", i, "_",
                         disease, "_", tissueName, ".csv")
        trenaEnsembleDF = do.call("rbind", trenaOutputOrganizedList)
        colnames(trenaEnsembleDF) = goalColNamesTrena
        print(paste0(":) please note we wrote out the file here for now:", outpath))
        write.csv(trenaEnsembleDF, outpath)
      }

      if (i %% numIterationsToPrintProgressForTrenaOutput == 0){
        print(paste0(":) please note we are on gene ", i, " of ", length(genesVec), ": ", targetGeneName))
      }
      targetGeneName = genesVec[i]

      if (useVarianceFilterToReduceCandidateTFsForTrenaToSpeedUpProcess){
        variance.filter = VarianceFilter(mtx.assay = allsamplesMatrix,
                                         targetGene = targetGeneName,
                                         varSize = varianceSizeFilterToUseForTrena) #0.25)
        tf.list <- getCandidates(variance.filter)
        length(tf.list$tf.vars)
        candidateTFs = tf.list$tfs
      } else {
        candidateTFs = tfsToUse # all TFs we have in our list
      }

      # solvers we will use
      ensemble.solver <- suppressWarnings(EnsembleSolver(mtx.assay = allsamplesMatrix,
                                                         targetGene = targetGeneName,
                                                         solverNames = solverNamesVec,
                                                         candidateRegulators = candidateTFs))


      #geneCutoff = 0.05)
      tbl.out <- suppressWarnings(run(ensemble.solver))
      targetGeneName = rep(targetGeneName, nrow(tbl.out))
      regNetworkDataSource = rep("TReNA Ensemble Solver", nrow(tbl.out))
      #regNetworkMethodology = rep("Combo of: betaLasso, lassoPValue, pearsonCoeff, rfScore, betaRidge, spearmanCoeff, xgboost")
      power = rep(powerEstimate, nrow(tbl.out))
      region = rep(bodyRegion, nrow(tbl.out))
      tissue = rep(tissueName, nrow(tbl.out))
      diseaseVec = rep(disease, nrow(tbl.out))
      geneExpressionDataUsed = rep(dataScaling, nrow(tbl.out))
      tbl.outDF = data.frame(cbind(targetGeneName, tbl.out, region, tissue, diseaseVec, power,
                                   geneExpressionDataUsed, regNetworkDataSource, regNetworkMethodology))
      colnames(tbl.outDF)[1] = "targetGene"
      colnames(tbl.outDF)[2] = "TF"
      trenaOutputOrganizedList[[i]] = tbl.outDF
    }, error = function(e) {
      print(paste0(":( error for index: ", i))
      #errorVec = c(errorVec, i)
    })
  }



  trenaEnsembleDF = do.call("rbind", trenaOutputOrganizedList)
  colnames(trenaEnsembleDF) = goalColNamesTrena
  # perhaps a filtering column as well, where here it is on default: geneValueCutoff = 0.10 (so only top 10% most important are shown)
  colnames(trenaEnsembleDF)[1] = "targetGene"
  colnames(trenaEnsembleDF)[2] = "TF"
  head(trenaEnsembleDF)

  write.csv(trenaEnsembleDF, trenaOutpathALLCSV) #outpathALL) #"D://trenaEnsembleDLPFC_NEWEREST.csv")
  save(trenaEnsembleDF, file = trenaOutpathALLRData) #"D://trenaEnsembleDLPFC_NEWEREST.RData")

  #filePathOfTrena: stitching together the csvs :)
  # trenaCSVs = list.files(filePathOfTrena, "*.csv")
  # trenaEnsembleDF = data.frame()
  # for (i in 1:length(trenaCSVs)){
  #   filePathName = paste0(filePathOfTrena, "//", trenaCSVs[i])
  #   trenaDF = read.csv(filePathName, header = TRUE)
  #   if (i == 1){
  #     trenaEnsembleDF = trenaDF
  #   } else {
  #   trenaEnsembleDF = rbind(trenaEnsembleDF, trenaDF)
  #   }
  # }
  dim(trenaEnsembleDF)
  trenaEnsembleDF = unique(trenaEnsembleDF)
  dim(trenaEnsembleDF)
  head(trenaEnsembleDF)

  if (speedUpTrenaEnsembleByEliminatingLassoRelatedModels){
    trenaEnsembleDF$Info = paste0("TReNA Ensemble Solver: (pearsonCoeff: ",round(trenaEnsembleDF$pearsonCoeff, numberOfRoundingDigits),
                                  "; rfScore: ",round(trenaEnsembleDF$rfScore, numberOfRoundingDigits),
                                  "; betaRidge: ",round(trenaEnsembleDF$betaRidge, numberOfRoundingDigits),
                                  "; spearmanCoeff: ",round(trenaEnsembleDF$spearmanCoeff, numberOfRoundingDigits),
                                  "; xgboost: ",round(trenaEnsembleDF$xgboost, numberOfRoundingDigits), ")"
    )
  } else { # please note that we include everything here
    trenaEnsembleDF$Info = paste0("TReNA Ensemble Solver: (betaLasso: ",
                                  round(trenaEnsembleDF$betaLasso, numberOfRoundingDigits),
                                  "; lassoPValue: ",round(trenaEnsembleDF$lassoPValue, numberOfRoundingDigits),
                                  "; pearsonCoeff: ",round(trenaEnsembleDF$pearsonCoeff, numberOfRoundingDigits),
                                  "; rfScore: ",round(trenaEnsembleDF$rfScore, numberOfRoundingDigits),
                                  "; betaRidge: ",round(trenaEnsembleDF$betaRidge, numberOfRoundingDigits),
                                  "; spearmanCoeff: ",round(trenaEnsembleDF$spearmanCoeff, numberOfRoundingDigits),
                                  "; xgboost: ",round(trenaEnsembleDF$xgboost, numberOfRoundingDigits), ")"
    )
  }
  head(trenaEnsembleDF)
  trenaEnsembleDF$Source = rep("TReNa", nrow(trenaEnsembleDF))
  colnames(trenaEnsembleDF)
  trenaCols = c(grep("\\bTF\\b", colnames(trenaEnsembleDF)),
                grep("Gene\\b", colnames(trenaEnsembleDF)),
                grep("\\bSource\\b", colnames(trenaEnsembleDF)),
                grep("\\bInfo\\b", colnames(trenaEnsembleDF)))

  trenaEnsembleDF = trenaEnsembleDF[,trenaCols]
  head(trenaEnsembleDF)
  colnames(trenaEnsembleDF) = c("TF", "RegulatedGene", "Source",    "Info")
  head(trenaEnsembleDF)
  trenaEnsembleDF$CombinedName = paste0(trenaEnsembleDF$TF, " || ", trenaEnsembleDF$RegulatedGene)
  #head(resultDF)
  head(trenaEnsembleDF)


  write.csv(trenaEnsembleDF, trenaOutpathALLCSV_Final) #outpathALL) #"D://trenaEnsembleDLPFC_NEWEREST.csv")
  save(trenaEnsembleDF, file = trenaOutpathALLRData_Final) #"D://trenaEnsembleDLPFC_NEWEREST.RData")

  print(dim(trenaEnsembleDF))
  print(head(trenaEnsembleDF))
  return(trenaEnsembleDF)
  }


findingModulePhenotypeCorrelationsAndPValues <- function(MEs,
                                                         datTraits, modTraitCorrelationOutputNameCSV,
                                                         datExpr, modTraitCorrPValOutputNameCSV,
                                                         modTraitCorrAndCorrPValueOutputNameRData,
                                                         modulePhenoAssociationFilePath,
                                                         performAdditionalKMeansStep,
                                                         pValueCutOffForSignificance,powerEstimate,
                                                         signifModTraitCorrelationOutputNameCSV,positiveCorrelationsOnly,
                                                         posCorsOnlyModuleTraitsOutputNameCSV, numOfRoundingDigits, dynamicColorsKmeans,
                                                         fullDiseaseName, tissueName, bodyRegion, dataScaling,
                                                         updatedModuleTraitPhenotypeDFFilePath){
  print(paste0(":) Please note that this function helps find Module-Trait Correlations that are Significant as well :)"))
  moduleTraitCor = cor(MEs, datTraits, use = "p") #, robustY = FALSE, maxPOutliers = 0.1);
  
  write.csv(data.frame(moduleTraitCor), modTraitCorrelationOutputNameCSV)
  
  nSamples = nrow(datExpr)
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
  write.csv(data.frame(moduleTraitPvalue), modTraitCorrPValOutputNameCSV)
  currentDate = as.character(now()) # current date and time str_replace_all(Sys.Date(), "-", "_")
  
  save(moduleTraitCor, moduleTraitPvalue, currentDate, file = modTraitCorrAndCorrPValueOutputNameRData)
  
  
  
  setwd(modulePhenoAssociationFilePath)
  if (performAdditionalKMeansStep){
    info = "WGCNA with kMeans Applied"
  } else {
    info = "WGCNA Only"
  }
  
  modTraitDF1 = melt(moduleTraitCor)
  modTraitDF1$newName = paste0(modTraitDF1[,1], "_", modTraitDF1[,2])
  
  colnames(modTraitDF1) = c("module", "phenotype", "correlation_r", "comboName")
  
  modTraitDF2 = melt(moduleTraitPvalue)
  modTraitDF2$newName = paste0(modTraitDF2[,1], "_", modTraitDF2[,2])
  modTraitDF2 = modTraitDF2[,c(4, 3)]
  colnames(modTraitDF2) = c("comboName", "corr_pValue")
  
  combinedModTraitDF = merge(modTraitDF1, modTraitDF2, by.x = "comboName", by.y = "comboName", sort = TRUE)
  head(combinedModTraitDF)
  rm(modTraitDF1)
  rm(modTraitDF2)
  dim(combinedModTraitDF) # [1] 594   5
  significantModTraits = combinedModTraitDF[which(combinedModTraitDF$corr_pValue < pValueCutOffForSignificance),]
  dim(significantModTraits) # [1] 86  5
  significantModTraits$powerUsed  = rep(powerEstimate, nrow(significantModTraits))
  significantModTraits$dateCreated = rep(currentDate, nrow(significantModTraits))
  
  write.csv(significantModTraits, signifModTraitCorrelationOutputNameCSV) #paste0(outputPath, "dlpfc_wgcna_with_kmeans_Power", powerEstimate, "_StatisticallySignificantModuleTraitCorrelations_.csv"))
  
  # please note that if we only want the positive correlations r (Between genes and modules), this will help us limit
  # down our dataframe to only the positive ones:
  if (positiveCorrelationsOnly == TRUE){
    significantModTraits = significantModTraits[which(significantModTraits$correlation_r > 0),]
    #write.csv(significantModTraits, paste0(outputPath, "wgcna_with_kmeans_Power", powerEstimate, "_StatisticallySignificantModuleTraitCorrelations_POSITIVEonly_.csv"))
    write.csv(significantModTraits, posCorsOnlyModuleTraitsOutputNameCSV)
  }
  numRoundingDigits = numOfRoundingDigits
  significantModTraits$roundedCorr = round(significantModTraits$correlation_r, numRoundingDigits)
  head(significantModTraits)
  uniqueModules = as.vector(unique(significantModTraits[,2])) #
  
  moduleNamesVec = unique(dynamicColorsKmeans)
  print(paste0(":) Please note that there are: ", length(uniqueModules), " total gene modules that were detected using WGCNA with power estimate: ", powerEstimate, " that have significant Module-Trait Associations"))
  
  simplePositivePhenotypesVec = c() # please note that this holds information on the positive phenotypes (positively-correlated traits)
  moreInfoPositivePhenotypesVec = c() # please note that this also includes the correlation (r) values
  numOfPositivePhenotypesVec = c() # please note that this holds the number of positive phenotypes for that module
  
  if (positiveCorrelationsOnly == FALSE){
    simpleNegativePhenotypesVec = c() # please note that this holds information on the negative phenotypes (negatively-correlated traits)
    moreInfoNegativePhenotypesVec = c() # please note that this also includes the correlation (r) values
    numOfNegativePhenotypesVec = c() # please note that this holds the number of negative phenotypes for that module
    totalNumPhenotypesVec = c() # please note that this holds the total number of phenotypes (pos and neg)
  }
  
  for (i in 1:length(uniqueModules)){
    currentMod = uniqueModules[i]
    
    print(paste0(":) please note i = ", i, " and current module is: ", currentMod))
    miniDF = significantModTraits[which(significantModTraits$module == currentMod),]
    posDF = miniDF[which(miniDF$correlation_r > 0),] # should be the same as miniDF, but just to be safe
    simplePosPheno = "" # please initialize as empty strings
    moreInfoPosPheno = "" # please initialize as empty strings
    numPos = nrow(posDF)
    if (numPos > 0){ # please note that we have some positive rows
      phenotypes = as.vector(posDF$phenotype)
      rVals = as.vector(posDF$roundedCorr)
      
      for (j in 1:nrow(posDF)){
        pheno = phenotypes[j]
        r = rVals[j]
        moreInfoVal = paste0(pheno, " (r = ", r, ")") # phenotype and r value :)
        if (j == 1){
          simplePosPheno = paste0(simplePosPheno,  pheno)
          moreInfoPosPheno = paste0(moreInfoPosPheno,  moreInfoVal)
        } else { # we need to add a comma to the previous
          simplePosPheno = paste0(simplePosPheno, ", ", pheno)
          moreInfoPosPheno = paste0(moreInfoPosPheno, ", ", moreInfoVal)
        }
        
      }
    }
    # please update and append info from this module :)
    simplePositivePhenotypesVec = c(simplePositivePhenotypesVec, simplePosPheno) # please note that this holds information on the positive phenotypes (positively-correlated traits)
    moreInfoPositivePhenotypesVec = c(moreInfoPositivePhenotypesVec, moreInfoPosPheno) # please note that this also includes the correlation (r) values
    numOfPositivePhenotypesVec = c(numOfPositivePhenotypesVec, numPos)
    # if we are also analyzing negative correlations :)
    if (positiveCorrelationsOnly == FALSE){
      simpleNegPheno = "" # please initialize as empty strings
      moreInfoNegPheno = "" # please initialize as empty strings
      negDF = miniDF[which(miniDF$correlation_r < 0),]
      numNeg = nrow(negDF)
      
      if (numNeg > 0){ # please note that we have some negative rows
        phenotypes = as.vector(negDF$phenotype)
        rVals = as.vector(negDF$roundedCorr)
        
        for (j in 1:nrow(negDF)){
          pheno = phenotypes[j]
          r = rVals[j]
          moreInfoVal = paste0(pheno, " (r = ", r, ")") # phenotype and r value :)
          if (j == 1){
            simpleNegPheno = paste0(simpleNegPheno,  pheno)
            moreInfoNegPheno = paste0(moreInfoNegPheno,  moreInfoVal)
          } else { # we need to add a comma to the previous
            simpleNegPheno = paste0(simpleNegPheno, ", ", pheno)
            moreInfoNegPheno = paste0(moreInfoNegPheno, ", ", moreInfoVal)
          }
          
        }
      }
      # please update and append info from this module :)
      simpleNegativePhenotypesVec = c(simpleNegativePhenotypesVec, simpleNegPheno) # please note that this holds information on the positive phenotypes (positively-correlated traits)
      moreInfoNegativePhenotypesVec = c(moreInfoNegativePhenotypesVec, moreInfoNegPheno) # please note that this also includes the correlation (r) values
      numOfNegativePhenotypesVec = c(numOfNegativePhenotypesVec, numNeg)
      totalNumPhenotypesVec = c(totalNumPhenotypesVec, (numPos + numNeg))
    }
    
  }
  if (positiveCorrelationsOnly == TRUE){
    updatedModuleTraitPhenotypeDF = data.frame(uniqueModules,
                                               simplePositivePhenotypesVec,
                                               moreInfoPositivePhenotypesVec,
                                               numOfPositivePhenotypesVec)
  } else {
    updatedModuleTraitPhenotypeDF = data.frame(uniqueModules,
                                               simplePositivePhenotypesVec,
                                               moreInfoPositivePhenotypesVec,
                                               numOfPositivePhenotypesVec,
                                               simpleNegativePhenotypesVec,
                                               moreInfoNegativePhenotypesVec,
                                               numOfNegativePhenotypesVec,
                                               totalNumPhenotypesVec)
  }
  updatedModuleTraitPhenotypeDF
  
  dateCreated = as.character(now())
  updatedModuleTraitPhenotypeDF$moduleName = str_replace_all(updatedModuleTraitPhenotypeDF$uniqueModules,
                                                             "ME", "")
  colnames(updatedModuleTraitPhenotypeDF)[grep("uniqueModules", colnames(updatedModuleTraitPhenotypeDF))] = "moduleEigenegene"
  updatedModuleTraitPhenotypeDF$powerUsed = rep(powerEstimate, nrow(updatedModuleTraitPhenotypeDF))
  updatedModuleTraitPhenotypeDF$disease = rep(fullDiseaseName, nrow(updatedModuleTraitPhenotypeDF))
  updatedModuleTraitPhenotypeDF$tissueName = rep(tissueName, nrow(updatedModuleTraitPhenotypeDF))
  updatedModuleTraitPhenotypeDF$bodyRegion = rep(bodyRegion, nrow(updatedModuleTraitPhenotypeDF))
  updatedModuleTraitPhenotypeDF$dataScaling = rep(dataScaling, nrow(updatedModuleTraitPhenotypeDF))
  updatedModuleTraitPhenotypeDF$info = rep(info, nrow(updatedModuleTraitPhenotypeDF))
  updatedModuleTraitPhenotypeDF$dateCreated = rep(dateCreated, nrow(updatedModuleTraitPhenotypeDF))
  
  print(head(updatedModuleTraitPhenotypeDF))
  write.csv(updatedModuleTraitPhenotypeDF, updatedModuleTraitPhenotypeDFFilePath)
  
  outputsToReturn = list()
  outputsToReturn$significantModTraits = significantModTraits
  outputsToReturn$updatedModuleTraitPhenotypeDF = updatedModuleTraitPhenotypeDF
  outputsToReturn$moduleTraitCor = moduleTraitCor
  outputsToReturn$moduleTraitPvalue = moduleTraitPvalue
  outputsToReturn$info = info
  return(outputsToReturn)
}



findingGenePhenotypeCorrelationsAndPValues <- function(datTraits, datExpr,
                                                       wgcnaWithKmeansGeneTraitSignificancePath,
                                                       wgcnaWithKmeansGeneTraitSigPValuePath,
                                                       bodyRegion, tissueName, fullDiseaseName, dataScaling,
                                                       geneAndEntrezIDMappingDF,
                                                       wgcnaWithKmeansGeneTraitOrganizedPath,
                                                       updatedGeneTraitPhenotypeDFFilePath,
                                                       wgcnaWithKmeansGeneTraitsEdgeTablePathCSV,
                                                       wgcnaWithKmeansGeneTraitsEdgeTablePathRData,
                                                       powerEstimate,signifGeneTraitCorrelationOutputNameCSV,
                                                       signifPositiveGeneTraitCorrelationOutputNameCSV,
                                                       positiveCorrelationsOnly, numberOfRoundingDigits,
                                                       dynamicColorsKmeans){
  
  
  print(paste0(":) Please note that this function helps find Gene-Trait Correlations that are Significant as well :)"))
  
  geneTraitSignificanceResultsList = pleaseGetTheGeneTraitSignificanceValuesForTheTraits(datTraits,datExpr,
                                                                                         wgcnaWithKmeansGeneTraitSignificancePath,
                                                                                         wgcnaWithKmeansGeneTraitSigPValuePath,
                                                                                         bodyRegion, tissueName, fullDiseaseName, dataScaling,
                                                                                         geneAndEntrezIDMappingDF)
  geneTraitSignificanceDF = data.frame(geneTraitSignificanceResultsList$geneTraitSignificanceDF)
  GSPvalueDF = data.frame(geneTraitSignificanceResultsList$GSPvalueDF)
  
  dim(geneTraitSignificanceDF)
  dim(GSPvalueDF)
  
  head(geneTraitSignificanceDF)
  
  geneTraitSignificanceCols = grep("GeneSignificance_", colnames(geneTraitSignificanceDF))
  
  write.csv(geneTraitSignificanceDF, wgcnaWithKmeansGeneTraitOrganizedPath) #"D://geneTraitSignificanceDFNeuroscience.csv")
  gene_names = row.names(geneTraitSignificanceDF)
  geneTraitSignifCorrelation = cbind(gene_names,
                                     geneTraitSignificanceDF[,geneTraitSignificanceCols])
  row.names(geneTraitSignifCorrelation) = row.names(geneTraitSignificanceDF)
  geneTraitSignifCorrelationMelt = melt(geneTraitSignifCorrelation)
  geneTraitSignifCorrelationMelt[,2] = str_replace_all(geneTraitSignifCorrelationMelt[,2], "GeneSignificance_", "")
  colnames(geneTraitSignifCorrelationMelt) = c("geneName", "phenotype", "geneTraitCorrelation")
  head(geneTraitSignifCorrelationMelt)
  gspSignificanceCols = grep("GeneSignificancePValue", colnames(GSPvalueDF))
  
  geneTraitSignifCorrPVal = cbind(gene_names, GSPvalueDF[,gspSignificanceCols])
  geneTraitSignifCorrPValMelt = melt(geneTraitSignifCorrPVal)
  geneTraitSignifCorrPValMelt[,2] = str_replace_all(geneTraitSignifCorrPValMelt[,2], "GeneSignificancePValue_", "")
  colnames(geneTraitSignifCorrPValMelt) = c("geneName", "phenotype", "geneTraitCorrPVal")
  head(geneTraitSignifCorrPValMelt)
  
  
  combinedGeneTraitDF = merge(geneTraitSignifCorrelationMelt, geneTraitSignifCorrPValMelt, by = c("geneName", "phenotype"))
  dim(combinedGeneTraitDF)
  head(combinedGeneTraitDF)
  
  naCorrs = which(is.na(combinedGeneTraitDF$geneTraitCorrelation))
  if (length(naCorrs) > 0){
    combinedGeneTraitDF = combinedGeneTraitDF[-naCorrs,]
  }
  dim(combinedGeneTraitDF)
  #library("lubridate")
  currentDate = as.character(now())
  write.csv(combinedGeneTraitDF, wgcnaWithKmeansGeneTraitsEdgeTablePathCSV)
  save(combinedGeneTraitDF, currentDate, file = wgcnaWithKmeansGeneTraitsEdgeTablePathRData)
  
  
  
  significantGeneTraitsDF = combinedGeneTraitDF[which(combinedGeneTraitDF$geneTraitCorrPVal < pValueCutOffForSignificance),]
  dim(significantGeneTraitsDF) # [1] 86  5
  significantGeneTraitsDF$powerUsed  = rep(powerEstimate, nrow(significantGeneTraitsDF))
  significantGeneTraitsDF$dateCreated = rep(currentDate, nrow(significantGeneTraitsDF))
  significantGeneTraitsDF$dataScaling  = rep(dataScaling, nrow(significantGeneTraitsDF))
  
  write.csv(significantGeneTraitsDF, signifGeneTraitCorrelationOutputNameCSV) #paste0(outputPath, "dlpfc_wgcna_with_kmeans_pow", powerEstimate, "_StatisticallySignificantModuleTraitCorrelations_.csv"))
  
  # please note that if we only want the positive correlations r (Between genes and modules), this will help us limit
  # down our dataframe to only the positive ones:
  if (positiveCorrelationsOnly == TRUE){
    significantGeneTraitsDF = significantGeneTraitsDF[which(significantGeneTraitsDF$geneTraitCorrelation > 0),]
    write.csv(significantGeneTraitsDF, signifPositiveGeneTraitCorrelationOutputNameCSV)
    
  }
  #write.csv(significantGeneTraitsDF, signifPositiveGeneTraitCorrelationOutputNameCSV)
  
  numRoundingDigits = numberOfRoundingDigits #3
  significantGeneTraitsDF$roundedCorr = round(significantGeneTraitsDF$geneTraitCorrelation, numRoundingDigits)
  
  head(significantGeneTraitsDF)
  uniqueGenes = as.vector(unique(significantGeneTraitsDF[,1])) #
  moduleNamesVec = unique(dynamicColorsKmeans)
  #as.vector(unique(significantGeneTraitsDF[,2]))
  print(paste0(":) Please note that there are: ", length(uniqueGenes), "  genes out of ",
               length(unique(combinedGeneTraitDF[,1])), " total genes that were detected using WGCNA with power estimate: ", powerEstimate, " (and kmeans after) that have significant Gene-Trait Associations"))
  
  simplePositivePhenotypesVec = c() # please note that this holds information on the positive phenotypes (positively-correlated traits)
  moreInfoPositivePhenotypesVec = c() # please note that this also includes the correlation (r) values
  numOfPositivePhenotypesVec = c() # please note that this holds the number of positive phenotypes for that module
  
  if (positiveCorrelationsOnly == FALSE){
    simpleNegativePhenotypesVec = c() # please note that this holds information on the negative phenotypes (negatively-correlated traits)
    moreInfoNegativePhenotypesVec = c() # please note that this also includes the correlation (r) values
    numOfNegativePhenotypesVec = c() # please note that this holds the number of negative phenotypes for that module
    totalNumPhenotypesVec = c() # please note that this holds the total number of phenotypes (pos and neg)
  }
  
  for (i in 1:length(uniqueGenes)){
    currentGene = uniqueGenes[i]
    
    print(paste0(":) please note i = ", i, " and current gene is: ", currentGene))
    miniDF = significantGeneTraitsDF[which(significantGeneTraitsDF$geneName == currentGene),]
    posDF = miniDF[which(miniDF$geneTraitCorrelation > 0),] # should be the same as miniDF, but just to be safe
    simplePosPheno = "" # please initialize as empty strings
    moreInfoPosPheno = "" # please initialize as empty strings
    numPos = nrow(posDF)
    if (numPos > 0){ # please note that we have some positive rows
      phenotypes = as.vector(posDF$phenotype)
      rVals = as.vector(posDF$roundedCorr)
      
      for (j in 1:nrow(posDF)){
        pheno = phenotypes[j]
        r = rVals[j]
        moreInfoVal = paste0(pheno, " (r = ", r, ")") # phenotype and r value :)
        if (j == 1){
          simplePosPheno = paste0(simplePosPheno,  pheno)
          moreInfoPosPheno = paste0(moreInfoPosPheno,  moreInfoVal)
        } else { # we need to add a comma to the previous
          simplePosPheno = paste0(simplePosPheno, ", ", pheno)
          moreInfoPosPheno = paste0(moreInfoPosPheno, ", ", moreInfoVal)
        }
        
      }
    }
    # please update and append info from this module :)
    simplePositivePhenotypesVec = c(simplePositivePhenotypesVec, simplePosPheno) # please note that this holds information on the positive phenotypes (positively-correlated traits)
    moreInfoPositivePhenotypesVec = c(moreInfoPositivePhenotypesVec, moreInfoPosPheno) # please note that this also includes the correlation (r) values
    numOfPositivePhenotypesVec = c(numOfPositivePhenotypesVec, numPos)
    # if we are also analyzing negative correlations :)
    if (positiveCorrelationsOnly == FALSE){
      simpleNegPheno = "" # please initialize as empty strings
      moreInfoNegPheno = "" # please initialize as empty strings
      negDF = miniDF[which(miniDF$geneTraitCorrelation < 0),]
      numNeg = nrow(negDF)
      
      if (numNeg > 0){ # please note that we have some negative rows
        phenotypes = as.vector(negDF$phenotype)
        rVals = as.vector(negDF$roundedCorr)
        
        for (j in 1:nrow(negDF)){
          pheno = phenotypes[j]
          r = rVals[j]
          moreInfoVal = paste0(pheno, " (r = ", r, ")") # phenotype and r value :)
          if (j == 1){
            simpleNegPheno = paste0(simpleNegPheno,  pheno)
            moreInfoNegPheno = paste0(moreInfoNegPheno,  moreInfoVal)
          } else { # we need to add a comma to the previous
            simpleNegPheno = paste0(simpleNegPheno, ", ", pheno)
            moreInfoNegPheno = paste0(moreInfoNegPheno, ", ", moreInfoVal)
          }
          
        }
      }
      # please update and append info from this module :)
      simpleNegativePhenotypesVec = c(simpleNegativePhenotypesVec, simpleNegPheno) # please note that this holds information on the positive phenotypes (positively-correlated traits)
      moreInfoNegativePhenotypesVec = c(moreInfoNegativePhenotypesVec, moreInfoNegPheno) # please note that this also includes the correlation (r) values
      numOfNegativePhenotypesVec = c(numOfNegativePhenotypesVec, numNeg)
      totalNumPhenotypesVec = c(totalNumPhenotypesVec, (numPos + numNeg))
    }
    
  }
  if (positiveCorrelationsOnly == TRUE){
    updatedGeneTraitPhenotypeDF = data.frame(uniqueGenes,
                                             simplePositivePhenotypesVec,
                                             moreInfoPositivePhenotypesVec,
                                             numOfPositivePhenotypesVec)
  } else {
    updatedGeneTraitPhenotypeDF = data.frame(uniqueGenes,
                                             simplePositivePhenotypesVec,
                                             moreInfoPositivePhenotypesVec,
                                             numOfPositivePhenotypesVec,
                                             simpleNegativePhenotypesVec,
                                             moreInfoNegativePhenotypesVec,
                                             numOfNegativePhenotypesVec,
                                             totalNumPhenotypesVec)
  }
  head(updatedGeneTraitPhenotypeDF)
  
  dateCreated = as.character(now())
  
  colnames(updatedGeneTraitPhenotypeDF)[grep("uniqueGenes", colnames(updatedGeneTraitPhenotypeDF))] = "geneName"
  updatedGeneTraitPhenotypeDF$powerUsed = rep(powerEstimate, nrow(updatedGeneTraitPhenotypeDF))
  updatedGeneTraitPhenotypeDF$disease = rep(fullDiseaseName, nrow(updatedGeneTraitPhenotypeDF))
  updatedGeneTraitPhenotypeDF$tissueName = rep(tissueName, nrow(updatedGeneTraitPhenotypeDF))
  updatedGeneTraitPhenotypeDF$bodyRegion = rep(bodyRegion, nrow(updatedGeneTraitPhenotypeDF))
  updatedGeneTraitPhenotypeDF$dataScaling = rep(dataScaling, nrow(updatedGeneTraitPhenotypeDF))
  #updatedGeneTraitPhenotypeDF$info = rep(info, nrow(updatedGeneTraitPhenotypeDF))
  updatedGeneTraitPhenotypeDF$dateCreated = rep(dateCreated, nrow(updatedGeneTraitPhenotypeDF))
  
  significantGeneNames = updatedGeneTraitPhenotypeDF[,1]
  row.names(updatedGeneTraitPhenotypeDF) = significantGeneNames
  updatedGeneTraitPhenotypeDF = updatedGeneTraitPhenotypeDF[,-1]
  head(updatedGeneTraitPhenotypeDF)
  updatedGeneTraitPhenotypeDF = merge(geneAndEntrezIDMappingDF, updatedGeneTraitPhenotypeDF, by = 0)
  rownames(updatedGeneTraitPhenotypeDF) = significantGeneNames
  updatedGeneTraitPhenotypeDF = updatedGeneTraitPhenotypeDF[,-1]
  colnames(updatedGeneTraitPhenotypeDF)[1] = "entrezID"
  head(updatedGeneTraitPhenotypeDF)
  print(head(updatedGeneTraitPhenotypeDF))
  write.csv(updatedGeneTraitPhenotypeDF, updatedGeneTraitPhenotypeDFFilePath)
  
  outputsToReturn <- list()
  outputsToReturn$updatedGeneTraitPhenotypeDF = updatedGeneTraitPhenotypeDF
  outputsToReturn$geneTraitSignificanceDF = geneTraitSignificanceDF
  outputsToReturn$GSPvalueDF = GSPvalueDF
  outputsToReturn$significantGeneTraitsDF = significantGeneTraitsDF
  return(outputsToReturn)
}





performingKMeansAfterWGCNA <- function(dynamicColorsPower,
                                       performAdditionalKMeansStep,
                                       datExpr, meg, nIterations){
  ############################################################
  
  # Please calculate initial Module Eigengenes
  # Please note that this part is from: https://github.com/juanbot/km2gcn
  #ALGORITHM SPECIFICATION
  #Step 1. Let D be the expression data in which dij in D represents the expression value for
  #sample i and gene j, being s samples and g genes in total.
  #Step 2. Construct the partition by the WGCNA process, let P_D={m_1, m_2, ..., m_n} be
  #that partition where m_k is the k-th module.
  #Step 3. Get the eigengenes for each module within the partition, E={e_1, e_2, ..., e_n}
  #Step 4. Set up the k-means clustering
  #Step 4.1. Set k to n
  #Step 4.2. Set the centroids C to the eigengenes E, thus C to E
  #Step 5. Run the algorithm and monitor its evolution
  #Step 5.1 Set iterations to 0
  #Step 5.2 Create a new partition P', given C with n modules such that, for each gene, 1 <=
  #		j <= g, g_j belongs to the module c_t in C such that a distance meassure d(g_j,c_t) is
  #		minimum.
  #		Step 5.3 Calculate eigengenes of P', giving a new E'
  #		Step 5.4 Evaluate the progress. If progress done, set iterations to iterations + 1 and
  #		C to E' and go to step 5.2
  #Step 5.5 Finish
  #
  #####
  
  moduleColors = dynamicColorsPower
  mods = unique(moduleColors)
  
  #Please Get rid of the MEs not found in mods
  mods = paste0("ME",mods)
  
  
  ####
  #Gather the current partition we start from
  partition.in.colors <- moduleColors
  if(performAdditionalKMeansStep){
    print("Please note that we are getting initial eigengenes for the K-Means Step after WGCNA.  Code is adapted from: https://github.com/juanbot/km2gcn  :)")
    if(sum(partition.in.colors == "grey") < min.genes.for.grey){
      eigengenes = moduleEigengenes(datExpr,partition.in.colors, excludeGrey=TRUE)
    } else {
      eigengenes = moduleEigengenes(datExpr,partition.in.colors, excludeGrey=F)
    }
    cat("We got",length(eigengenes$eigengenes)," eigengene vectors\n")
    print(head(eigengenes$eigengenes))
    
    # This variable is fixed and used as a reference to indicate the
    #modules used (they are colours but the position within the vector is
    #also relevant)
    centroid.labels <- substring(names(eigengenes$eigengenes),3)
    print("Module colors are")
    print(head(centroid.labels))
    
    k <- length(eigengenes$eigengenes)
    #Centroids must be a matrix with as many colums as centroids,
    #as many rows as samples
    centroids <- createCentroidMatrix(eigengenes$eigengenes)
    
    
    print("We have generated centroids")
    print(head(centroids))
    
    
    #Step 5
    #For storage of all the partitions created during the iterations
    partitions <- list()
    #A partition will be a list of as much elements as genes and for the
    #i-th position it stores the index of the module the ith gene belongs
    #to, and the color can be found in "centroid.labels"
    print("From partition ")
    print(head(partition.in.colors))
    print("We create ")
    new.partition <- match(partition.in.colors, centroid.labels)
    print(head(new.partition))
    names(new.partition) <- centroid.labels[new.partition]
    partitions[[1]] <- new.partition
    
    exchangesVec = c()
    #Launch the iterations
    exchanged.genes = meg + 1
    iteration = 1
    while(exchanged.genes > meg & iteration <= nIterations){
      print(paste0(":) Starting partition ",iteration))
      print(paste0("# of centroids before getting new partition (",ncol(centroids), " gene co-expression modules)"))
      
      new.partition <- apply(datExpr,MARGIN=2,getBestModuleCor,centroids=centroids,
                             signed=(netType == "signed"))
      
      partitions[[iteration + 1]] <- new.partition
      #Get the control values for the new partition
      exchanged.gene.count <- length(getExchangedGenes(partitions[[iteration]],
                                                       partitions[[iteration + 1]]))
      cat("A total of ", exchanged.gene.count,
          " genes moved to another partition\n")
      
      new.partition.in.colors <- centroid.labels[unlist(new.partition)]
      centroids <- getNewCentroids(datExpr,new.partition.in.colors,centroid.labels,min.genes.for.grey)
      
      exchanged.genes = exchanged.gene.count
      exchangesVec = c(exchangesVec, exchanged.gene.count)
      iteration = iteration + 1
    }
    # https://github.com/juanbot/CoExpNets/blob/master/R/main.R
    cat("We finish with",iteration,"iterations\n")
    cat("Last number of gene changes where",exchanged.genes,"\n")
    
    # Convert numeric labels into colors
    dynamicColorsKmeans = new.partition.in.colors # labels2colors(new.partition.in.colors) #assignments44)
    print("Please note we are done with the K-Means step after WGCNA")
  } else {
    print(paste0("Please note that since perform additional K-Means Step is ", performAdditionalKMeansStep, " we only performed WGCNA to get our final gene co-expression modules"))
    dynamicColorsKmeans = dynamicColorsPower
  }
  print(":) Please note we return dynamicColorsPower, which has the final module assignments for our respective genes :)")
  return(dynamicColorsKmeans)
}



updatingFilePathInfoDF <- function(filePathInfoDF,
                                   numSamples, numGenes, powerEstimate,
                                   MEsFilePath, moduleAssignmentsFilePath,
                                   geneModuleMembershipFilePath,
                                   filePathInfoDF_FileNameRData,
                                   filePathInfoDF_FileNameCSV){
  
  filePathInfoDF = rbind(filePathInfoDF,
                         c("numSamples:", numSamples,
                           "integer greater than or equal to 1",
                           paste0("Please note this numSamples is the # of data samples in the gene expression dataset.")))
  
  filePathInfoDF = rbind(filePathInfoDF,
                         c("numGenes:", numGenes,
                           "integer greater than or equal to 1", paste0("Please note this numGenes is the # of genes in the gene expression dataset.")))
  
  filePathInfoDF = rbind(filePathInfoDF,
                         c("powerEstimate:", powerEstimate,
                           "integer greater than or equal to 1 and less than or equal to maxNumPowers",
                           paste0("Please note this powerEstimate is the power between 1 and ", maxNumPowers, " that WGCNA recommends for us to use based on its analysis of Soft-Thresholding and its algorithm.")))
  
  
  
  filePathInfoDF = rbind(filePathInfoDF,
                         c("MEsFilePath:", MEsFilePath,
                           "File Path",
                           paste0("Please note this MEsFilePath is the File Path for the Module Eigengenes from the Co-Expression Gene Modules")))
  
  filePathInfoDF = rbind(filePathInfoDF,
                         c("moduleAssignmentsFilePath:", moduleAssignmentsFilePath,
                           "File Path",
                           paste0("Please note this moduleAssignmentsFilePath is the File Path for the Module Assignments for the Genes from the Co-Expression Gene Modules")))
  
  filePathInfoDF = rbind(filePathInfoDF,
                         c("geneModuleMembershipFilePath:", geneModuleMembershipFilePath,
                           "File Path",
                           paste0("Please note this geneModuleMembershipFilePath is the File Path for the Gene Module Memberships for the Genes from the Co-Expression Gene Modules")))
  
  filePathInfoDF = unique(filePathInfoDF)
  print(":) Please note we have finished updating filePathInfoDF and have save results here:")
  save(filePathInfoDF, file = filePathInfoDF_FileNameRData)
  write.csv(filePathInfoDF, filePathInfoDF_FileNameCSV)
  print(filePathInfoDF_FileNameRData)
  print(filePathInfoDF_FileNameCSV)
  
  
  return(filePathInfoDF)
}



findingSNPsBreakingTFBindingSitesFromGWASDataFunction <- function(gwasDataFile,
                                                                  pValThreshForSNPs,
                                                                  contextForGWASDataSet,
                                                                  updatedGWASDataSetFileNameALL_CSV,
                                                                  updatedGWASDataSetFileNameALL_RData,
                                                                  updatedGWASDataSetFileNameSNPandP_CSV,
                                                                  updatedGWASDataSetFileNameSNPandP_RData,
                                                                  motifbreakrIntervalLength,
                                                                  filePathOfMotifBreakR_Step1,
                                                                  filePathOfMotifBreakR_Step2,
                                                                  versionUsed, genomeUsed,fullDiseaseName,
                                                                  motifBreakrThreshold,
                                                                  methodToUseMotifBreakr,
                                                                  motifBackgroundProbabilitiesVec,
                                                                  motifBreakRParallelComputing_BPPARAM,
                                                                  snpDFPath,
                                                                  snpDFPathRData){
  print(paste0(":) Please note that this function looks at the ", fullDiseaseName, " SNPs (with GWAS P-Value below ",
               pValThreshForSNPs, ") that impact Transcription Factor (TF) binding to TF Binding Sites.  Please note this info is from ", 
               contextForGWASDataSet, " and this package runs MotifBreakR"))
  
  
  # reading in SNPs file:
  gwasSNPsDF = read.csv(gwasDataFile,header = TRUE)
  # SNP, MarkerName, MARKERNAME, snp, snpid, SNPID, SNPId, RS, RSID, RS_NUMBER, RS_NUMBERS
  gwasColsUpperCase = toupper(colnames(gwasSNPsDF))
  # SNP, MARKERNAME, SNPID, RS, RSID, RS_NUMBER, RS_NUMBERS, RSNUMBER,
  # RS_ID, MARKER_NAME, SNP_ID
  # Please note the names of some potential gwas SNP columns
  potentialSNPIDColumnNamesVec = c("SNP", "MARKERNAME", "SNPID", "RS", "RSID", "RS_NUMBER",
                                   "RS_NUMBERS", "RSNUMBER",  "RS_ID", "MARKER_NAME",  "SNP_ID")
  SNPIDcolName = intersect(gwasColsUpperCase, potentialSNPIDColumnNamesVec)
  # please note that \\b is for an exact match once we find the name of the column :)
  SNPIDcol = grep(paste0("\\b", SNPIDcolName, "\\b"), gwasColsUpperCase)
  SNPIDcol
  
  potentialSNPPvalueColumnNames = c('P', 'PVALUE', 'P_VALUE', 'PVAL',
                                    'P_VAL', 'GC_PVALUE', "P-VALUE",
                                    'P-VAL', 'GC-PVALUE', 'SNP-P',
                                    'SNP-PVALUE')
  
  SNPPvalcolName = intersect(gwasColsUpperCase, potentialSNPPvalueColumnNames)
  # please note that \\b is for an exact match once we find the name of the column :)
  SNPPvalcol = grep(paste0("\\b",SNPPvalcolName,"\\b"), gwasColsUpperCase)
  SNPPvalcol
  # Please ensure that the Pvalue column is numeric
  # dataFrame[rows, columns]
  gwasSNPsDF[,SNPPvalcol] = as.numeric(gwasSNPsDF[,SNPPvalcol])
  
  # please note that updatedGWASDF contains all of the original columns as well as the SNPs that have met the threshold requirements
  updatedGWASDF = unique(gwasSNPsDF[which(gwasSNPsDF[,SNPPvalcol] < pValThreshForSNPs),])
  dim(updatedGWASDF)
  print(paste0(":) Please note that we found ", nrow(updatedGWASDF), " unique SNPs for the ", fullDiseaseName, " GWAS Dataset, ",
               " which have a p-value (for the Beta coefficient) less than the threshold ye specified: ",
               pValThreshForSNPs))
  updatedGWASDF$context = rep(contextForGWASDataSet, nrow(updatedGWASDF))
  updatedGWASDF$pValThresholdForSNPs = rep(pValThreshForSNPs, nrow(updatedGWASDF))
  
  write.csv(updatedGWASDF, updatedGWASDataSetFileNameALL_CSV)
  save(updatedGWASDF, pValThreshForSNPs, contextForGWASDataSet, file = updatedGWASDataSetFileNameALL_RData)
  
  rm(gwasSNPsDF)
  # plesae note that this just has the SNPs and the Pvalues
  snpsToUseDF = unique(updatedGWASDF[,c(SNPIDcol, SNPPvalcol)]) # SNP and P <-- p-value
  snpsToUseDF = data.frame(snpsToUseDF)
  colnames(snpsToUseDF) = c("SNP", "P")
  head(snpsToUseDF)
  dim(snpsToUseDF)
  infoSNPs = paste0("SNPs column and P-value only for ", fullDiseaseName, " GWAS")
  snpsToUseDF[]
  write.csv(updatedGWASDF, updatedGWASDataSetFileNameSNPandP_CSV)
  save(updatedGWASDF, infoSNPs, pValThreshForSNPs,
       contextForGWASDataSet, file = updatedGWASDataSetFileNameSNPandP_RData)
  
  
  rsIDsVec = unique(as.vector(snpsToUseDF$SNP))
  print(paste0(":) Please note that we have found: ", length(rsIDsVec), " SNPs"))
  
  rsIDsVec = rsIDsVec[grep("rs", rsIDsVec)]
  length(rsIDsVec) # 94839
  
  
  print(paste0(":) Please note that of those SNPs, we found ", length(rsIDsVec), " SNPs have rsIDs, which is great! :)"))
  data(motifbreakR_motif)
  
  
  snpsFound = 0
  intervalStart = 1
  intervalEnd = motifbreakrIntervalLength + intervalStart - 1
  while (intervalStart < length(rsIDsVec)){
    
    print(paste0(":) please note this SNP interval: ", intervalStart, " to ", intervalEnd))
    snpsUsed = rsIDsVec[intervalStart:intervalEnd]
    snps.mb <- snps.from.rsid(rsid = snpsUsed,
                              dbSNP = snpsDataBase, #SNPlocs.Hsapiens.dbSNP144.GRCh37,
                              search.genome = BSgenome.Hsapiens.UCSC.hg19)
    
    pValName = str_replace_all(pValThreshForSNPs, "-", "minus")
    
    outputFile = paste0(filePathOfMotifBreakR_Step1, intervalStart, "_to_",
                        intervalEnd, "_", fullDiseaseName, "_GWASPVal_",  pValName, ".RData")
    print(paste0(":) please note we write out our output for SNPs ", intervalStart,
                 " to ", intervalEnd, " here: ", outputFile))
    snpsFound = snpsFound + dim(data.frame(snps.mb))[1] # 4756
    print(paste0(":) please note that we found ", snpsFound, " SNPs here!"))
    save(snps.mb, versionUsed, genomeUsed, snpsUsed, pValThreshForSNPs,fullDiseaseName,
         intervalStart, intervalEnd,
         file = outputFile)
    
    print(paste0(":) Please note that we are done with the first part (snps.mb) and now are running motifbreakR :)"))
    
    results <- motifbreakR(snpList = snps.mb, filterp = TRUE,
                           pwmList = motifbreakR_motif,
                           threshold = motifBreakrThreshold, #1e-4,
                           method = methodToUseMotifBreakr, #ic",
                           bkg = motifBackgroundProbabilitiesVec, #c(A=0.25, C=0.25, G=0.25, T=0.25),
                           BPPARAM =  motifBreakRParallelComputing_BPPARAM)#BiocParallel::SnowParam(numWorkersMotifbreakR))
    filePathOfMotifBreakR_Step2
    resultsFilePath = paste0(filePathOfMotifBreakR_Step2, "motifbreakr_results_", fullDiseaseName,"_",
                             intervalStart, "_to_", intervalEnd, ".RData")
    resultsDf = data.frame(results)
    save(results, resultsDf, intervalStart, methodToUseMotifBreakr, versionUsed,
         motifBreakrThreshold, intervalEnd, file = resultsFilePath)
    
    print(paste0(":) Please note that we have finished saving motifbreakR results for this interval (",
                 intervalStart, " to ", intervalEnd, " for ", fullDiseaseName, " here: ", versionUsed))
    # incrementing for the next iteration :)
    intervalStart = intervalEnd + 1
    intervalEnd = motifbreakrIntervalLength + intervalStart - 1
    if (intervalEnd > length(rsIDsVec)){
      intervalEnd = length(rsIDsVec)
      print(paste0(":) please note we are on the last interval, which will be smaller: ",
                   intervalStart, " to ", intervalEnd))
    }
    print("               ")
    print(":) ########################## :)")
    
    
  }
  
  print(paste0(":) please note that we found: ", snpsFound, " SNPs, out of ", length(rsIDsVec)))
  
  
  # please write code to combine all into 1 dataframe :)
  
  #######################################################
  
  
  resultsFilePathsToUseVec = list.files(paste0(filePathOfMotifBreakR_Step2))#, "motifbreakr_results_"))
  
  motifBreakRSNPsDF = data.frame()
  for (i in 1:length(resultsFilePathsToUseVec)){
    print(paste0(":) i = ", i))
    resultPath = paste0(filePathOfMotifBreakR_Step2, resultsFilePathsToUseVec[i])
    snpMotifBreakRDFPart = load(resultPath) #read.csv(resultPath, header = TRUE) #, row.names = 1)
    resultsDf = data.frame(resultsDf)
    if (i == 1){
      motifBreakRSNPsDF = resultsDf
    } else {
      motifBreakRSNPsDF = rbind(motifBreakRSNPsDF, resultsDf)
    }
  }
  motifBreakRSNPsDF = unique(motifBreakRSNPsDF)
  motifBreakRSNPsDF = as.data.frame(motifBreakRSNPsDF)
  str(motifBreakRSNPsDF)
  colnames(motifBreakRSNPsDF)[grep("seqnames", colnames(motifBreakRSNPsDF))] = "chromosome"
  #colnames(motifBreakRSNPsDF)
  head(motifBreakRSNPsDF)
  str(motifBreakRSNPsDF)
  colNum = grep("motifPos", colnames(motifBreakRSNPsDF))
  motifBreakRSNPsDF= motifBreakRSNPsDF[,-colNum]
  head(motifBreakRSNPsDF)
  #write.csv(motifBreakRSNPsDF, snpDFPath) # = paste0(filePathOfMotifBreakR_Final, "motifbreakr_FinalResultsForSNPsWithP_", pValName,".csv")
  #save(motifBreakRSNPsDF, file = snpDFPathRData)# #snpDFPathRData = paste0(filePathOfMotifBreakR_Final, "motifbreakr_FinalResultsForSNPsWithP_", pValName,".RData")
  dim(motifBreakRSNPsDF)
  
  updatedGWASDF = read.csv(updatedGWASDataSetFileNameALL_CSV, header = TRUE)
  colnames(updatedGWASDF)[grep("\\bSNP\\b", colnames(updatedGWASDF))] = "SNP_id"
  head(updatedGWASDF)
  
  motifBreakRSNPsDF = motifBreakRSNPsDF %>% inner_join(updatedGWASDF)
  head(motifBreakRSNPsDF)
  
  write.csv(motifBreakRSNPsDF, snpDFPath) # = paste0(filePathOfMotifBreakR_Final, "motifbreakr_FinalResultsForSNPsWithP_", pValName,".csv")
  save(motifBreakRSNPsDF, file = snpDFPathRData)# #snpDFPathRData = paste0(filePathOfMotifBreakR_Final, "motifbreakr_FinalResultsForSNPsWithP_", pValName,".RData")
  head(motifBreakRSNPsDF)
  print(paste0(":) please note that motifBreakRSNPsDF has these dimensions: ", dim(motifBreakRSNPsDF)))
  print(head(motifBreakRSNPsDF))
  return((motifBreakRSNPsDF))
}



gettingInformationForGeneExpressionGeneRegulatoryNetworks <- function(tfsInGeneExpressionDataRData,
                                                                      entrezMappingFilePath,
                                                                      tfsInGeneExpressionDataCSV,
                                                                      filePathOfGeneExpressionDataUsedByWGCNA,
                                                                      inputGeneExpressionDataFilePath,
                                                                      geneToEntrezIDPathForGeneExpressionDataCSV,
                                                                      geneToEntrezIDPathForGeneExpressionDataRData,
                                                                      tfsDF,
																	  filePathOfTFsToUse){
  print(":) Please note that this function gets information on the geneAndEntrezIDMappingDF, tfsToUse, and allsamples (Gene expression data), which will be used by RTN, GENIE3, and TReNA Ensemble Solver :)")
  print(paste0(":) Please note that an RData File with this information will be available here: ", tfsInGeneExpressionDataRData))
  geneInfoFilePath = entrezMappingFilePath
  
  head(infoDF)
  
  allsamples <- read.csv(filePathOfGeneExpressionDataUsedByWGCNA,
                         header = TRUE, row.names = 1)
  
  dim(allsamples)
  geneExpressionSamplesVec = row.names(allsamples)
  annot = colnames(allsamples)
  if (is.null(filePathOfTFsToUse)){
    print("please note that we use adsnpheno TFs list to find set of TFs here:")
    geneAndEntrezIDMappingDF = geneExpressionDataSetToEntrezIDMapping(inputGeneExpressionDataFilePath,
                                                                      geneToEntrezIDPathForGeneExpressionDataCSV,
                                                                      geneToEntrezIDPathForGeneExpressionDataRData)
    tfsToUse = identifyingTFsInGeneExpressionDataSet(geneAndEntrezIDMappingDF,
                                                     allsamples, tfsDF)
    write.csv(tfsToUse, tfsInGeneExpressionDataCSV)
    save(tfsToUse, annot, allsamples, geneExpressionSamplesVec, file = tfsInGeneExpressionDataRData)
    
  } else if (useTFsFoundByADSNPheno){
    print("please note that we use adsnpheno TFs list to find set of TFs here:")
    geneAndEntrezIDMappingDF = geneExpressionDataSetToEntrezIDMapping(inputGeneExpressionDataFilePath,
                                                                      geneToEntrezIDPathForGeneExpressionDataCSV,
                                                                      geneToEntrezIDPathForGeneExpressionDataRData)
    tfsToUse = identifyingTFsInGeneExpressionDataSet(geneAndEntrezIDMappingDF,
                                                     allsamples, tfsDF)
    write.csv(tfsToUse, tfsInGeneExpressionDataCSV)
    save(tfsToUse, annot, allsamples, geneExpressionSamplesVec, file = tfsInGeneExpressionDataRData)
    
  } else {
    print(paste0("please note that we use custom list of TFs ye gave: (", tfsUsed, ") to find set of TFs here:"))
    
    tfsToUse = read.csv(filePathOfTFsToUse, header = TRUE)
    tfsToUse = as.character(tfsToUse[,grep(colnames(tfsToUse), "tfName")])
    print(paste0("please note that your custom list had: ", length(tfsToUse), "TFs"))#
    tfsToUse = intersect(row.names(allsamples), tfsToUse)
    print(paste0("of those TFs ye gave, we matched by the tfName and the corresponding geneName in the gene expression data, please note that we found ", length(tfsToUse), " matching TFs in the gene expression dataset."))#
    save(tfsToUse, annot, allsamples, geneExpressionSamplesVec, file = tfsInGeneExpressionDataRData)
    
  }
  outputsList = list()
  outputsList$tfsToUse = tfsToUse
  outputsList$allsamples = allsamples
  outputsList$geneExpressionSamplesVec = geneExpressionSamplesVec
  outputsList$tfsInGeneExpressionDataRData = tfsInGeneExpressionDataRData
  return(outputsList)
}



##########################################################################################
# PLEASE NOTE THIS DEMO FOR 1 CHROMOSOME:
# user can have the option to create an interactionsDF or to run scgrnom to get it
# caseInfo


pleaseBuildChromatinGeneRegulatoryNetwork <- function(numOfCoresToUseForScgrnom,
                                                      haveChromatinInteractionRegulatoryNetwork,
                                                      haveEnhancerAndPromoterInteractionsData,
                                                      haveEnhancersAndPromotersHICData,
                                                      bodyRegion, otherInfo, folderName,
                                                      chromatinInteraction_Step2_FilePath,
                                                      chromatinInteraction_Step2_FilePath_InitialRDataObjects,
                                                      chromatinInteraction_Step2_FilePath_CleanedRDataObjects,
                                                      chromatinInteraction_Step2_FilePath_CleanedCSVObjects,
                                                      includeTheIndividualTFsInGroup,
                                                      pleasePrintProgressAfterNumIterations_GetTF,
                                                      chromatinInteraction_Step1ADSNPheno_FilePath,
                                                      maxBasesUpstreamForPromoterLength, enhancerChromatinFilePath,
                                                      geneIDColumnNameEnhancerDF, enhancerRegionCol,
                                                      pleaseOrganizeEnhancerRegions,
                                                      chromatinInteractionRegulatoryNetworkFilePath,
                                                      enhAndPromoterInteractionsFilePath,
                                                      chromatinInteractionDataFilePath,
                                                      enhancerAndPromoterInteraction){
  
  print(paste0(":) Please note that this function helps us create a Gene Regulatory Network using Chromatin (Epigenetics) data on interactions between promoters and enhancers near target genes."))
  
  if (haveChromatinInteractionRegulatoryNetwork){ # please we already have this network (from another external source) and therefore do not need to do scgrnom at all.
    chromRegulatoryNetwork = read.csv(chromatinInteractionRegulatoryNetworkFilePath, header = TRUE)
    caseInfo = "haveChromatinRegulatoryNetwork_Done"
  } else if (haveEnhancerAndPromoterInteractionsData){ # please do not perform step 1 as we already have this interactions data
    interactionsDF_ALL = read.csv(enhAndPromoterInteractionsFilePath, header = TRUE)
    caseInfo = "haveEnhancerAndPromoterInteractionsData_Step1_Done"
    
  }  else if(haveEnhancersAndPromotersHICData){ # please perform scgrnom part 1 and part 2
    
    enhancer <- processEnhancerFile(chromatinInteractionDataFilePath)
    
    hic_interaction = processHICInteractionFile(chromatinInteractionDataFilePath)
    save(enhancer, hic_interaction,
         file = paste0(enhancerAndPromoterInteraction, "chromatinInteractionData.RData"))
    print(paste0(":) Please note that the enhancer and promoter interaction data for ", folderName,
                 " for scgrnom has been stored here: ", enhancerAndPromoterInteraction))
    
    caseInfo = "haveHIC_PromoterAndEnhancerData"
    
  } else if (haveEnhancersInteractionData){ # we perform adsnpheno for getting the promoter info
    caseInfo = "haveOnlyEnhancerData"
  } else {
    print(paste0(":(  Please note that ye did not meet any of the criteria for the Chromatin Interaction Regulatory Network.  At a minimum, please find enhancer location data for the region of interest ",
                 bodyRegion, "in the form: (chr, start, end) please.  Thank you! :)"))
    caseInfo = "ERROR! :("
  }
  
  
  print(paste0(":) Please note that the caseInfo is: ", caseInfo))
  
  # Please note that since the chromosomes are usually smaller in size the higher the number of the chromosome,
  # we start with the later chromosomes to ensure those work :)
  chromosomeProblems = c() # please note this stores the chromosomes that gave us problems
  
  for (i in 1:length(chromosomeVec)){
    if (caseInfo == "haveChromatinRegulatoryNetwork_Done"){
      print(paste0("Please exit loop since caseInfo: ", caseInfo))
      break
    }
    
    tryCatch(
      {
        chromName = chromosomeVec[(length(chromosomeVec) - i) + 1] # we start with the later chromosomes to ensure those work :)
        
        print(paste0(":) please note that i = ", i, " of ", length(chromosomeVec),
                     " and the corresponding chromosome we are getting interaction information for is: ", chromName))
        
        if (caseInfo == "haveEnhancerAndPromoterInteractionsData_Step1_Done"){
          df_chrInteractions =  pleaseGetInteractionsDFForChromosome_haveEnhancerAndPromoterInteractionsData_Step1_Done(chromName, caseInfo, enhAndPromoterInteractionsFilePath)
        } else if (caseInfo == "haveHIC_PromoterAndEnhancerData") {
          df_chrInteractions =  pleaseGetInteractionsDFForChromosome_haveHIC_PromoterAndEnhancerData(chromName, chromatinInteraction_Step1_FilePath,
                                                                                                     hic_interaction, enhancers, bodyRegion, caseInfo)
          
        } else if (caseInfo == "haveOnlyEnhancerData") {
          df_chrInteractions = pleaseGetInteractionsDFForChromosome_haveOnlyEnhancerData(chromName, chromatinInteraction_Step1ADSNPheno_FilePath,
                                                                                         bodyRegion,maxBasesUpstreamForPromoterLength, enhancerChromatinFilePath, geneIDColumnNameEnhancerDF, enhancerRegionCol,
                                                                                         pleaseOrganizeEnhancerRegions, caseInfo)
        }
        
        print(paste0(":) Please note we used the respective interactionFunction for caseInfo: ", caseInfo, " for chromosome: ", chromName, " and have our dataset! :)"))
        print(paste0("Please note the dimensions: df_chrInteractions = ", dim(df_chrInteractions)))# now we call the getTF function by scGRNom to get the regulatory information for this chromosome ", chromName, " for ", bodyRegion))
        
        print(paste0("Please note that now we call the getTF function by scGRNom to get the regulatory information for this chromosome ", chromName, " for ", bodyRegion))
        organizedGetTFInfoForChromDF_ForChrom = pleaseGetTFsFromPromoterAndInteractionDataStep2(df_chrInteractions,
                                                                                                chromName, numOfCoresToUseForScgrnom,
                                                                                                bodyRegion, otherInfo, folderName, chromatinInteraction_Step2_FilePath,
                                                                                                chromatinInteraction_Step2_FilePath_InitialRDataObjects,
                                                                                                chromatinInteraction_Step2_FilePath_CleanedRDataObjects,
                                                                                                chromatinInteraction_Step2_FilePath_CleanedCSVObjects,
                                                                                                includeTheIndividualTFsInGroup,
                                                                                                pleasePrintProgressAfterNumIterations_GetTF)
      },
      error = function(cond){
        message(paste0("Please note we had an error for this chromosome: ", chromName))
        message(cond)
        chromosomeProblems = c(chromosomeProblems, chromName)
        # next
      })
  }
  
  scgrnGetTFsForEachChromosomesChromatinNetworkVec = list.files(chromatinInteraction_Step2_FilePath_CleanedCSVObjects, pattern = "scgrnomGetTFs_")
  scgrnGetTFsForEachChromosomesChromatinNetworkVec = paste0(chromatinInteraction_Step2_FilePath_CleanedCSVObjects, scgrnGetTFsForEachChromosomesChromatinNetworkVec)
  print(":) Please note the files are here (organized for TF to TG interactions for Chromatin Regulatory Network for a given Chromosome: ")
  print(scgrnGetTFsForEachChromosomesChromatinNetworkVec)
  return(scgrnGetTFsForEachChromosomesChromatinNetworkVec)
}



# Building Transcriptional Regulatory Network Using RTN
# RTN:
pleaseGet_RTN_GeneRegulatoryNetwork <- function(tfsInGeneExpressionDataRData,
                                                numRTNPermutations,
                                                oldRDataPart,
                                                newRDataPart,
                                                numRoundingDigits,
                                                RTN_rDataPart,
                                                csvPartRTNFinal,
                                                csvPartRTNFinal_NoDate){
  print(paste0(":) please note the # of rows in gene expression data set is: ", nrow(allsamples)))
  
  load(tfsInGeneExpressionDataRData) # loading info on TFsToUse
  rtni <- tni.constructor(expData = as.matrix(allsamples),
                          regulatoryElements = tfsToUse)#,
  # -Preprocessing for input data...
  # --Checking 'regulatoryElements' in 'rowAnnotation'...
  # --Checking 'expData'...
  #   -Preprocessing complete!
  #rowAnnotation = annot)
  
  rtni <- tni.permutation(rtni, nPermutations = numRTNPermutations)
  # -Performing permutation analysis...
  # --For # regulons...
  infoRTN = "old (after performing tni.permutation only)"
  save(rtni, infoRTN, numRTNPermutations, file = oldRDataPart)
  print(paste0(":) Please note we wrote the rdata output for RTN here for ", numRTNPermutations, " permutations: ", newRDataPart))
  rtni <- tni.bootstrap(rtni)
  rtni <- tni.dpi.filter(rtni)
  infoRTN = "final (after performing tni.permutation, bootstrap, and dpi.filter)"
  
  save(rtni, infoRTN, file = newRDataPart)
  
  tni.regulon.summary(rtni)
  regulons <- tni.get(rtni, what = "regulons.and.mode", idkey = "ID")
  regulons
  #numRoundingDigits = 4
  regulonsNames = names(regulons)
  regulonsDF_list = list()
  counter = 1
  for (i in 1:length(regulonsNames)){
    regulon = regulonsNames[i]
    regulatedGenes = regulons[[regulon]] #[regulon]]
    if (length(regulatedGenes) > 0){
      
      regulatedGenes
      mini = data.frame(regulatedGenes)
      targetGene = row.names(mini)
      mutualInfo = mini[,1]
      mutualInfoName = paste0("RTN (Mutual Info: ",
                              round(mutualInfo,
                                    numRoundingDigits), ")")
      dfToAdd = data.frame(cbind(regulon, targetGene, "RTN", mutualInfoName))
      colnames(dfToAdd) = c("TF", "RegulatedGene", "Source", "Info")
      regulonsDF_list[[counter]] = dfToAdd
      counter = counter + 1
    }
  }
  
  rtnRegulonsDF = do.call("rbind", regulonsDF_list)
  rtnRegulonsDF$CombinedName = paste0(rtnRegulonsDF$TF, " || ", rtnRegulonsDF$RegulatedGene)
  head(rtnRegulonsDF)
  dim(rtnRegulonsDF)
  print(paste0(":) Please note we wrote the csv output for RTN here: ", newCsvPart))
  grnMethod = "RTN"
  infoRTN = "final (after performing tni.permutation, bootstrap, and dpi.filter, and adjusting to make Data Frame)"
  
  save(rtni, infoRTN, rtnRegulonsDF, allsamples, tfsToUse, grnMethod, numRTNPermutations, file = RTN_rDataPart)
  
  print(paste0(":) Please note we have finished with RTN :) and wrote out results here: ", RTN_rDataPart))
  
  write.csv(rtnRegulonsDF, csvPartRTNFinal)
  write.csv(rtnRegulonsDF, csvPartRTNFinal_NoDate)
  
  print(paste0(":) Please note we have finished running RTN for TF-TG relationships and we have our final output in our rtnRegulonsDF"))
  print(paste0(" Please note dimensions: ", RTN_rDataPart))
  print(head(RTN_rDataPart))
  return(rtnRegulonsDF)
  
}


# Please note that the code below on Decision Curve Analysis is from the Memorial
# Sloan Kettering Cancer Center:
dca <- function(data, outcome, predictors, xstart=0.01, xstop=0.99, xby=0.01,
                ymin=-0.05, probability=NULL, harm=NULL,graph=TRUE, intervention=FALSE,
                interventionper=100, smooth=FALSE,loess.span=0.10, context = "") {

  print(paste0(":) Please note that this code on DCA (Decision Curve Analysis) is directly from the Memorial Sloan Kettering Cancer Center (MSKCC): https://www.mskcc.org/departments/epidemiology-biostatistics/biostatistics/decision-curve-analysis"))
  # LOADING REQUIRED LIBRARIES
  require(stats)

  # data MUST BE A DATA FRAME
  if (class(data)!="data.frame") {
    stop("Input data must be class data.frame")
  }

  #ONLY KEEPING COMPLETE CASES
  data=data[complete.cases(data[append(outcome,predictors)]),append(outcome,predictors)]

  # outcome MUST BE CODED AS 0 AND 1
  if (max(data[[outcome]])>1 | min(data[[outcome]])<0) {
    stop("outcome cannot be less than 0 or greater than 1")
  }
  # xstart IS BETWEEN 0 AND 1
  if (xstart<0 | xstart>1) {
    stop("xstart must lie between 0 and 1")
  }

  # xstop IS BETWEEN 0 AND 1
  if (xstop<0 | xstop>1) {
    stop("xstop must lie between 0 and 1")
  }

  # xby IS BETWEEN 0 AND 1
  if (xby<=0 | xby>=1) {
    stop("xby must lie between 0 and 1")
  }

  # xstart IS BEFORE xstop
  if (xstart>=xstop) {
    stop("xstop must be larger than xstart")
  }

  #STORING THE NUMBER OF PREDICTORS SPECIFIED
  pred.n=length(predictors)

  #IF probability SPECIFIED ENSURING THAT EACH PREDICTOR IS INDICATED AS A YES OR NO
  if (length(probability)>0 & pred.n!=length(probability)) {
    stop("Number of probabilities specified must be the same as the number of predictors being checked.")
  }

  #IF harm SPECIFIED ENSURING THAT EACH PREDICTOR HAS A SPECIFIED HARM
  if (length(harm)>0 & pred.n!=length(harm)) {
    stop("Number of harms specified must be the same as the number of predictors being checked.")
  }

  #INITIALIZING DEFAULT VALUES FOR PROBABILITES AND HARMS IF NOT SPECIFIED
  if (length(harm)==0) {
    harm=rep(0,pred.n)
  }
  if (length(probability)==0) {
    probability=rep(TRUE,pred.n)
  }


  #CHECKING THAT EACH probability ELEMENT IS EQUAL TO YES OR NO,
  #AND CHECKING THAT PROBABILITIES ARE BETWEEN 0 and 1
  #IF NOT A PROB THEN CONVERTING WITH A LOGISTIC REGRESSION
  for(m in 1:pred.n) {
    if (probability[m]!=TRUE & probability[m]!=FALSE) {
      stop("Each element of probability vector must be TRUE or FALSE")
    }
    if (probability[m]==TRUE & (max(data[predictors[m]])>1 | min(data[predictors[m]])<0)) {
      stop(paste(predictors[m],"must be between 0 and 1 OR sepcified as a non-probability in the probability option",sep=" "))
    }
    if(probability[m]==FALSE) {
      model=NULL
      pred=NULL
      model=glm(data.matrix(data[outcome]) ~ data.matrix(data[predictors[m]]), family=binomial("logit"))
      pred=data.frame(model$fitted.values)
      pred=data.frame(pred)
      names(pred)=predictors[m]
      data=cbind(data[names(data)!=predictors[m]],pred)
      print(paste(predictors[m],"converted to a probability with logistic regression. Due to linearity assumption, miscalibration may occur.",sep=" "))
    }
  }

  # THE PREDICTOR NAMES CANNOT BE EQUAL TO all OR none.
  if (length(predictors[predictors=="all" | predictors=="none"])) {
    stop("Prediction names cannot be equal to all or none.")
  }

  #########  CALCULATING NET BENEFIT   #########
  N=dim(data)[1]
  event.rate=colMeans(data[outcome])

  # CREATING DATAFRAME THAT IS ONE LINE PER THRESHOLD PER all AND none STRATEGY
  nb=data.frame(seq(from=xstart, to=xstop, by=xby))
  names(nb)="threshold"
  interv=nb

  nb["all"]=event.rate - (1-event.rate)*nb$threshold/(1-nb$threshold)
  nb["none"]=0

  # CYCLING THROUGH EACH PREDICTOR AND CALCULATING NET BENEFIT
  for(m in 1:pred.n){
    for(t in 1:length(nb$threshold)){
      # COUNTING TRUE POSITIVES AT EACH THRESHOLD
      tp=mean(data[data[[predictors[m]]]>=nb$threshold[t],outcome])*sum(data[[predictors[m]]]>=nb$threshold[t])
      # COUNTING FALSE POSITIVES AT EACH THRESHOLD
      fp=(1-mean(data[data[[predictors[m]]]>=nb$threshold[t],outcome]))*sum(data[[predictors[m]]]>=nb$threshold[t])
      #setting TP and FP to 0 if no observations meet threshold prob.
      if (sum(data[[predictors[m]]]>=nb$threshold[t])==0) {
        tp=0
        fp=0
      }

      # CALCULATING NET BENEFIT
      nb[t,predictors[m]]=tp/N - fp/N*(nb$threshold[t]/(1-nb$threshold[t])) - harm[m]
    }
    interv[predictors[m]]=(nb[predictors[m]] - nb["all"])*interventionper/(interv$threshold/(1-interv$threshold))
  }

  # CYCLING THROUGH EACH PREDICTOR AND SMOOTH NET BENEFIT AND INTERVENTIONS AVOIDED
  for(m in 1:pred.n) {
    if (smooth==TRUE){
      lws=loess(data.matrix(nb[!is.na(nb[[predictors[m]]]),predictors[m]]) ~ data.matrix(nb[!is.na(nb[[predictors[m]]]),"threshold"]),span=loess.span)
      nb[!is.na(nb[[predictors[m]]]),paste(predictors[m],"_sm",sep="")]=lws$fitted

      lws=loess(data.matrix(interv[!is.na(nb[[predictors[m]]]),predictors[m]]) ~ data.matrix(interv[!is.na(nb[[predictors[m]]]),"threshold"]),span=loess.span)
      interv[!is.na(nb[[predictors[m]]]),paste(predictors[m],"_sm",sep="")]=lws$fitted
    }
  }

  # PLOTTING GRAPH IF REQUESTED
  if (graph==TRUE) {
    require(graphics)

    # PLOTTING INTERVENTIONS AVOIDED IF REQUESTED
    if(intervention==TRUE) {
      # initialize the legend label, color, and width using the standard specs of the none and all lines
      legendlabel <- NULL
      legendcolor <- NULL
      legendwidth <- NULL
      legendpattern <- NULL

      #getting maximum number of avoided interventions
      ymax=max(interv[predictors],na.rm = TRUE)

      #INITIALIZING EMPTY PLOT WITH LABELS
      plot(x=nb$threshold, y=nb$all, type="n" ,xlim=c(xstart, xstop), ylim=c(ymin, ymax), xlab="Threshold probability", ylab=paste("Net reduction in interventions per",interventionper,"patients"))

      #PLOTTING INTERVENTIONS AVOIDED FOR EACH PREDICTOR
      for(m in 1:pred.n) {
        if (smooth==TRUE){
          lines(interv$threshold,data.matrix(interv[paste(predictors[m],"_sm",sep="")]),col=m,lty=2)
        } else {
          lines(interv$threshold,data.matrix(interv[predictors[m]]),col=m,lty=2)
        }

        # adding each model to the legend
        legendlabel <- c(legendlabel, predictors[m])
        legendcolor <- c(legendcolor, m)
        legendwidth <- c(legendwidth, 1)
        legendpattern <- c(legendpattern, 2)
      }
    } else {
      # PLOTTING NET BENEFIT IF REQUESTED

      # initialize the legend label, color, and width using the standard specs of the none and all lines
      legendlabel <- c("None", "All")
      legendcolor <- c(17, 8)
      legendwidth <- c(2, 2)
      legendpattern <- c(1, 1)

      #getting maximum net benefit
      ymax=max(nb[names(nb)!="threshold"],na.rm = TRUE)

      # inializing new benfit plot with treat all option
      plot(x=nb$threshold, y=nb$all, type="l", col=8, lwd=2 ,xlim=c(xstart, xstop), ylim=c(ymin, ymax),
           xlab="Threshold probability", ylab="Net benefit", main = context)
      # adding treat none option
      lines(x=nb$threshold, y=nb$none,lwd=2)
      #PLOTTING net benefit FOR EACH PREDICTOR
      for(m in 1:pred.n) {
        if (smooth==TRUE){
          lines(nb$threshold,data.matrix(nb[paste(predictors[m],"_sm",sep="")]),col=m,lty=2)
        } else {
          lines(nb$threshold,data.matrix(nb[predictors[m]]),col=m,lty=2)
        }
        # adding each model to the legend
        legendlabel <- c(legendlabel, predictors[m])
        legendcolor <- c(legendcolor, m)
        legendwidth <- c(legendwidth, 1)
        legendpattern <- c(legendpattern, 2)
      }
    }
    # then add the legend
    legend("topright", legendlabel, cex=0.8, col=legendcolor, lwd=legendwidth, lty=legendpattern)

  }

  #RETURNING RESULTS
  results=list()
  results$N=N
  results$predictors=data.frame(cbind(predictors,harm,probability))
  names(results$predictors)=c("predictor","harm.applied","probability")
  results$interventions.avoided.per=interventionper
  results$net.benefit=nb
  results$interventions.avoided=interv

  return(results)

}





pleaseGetInitialFilePathDerivationsAndInfoObjectsForWGCNA <- function(outputPathNameADSNPhenoOutputs, inputGeneExpressionDataFilePath,
disease, tissueName,
bodyRegion, tfsUsed,
performAdditionalKMeansStep,
phenotypesFilePath, log2transformInputData,
scaleInputData, includeTimeStamp, numberOfRoundingDigits, netType){
print(":) Please note that this function takes in initial user input and helps us store the file paths and other key variables that we may need for our analysis.") 
numOfRoundingDigits = numberOfRoundingDigits

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

if(isTRUE(performAdditionalKMeansStep)){
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
if (isTRUE(performAdditionalKMeansStep)){
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

diseaseTraits = read.csv(phenotypesFilePath, header = TRUE)

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
#if (includeTimeStamp == "TRUE"){
  outputAddOnWithDate = paste0(outputAddOn, returnCurrentDateToFilePath())
  updatedOutputAddOn_withDate = paste0(finalGeneExpressionPreparationFilePath, "//finalGeneExpressionDataExpr_",
                                       disease, "_", bodyRegion, outputAddOnWithDate, ".csv" )

  updatedOutputAddOnAndDatExpr_withDate = paste0(finalGeneExpressionPreparationFilePath, "//finalGeneExpressionDatExpr_",
                                                 disease, "_", bodyRegion, outputAddOnWithDate, ".csv" )

  geneExpressionOutputNameCSV_withDate = str_replace_all(updatedOutputAddOn_withDate, " ", "")
geneExpressionOutputNameCSV_withDate = str_replace_all(geneExpressionOutputNameCSV_withDate,  " ", "")
updatedOutputAddOnAndDatExpr_withDate = str_replace_all(updatedOutputAddOnAndDatExpr_withDate,  " ", "")
geneExpressionOutputNameRData_withDate = str_replace_all(geneExpressionOutputNameCSV_withDate, ".csv", ".RData") #paste0("_", diseaseName, "_", bodyRegion, outputAddOn, ".RData")

#}

# please get rid of any possible spaces in the filePath
inputGeneExpressionDataFilePath = str_replace_all(inputGeneExpressionDataFilePath,  " ", "")
geneToEntrezIDPathForGeneExpressionDataCSV = str_replace_all(geneToEntrezIDPathForGeneExpressionDataCSV,  " ", "")
geneToEntrezIDPathForGeneExpressionDataRData = str_replace_all(geneToEntrezIDPathForGeneExpressionDataRData,  " ", "")
geneExpressionOutputNameCSV = str_replace_all(updatedOutputAddOn, " ", "")
geneExpressionOutputNameRData = str_replace_all(geneExpressionOutputNameCSV, ".csv", ".RData") #paste0("_", diseaseName, "_", bodyRegion, outputAddOn, ".RData")

geneExpressionOutputNameCSV = str_replace_all(geneExpressionOutputNameCSV,  " ", "")
originalGeneFilteredCSV = str_replace_all(originalGeneFilteredCSV,  " ", "")
originalGeneFilteredRData = str_replace_all(originalGeneFilteredRData,  " ", "")
geneExpressionOutputNameRData = str_replace_all(geneExpressionOutputNameRData,  " ", "")
updatedOutputAddOnAndDatExpr = str_replace_all(updatedOutputAddOnAndDatExpr,  " ", "")
print(paste0("geneExpressionOutputNameCSV_withDate: ", geneExpressionOutputNameCSV_withDate))
print(paste0("geneExpressionOutputNameRData_withDate: ", geneExpressionOutputNameRData_withDate))

print(paste0("geneExpressionOutputNameCSV: ", geneExpressionOutputNameCSV))
print(paste0("geneExpressionOutputNameRData: ", geneExpressionOutputNameRData))

region = folderName
regionName = folderName


# file path for the Topological Overlap Matrix Info
tomFilePath = paste0(wgcnaAndKMeansOutputPath, "//TopologicalOverlapMatrix")
print(paste0("tomFilePath: ", tomFilePath))
dir.create(tomFilePath, showWarnings = FALSE)


initialRDataFilePath_FilePathDerivationsAndInfoObjectsForWGCNA = paste0(mainOutputFolder, "//initialRDataFilePath_FilePathDerivationsAndInfoObjectsForWGCNA.RData")


# please note that we also save the data objects from WGCNA here so we can call them in future functions without needing to know the power :)
wgcnaWithKmeansAllObjectsExceptTOM_SimpleFileName = paste0(mainOutputFolder, wgcnaToAdd, bodyRegion, "_wgcnaSimplePathALLKeyObjectsExceptTOM",outputAddOn,".RData")



#if (includeTimeStamp == "TRUE"){		 
save(fullDiseaseName, folderName, mainOutputFolder,
region, regionName, tomFilePath,
wgcnaWithKmeansAllObjectsExceptTOM_SimpleFileName,
filePathInfoDF_FileNameRData,
filePathInfoDF_FileNameCSV,
tfsInGeneExpressionDataCSV,
tfsInGeneExpressionDataRData,
wgcnaAndKMeansOutputPath,
wgcnaToAdd,
adsnphenoDirectory,
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
numberOfRoundingDigits, numOfRoundingDigits,
file = initialRDataFilePath_FilePathDerivationsAndInfoObjectsForWGCNA)
#} else {
#save(fullDiseaseName, folderName, mainOutputFolder,
#region, regionName, tomFilePath,
#wgcnaWithKmeansAllObjectsExceptTOM_SimpleFileName,
#filePathInfoDF_FileNameRData,
#filePathInfoDF_FileNameCSV,
#tfsInGeneExpressionDataCSV,
#tfsInGeneExpressionDataRData,
#wgcnaAndKMeansOutputPath,
#wgcnaToAdd,
#adsnphenoDirectory,
#datTraitsRDataFilePath,
#powerRelatedFilePathsRData,
#geneExpressionPreparationFilePath,
#finalGeneExpressionPreparationFilePath,
#entrezMappingFilePath,
#softThresholdPowerFilePath,
#tomFilePath, wgcnaFilePath,
#wgcnaWithKMeansAfterFilePath,
#pathForWGCNA, MEsFilePath,
#moduleAssignmentsFilePath,
#geneModuleMembershipFilePath,
#allObjectsFilePath,
#phenotypeAssociationFilePath,
#modulePhenoAssociationFilePath,
#genePhenoAssociationFilePath,
#diseaseTraits, newDate, outputAddOn, dataScaling,
#dataScalingOutputMini, geneToEntrezIDPathForGeneExpressionDataRData,
#geneToEntrezIDPathForGeneExpressionDataCSV,
#originalGeneFilteredCSV,
#originalGeneFilteredRData, updatedOutputAddOn,
#filePathOfGeneExpressionDataUsedByWGCNA,
#updatedOutputAddOnAndDatExpr,
#sftThresholdingOutputNameCSV,
#sftThresholdingOutputNameRData,
#filePathInfoDF, inputGeneExpressionDataFilePath,
#geneExpressionOutputNameCSV,
#geneExpressionOutputNameRData,
#numberOfRoundingDigits, numOfRoundingDigits,
#file = initialRDataFilePath_FilePathDerivationsAndInfoObjectsForWGCNA)


#}


print(paste0(":) Please note that we wrote a lot of the key files for our WGCNA Gene Co-Expression Network Analysis to an RData File :)"))
print(paste0("initialRDataFilePath_FilePathDerivationsAndInfoObjectsForWGCNA:", initialRDataFilePath_FilePathDerivationsAndInfoObjectsForWGCNA))
return(initialRDataFilePath_FilePathDerivationsAndInfoObjectsForWGCNA)

}



pleaseRunSourceFilesFromDataAndCodeFolders <- function(dataFolder, codeFolder){
print(paste0(":) Please note that this function takes in the dataFolder (", dataFolder, ") and the codeFolder (", codeFolder, 
") and then loads generalDataForPipeline, Python Scripts, packages, and functions we need :)"))

generalDataForPipeline = paste0(dataFolder, "generalDataForPipeline.RData")
load(generalDataForPipeline)
pathForPythonCode = paste0(codeFolder, "setupScripts//BuildingFullAndSubNetworksPythonCode.py")
filePathDerivationsSourceCode = paste0(codeFolder, "setupScripts//filePathDerivations.R")
packagesNeededSourceCode = paste0(codeFolder, "setupScripts//packagesNeeded.R")
functionsNeededSourceCode = paste0(codeFolder, "setupScripts//functionsNeeded.R")
userInputsSourceCode = paste0(codeFolder, "pythonCode//BuildingFullAndSubNetworksPythonCode.py")
generalDataForPipeline = paste0(dataFolder, "generalDataForPipeline.RData")

source(packagesNeededSourceCode)
source(functionsNeededSourceCode)
load(generalDataForPipeline)
source(filePathDerivationsSourceCode)
source_python(pathForPythonCode)

print(paste0(":) Please note: packagesNeededSourceCode: ", packagesNeededSourceCode))
print(paste0(":) Please note: functionsNeededSourceCode: ", functionsNeededSourceCode))
print(paste0(":) Please note: generalDataForPipeline: ", generalDataForPipeline))
print(paste0(":) Please note: filePathDerivationsSourceCode: ", filePathDerivationsSourceCode))
print(paste0(":) Please note: pathForPythonCode: ", pathForPythonCode))

fileNamesList = list()
fileNamesList$packagesNeededSourceCode = packagesNeededSourceCode
fileNamesList$functionsNeededSourceCode = functionsNeededSourceCode
fileNamesList$generalDataForPipeline = generalDataForPipeline
fileNamesList$filePathDerivationsSourceCode = filePathDerivationsSourceCode

fileNamesList$pathForPythonCode = pathForPythonCode
return(fileNamesList)

}



pleaseCreateModuleEnrichmentsFileNames <- function(regionName, pathForWGCNA, performMESHEnrichment, performMolSigDBEnrichment){
print(paste0(":) Please note that this function creates file names and file paths for the Module Enrichments :)"))
library(DOSE)
library(meshes)
library(MeSH.Hsa.eg.db)
library(clusterProfiler)
# "NeuroscienceProjectAxons"
#regionName = folderName

outputNamesList = list() # please create a list to return everything
moduleEnrichmentsFilePath = paste0(pathForWGCNA, "//ModuleEnrichments")
dir.create(moduleEnrichmentsFilePath, showWarnings = FALSE)
outputNamesList$moduleEnrichmentsFilePath = moduleEnrichmentsFilePath


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
  dir.create(snpsOutputPath, showWarnings = FALSE)
  filePathOfGWAS  = paste0(snpsOutputPath, "//", "GWAS_Data")
  dir.create(filePathOfGWAS, showWarnings = FALSE)
  filePathOfMotifBreakR_Step1  = paste0(snpsOutputPath, "//", "MotifbreakR_Part1//")
  filePathOfMotifBreakR_Step2  = paste0(snpsOutputPath, "//", "MotifbreakR_Part2_BrokenTFBS//")
  filePathOfMotifBreakR_Final  = paste0(snpsOutputPath, "//", "MotifbreakR_FinalResults//")

  pValName = str_replace(pValThreshForSNPs, "-", "Minus")

  updatedGWASDataSetFileNameALL_CSV = paste0(filePathOfGWAS, "//", "allGWAS_Data_", fullDiseaseName, "_SNPsWithPvalBelow_", pValName,".csv")
  updatedGWASDataSetFileNameSNPandP_CSV = paste0(filePathOfGWAS, "//", "SNPsAndP_GWAS_", fullDiseaseName, "_SNPsWithPBelow_", pValName,".csv")
  updatedGWASDataSetFileNameALL_RData = paste0(filePathOfGWAS, "//", "allGWAS_Data_", fullDiseaseName, "_SNPsWithPvalBelow_", pValName,".RData")
  updatedGWASDataSetFileNameSNPandP_RData = paste0(filePathOfGWAS, "//", "SNPsAndP_GWAS_", fullDiseaseName, "_SNPsWithPBelow_", pValName,".RData")

dir.create(filePathOfMotifBreakR_Step1)
dir.create(filePathOfMotifBreakR_Step2)
dir.create(filePathOfMotifBreakR_Final)


snpDFPath = paste0(filePathOfMotifBreakR_Final, "motifbreakr_FinalResultsForSNPsWithP_", pValName,".csv")
snpDFPathRData = paste0(filePathOfMotifBreakR_Final, "motifbreakr_FinalResultsForSNPsWithP_", pValName,".RData")

  snpsOutputPathRDataObjects = paste0(snpsOutputPath, "//", "SNPs_Analysis_RDataObjectsFilePaths.RData")

  save(fullDiseaseName, mainOutputFolder, snpsOutputPath, filePathOfGWAS, filePathOfMotifBreakR_Step1,
       filePathOfMotifBreakR_Step2, filePathOfMotifBreakR_Final, pValName,
       updatedGWASDataSetFileNameALL_CSV, updatedGWASDataSetFileNameSNPandP_CSV,
       updatedGWASDataSetFileNameALL_RData, updatedGWASDataSetFileNameSNPandP_RData,
	   snpDFPath, snpDFPathRData,
       file = snpsOutputPathRDataObjects)

  print(paste0(":) Please note that the SNPs and MotifbreakR file paths info (file path names) is stored here: ", snpsOutputPathRDataObjects))
  return(snpsOutputPathRDataObjects)

}



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


dir.create(chromatinInteraction_Step2_FilePath)
dir.create(chromatinInteraction_Step2_FilePath_InitialRDataObjects)
dir.create(chromatinInteraction_Step2_FilePath_CleanedRDataObjects)
dir.create(chromatinInteraction_Step2_FilePath_CleanedCSVObjects)

chromatinInteractionFilesRData = paste0(enhancerAndPromoterInteraction, folderName,"_chromatinInteractionFilePathNames.RData") 

if (haveEnhancersAndPromotersHICData){
  chromatinInteraction_Step1_FilePath = paste0(enhancerAndPromoterInteraction, "Chromatin_scgrnom_Step1//")
  dir.create(chromatinInteraction_Step1_FilePath)
  
  save(fullDiseaseName, mainOutputFolder, enhancerAndPromoterInteraction, folderName,
chromatinInteraction_Step2_FilePath, 
chromatinInteraction_Step2_FilePath_InitialRDataObjects, 
chromatinInteraction_Step2_FilePath_CleanedRDataObjects,
chromatinInteraction_Step2_FilePath_CleanedCSVObjects,
chromatinInteraction_Step1_FilePath,
file = chromatinInteractionFilesRData)
} else { #if (haveEnhancersInteractionData}{
  chromatinInteraction_Step1ADSNPheno_FilePath = paste0(enhancerAndPromoterInteraction, "Chromatin_adsnpheno_Step1//")
  dir.create(chromatinInteraction_Step1ADSNPheno_FilePath)
   save(fullDiseaseName, mainOutputFolder, enhancerAndPromoterInteraction, folderName,
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


dir.create(finalFullGeneRegulatoryNetworkOutputPath)
geneRegulatoryPathwayFileNamesRData = paste0(finalFullGeneRegulatoryNetworkOutputPath, "finalFullGRNRDataObjectsForFilePaths.RData")

save(finalFullNetworkChromatinAndGeneExpressCSV,
fullNetALLOrganizedCSV_Done,
fullNetALLOrganizedCSV,
finalFullGeneRegulatoryNetworkOutputPath,
trenaOutpathALLCSV,
trenaOutpathALLRData,
filePathOfTrena,
trenaOutpathALLRData_Final,
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
filePathOfAll4GRNs,
grnResultsCombinedFileNameCSV,
grnResultsComboForPythonInputFileNameCSV,
grnResultsCombinedFileNameRData,
filePathOfGenie3,
filePathOfRTN, file = geneRegulatoryPathwayFileNamesRData)

print(paste0(":) Please note that the information on the Gene Expression Regulatory Network File Path Names is stored in this RData Object: geneRegulatoryPathwayFileNamesRData = ", geneRegulatoryPathwayFileNamesRData))
return(geneRegulatoryPathwayFileNamesRData)

}