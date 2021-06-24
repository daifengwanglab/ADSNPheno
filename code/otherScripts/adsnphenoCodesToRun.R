
# source("F:\\organizedAlzheimers\\saniyaGithubCodeJune2021\\code\\packagesNeeded.R")
# load("F:\\organizedAlzheimers\\saniyaGithubCodeJune2021\\data\\dataDefaults\\generalDataForPipeline.RData")
# source("F:\\organizedAlzheimers\\saniyaGithubCodeJune2021\\code\\functionsNeeded.R")
# source("F:\\organizedAlzheimers\\saniyaGithubCodeJune2021\\code\\otherDefaultsNeeded.R")
#
#
# source("F:\\organizedAlzheimers\\saniyaGithubCodeJune2021\\code\\ltlCode\\inputsForLTLDemo.R")
#
# source("F:\\organizedAlzheimers\\saniyaGithubCodeJune2021\\code\\filePathDerivations.R")


# :) Please note that this R file by Saniya has all the parameters we need


#source("F:\\organizedAlzheimers\\ADSNPheno\\code\\InformationFromUserInput.R")

print(paste0(":) please note body region = ", bodyRegion))

# Step 2:  How to construct co-expression network and modules (WGCNA and K-means too).  Define parameters

newDate = str_replace_all(Sys.Date(), "-", "_")

numOfRoundingDigits = numberOfRoundingDigits



#####################################
#  Gene Co-Expression Analysis with WGCNA


############# WGCNA CLUSTERING
geneAndEntrezIDMappingDF = geneExpressionDataSetToEntrezIDMapping(inputGeneExpressionDataFilePath,
                                                                  geneToEntrezIDPathForGeneExpressionDataCSV,
                                                                  geneToEntrezIDPathForGeneExpressionDataRData)

datExpr = creatingDatExprFromGeneExpressionDataset(inputGeneExpressionDataFilePath,
                                                   geneToEntrezIDPathForGeneExpressionDataCSV,
                                                   geneToEntrezIDPathForGeneExpressionDataRData,
                                                   geneExpressionOutputNameCSV_withDate,
                                                   geneExpressionOutputNameCSV,
                                                   originalGeneFilteredCSV,
                                                   originalGeneFilteredRData,
                                                   geneExpressionOutputNameRData,
                                                   geneExpressionOutputNameRData_withDate,
                                                   updatedOutputAddOnAndDatExpr_withDate,
                                                   updatedOutputAddOnAndDatExpr)
dim(datExpr)


# please get the traits file:
datTraits = pleaseGetDatTraitsPhenotypesFile(geneExpressionOutputNameRData, phenotypesFilePath, datTraitsRDataFilePath)
numSamples = dim(datExpr)[1]
numGenes = dim(datExpr)[2]

powerEstimate = selectingWGCNAPower(datExpr = datExpr, maxNumPowers = maxNumPowers,
                                    netType = netType, ourPower = ourPower, useRecommendedPower = useRecommendedPower,
                                    log2transformInputData = log2transformInputData, scaleInputData = scaleInputData)


saveWGCNAPowerFilePathsInformation(powerRelatedFilePathsRData, powerEstimate,
                                   folderName, performAdditionalKMeansStep,
                                   tomFilePath, fullDiseaseName,
                                   bodyRegion, outputAddOn, allObjectsFilePath,
                                   genePhenoAssociationFilePath)

load(powerRelatedFilePathsRData) # please note that this stores the file path information related to the WGCNA power used


############ Please define additional variables based on the power estiamte:
# FILE PATH VARIABLES:

filePathInfoDF = updatingFilePathInfoDF(filePathInfoDF,
                                   numSamples, numGenes, powerEstimate,
                                   MEsFilePath, moduleAssignmentsFilePath,
                                   geneModuleMembershipFilePath,
                                   filePathInfoDF_FileNameRData,
                                   filePathInfoDF_FileNameCSV)


phenotypeAssociationFilePath
modulePhenoAssociationFilePath
genePhenoAssociationFilePath
geneModMembershipOutputNameCSV
mmpOutputNameCSV
modTraitCorrelationOutputNameCSV
modTraitCorrPValOutputNameCSV
modTraitCorrAndCorrPValueOutputNameRData
signifModTraitCorrelationOutputNameCSV
updatedModuleTraitPhenotypeDFFilePath
wgcnaWithKmeansAllObjectsExceptTOM
wgcnaWithKmeansAllObjectsIncludingTOM
wgcnaWithKmeansAllObjectsExceptTOMWithDate
wgcnaWithKmeansAllObjectsIncludingTOMWithDate
wgcnaWithKmeansGeneTraitSignificancePath
wgcnaWithKmeansGeneTraitSigPValuePath
wgcnaWithKmeansGeneModuleMembershipPath
wgcnaWithKmeansGeneModMemberPValuePath
wgcnaWithKmeansGeneTraitsEdgeTablePathCSV
wgcnaWithKmeansGeneTraitsEdgeTablePathRData
signifGeneTraitCorrelationOutputNameCSV
signifPositiveGeneTraitCorrelationOutputNameCSV
updatedGeneTraitPhenotypeDFFilePath
wgcnaPowerOutputNameRData
moduleCountsOutputName
wgcnaOutputName

#################################################################
print("############")
print(paste0(":) Please note that the powerEstimate we use will be: ", powerEstimate))
# a check here to see if we wanted to use recommendations or our own power
adjacencyPower = adjacency(datExpr, power = powerEstimate, type = netType)
# "Please note that the Soft-Threshold Power Selected is: 11"
# --> NEW "Please note that the Soft-Threshold Power Selected is: 10"
# Please Turn adjacency into topological overlap :) matrix (TOM)
TOMPower = TOMsimilarity(adjacencyPower);

rm(adjacencyPower)
dissTOMPower = 1 - TOMPower

save(dissTOMPower,  file = dissTOMPowerOutputNameRData)




#write.csv(TOMPower, file = tomPowerOutputName)

geneTreePower = hclust(as.dist(dissTOMPower), method = "average");

# Module identification using dynamic tree cut:
# We like large modules, so we set the minimum module size relatively high:
dynamicModsPower = cutreeDynamic(dendro = geneTreePower, distM = dissTOMPower,
                                 deepSplit = 2, pamRespectsDendro = FALSE,
                                 minClusterSize = minModuleSize);

dynamicColorsPower = labels2colors(dynamicModsPower)



#####################
MEList = moduleEigengenes(datExpr, colors = dynamicColorsPower)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);

save(dissTOMPower, MEList, MEs, MEDiss, geneTreePower, powerEstimate,  dynamicModsPower, dynamicColorsPower,
     file = wgcnaPowerOutputNameRData)

write.csv(data.frame(table(dynamicColorsPower)), moduleCountsOutputName) #    "D:\\initialDLPFC_modules.csv")

#####################
# Cluster module eigengenes
diseaseName = disease #'Alzheimers'
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = paste("Initial Clustering of Module Eigengenes (", diseaseName, "; ", bodyRegion,")"),
     xlab = paste0(netType, " Network with Power: ", powerEstimate), sub = "")

initialModuleCounts = data.frame(table(dynamicColorsPower))
print(initialModuleCounts)



########################################################
if (treeCutHeightForMergingCloseModules > 0){
  MEDissThres = treeCutHeightForMergingCloseModules
  # Plot the cut line into the dendrogram
  abline(h=MEDissThres, col = "red")
  # Call an automatic merging function
  merge = mergeCloseModules(datExpr, dynamicColorsPower, cutHeight = MEDissThres, verbose = 3)
  # The merged module colors
  dynamicColorsPower = merge$colors;
  # Eigengenes of the new merged modules:
  MEs = merge$newMEs;

  print(paste0("Please note that since treeCutHeightForMergingCloseModules is: ", treeCutHeightForMergingCloseModules, " we therefore merged close modules togehter."))
} else {
  print(paste0(":) Please note that we did not merge any close modules since treeCutHeightForMergingCloseModules is 0"))
}

############################################################


dynamicColorsKmeans = performingKMeansAfterWGCNA(dynamicColorsPower, performAdditionalKMeansStep,
                                       datExpr, meg, nIterations)

data.frame(table(dynamicColorsKmeans))
dim(datExpr)
MEList = moduleEigengenes(datExpr, colors = dynamicColorsKmeans)
MEs = MEList$eigengenes

write.csv(MEs, MEsWithKmeansOutputName)
gene_names = colnames(datExpr)
clusterInfo = data.frame(module = dynamicColorsKmeans)
rownames(clusterInfo)= gene_names
dim(clusterInfo)
clusterInfo = merge(clusterInfo, geneAndEntrezIDMappingDF, by = 0)
gene_names = clusterInfo[,1] # newer Chanu :)
rownames(clusterInfo)= clusterInfo[,1]#gene_names
head(clusterInfo)
clusterInfo = clusterInfo[,-1]
colnames(clusterInfo) = c("module", "entrezID")
head(clusterInfo)
clusterInfo$bodyRegion = rep(bodyRegion, nrow(clusterInfo))
clusterInfo$disease = rep(fullDiseaseName, nrow(clusterInfo))

write.csv(clusterInfo, geneAssignmentsWithKmeansOutputNameCSV)
save(clusterInfo, file = geneAssignmentsWithKmeansOutputNameRData)
data.frame(table(clusterInfo))

write.csv(data.frame(table(clusterInfo$module)), moduleCountsWithKmeansOutputNameCSV)


###############################################
# Finding Gene Module Membership:

dim(datExpr)
nSamples = nrow(datExpr)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
geneModMembershipOutputNameCSV = paste0(wgcnaAndKMeansOutputPath, "//WGCNA_withKMeans_",netType, "_pow",powerEstimate, "_geneModuleMembershipCorr", outputAddOn, ".csv")
mmpOutputNameCSV = paste0(wgcnaAndKMeansOutputPath, "//WGCNA_withKMeans_",netType, "_pow",powerEstimate, "_MMPvalues", outputAddOn, ".csv")

geneModuleMembership = merge(geneAndEntrezIDMappingDF, geneModuleMembership, by = 0)
rownames(geneModuleMembership)=gene_names
geneModuleMembership = geneModuleMembership[,-1]
colnames(geneModuleMembership)[1] = "entrezID"
MMPvalue = merge(geneAndEntrezIDMappingDF, MMPvalue, by = 0)
colnames(MMPvalue)[1] = "entrezID"
rownames(MMPvalue)=gene_names
MMPvalue = MMPvalue[,-1]
head(MMPvalue)
head(clusterInfo)


write.csv(geneModuleMembership, geneModMembershipOutputNameCSV)
write.csv(MMPvalue, mmpOutputNameCSV)


########################################################################
# Finding Hub Gene in Each Gene Module:
uniqueModules = unique(dynamicColorsKmeans)
print(paste0(":) Please note that there are: ", length(uniqueModules), " total unique gene co-expression modules here!"))
hubs  = chooseTopHubInEachModule(datExpr, dynamicColorsKmeans)
head(dynamicColorsKmeans)
hubsDF = data.frame(hubs)
hubsDF$power = rep(powerEstimate, nrow(hubsDF))
hubsDF$region = rep(bodyRegion, nrow(hubsDF))
head(hubsDF)
dim(hubsDF)
write.csv(hubsDF, hubsCSV)

##########################################
# Finding the Co-Expression Network for Each Gene Module
allModulesCoExpressionNetworkModularTOMResultsDF = data.frame() # only for modular edges
allNetworkModulesTomSummaryStatisticsDF = data.frame()
uniqueModules = unique(dynamicColorsKmeans)
TOMPower = 1 - dissTOMPower
for (moduleNum in 1:length(uniqueModules)){
  moduleToUse = uniqueModules[moduleNum]
  print(paste0(":) please note: moduleNum = ", moduleNum, " of ", length(uniqueModules),
               " and current gene module for ", bodyRegion, " is: ", moduleToUse))
  moduleCoExpressionResultsList1 = pleaseGetCoexpressionNetworkForAWGCNAModuleAfterKmeans(module = moduleToUse)
  moduleInteractionsDF = data.frame(moduleCoExpressionResultsList1$moduleInteractionsDF)
  summaryStatisticsOnPercentilesDF = data.frame(moduleCoExpressionResultsList1$summaryStatisticsOnPercentilesDF)

  if (moduleNum == 1){
    allNetworkModulesCoExpressionNetworkDF = moduleInteractionsDF
    allNetworkModulesTomSummaryStatisticsDF = summaryStatisticsOnPercentilesDF
  } else {
    allNetworkModulesCoExpressionNetworkDF = rbind(allNetworkModulesCoExpressionNetworkDF,
                                                   moduleInteractionsDF)
    allNetworkModulesTomSummaryStatisticsDF = rbind(allNetworkModulesTomSummaryStatisticsDF,
                                                    summaryStatisticsOnPercentilesDF)

  }
  print("*********************************************************** \n")
}
print(paste0(":) please note we are writing out the Co-Expression Network for ALL ", length(uniqueModules), " here:", allNetworkModulesCombinedTomSummaryStatisticsFileName))
write.csv(allNetworkModulesTomSummaryStatisticsDF, allNetworkModulesCombinedTomSummaryStatisticsFileName)


##################################################################################################################################
# Finding Significant Gene-Phenotype Relationships
source(functionsNeededSourceCode)
saveWGCNAPowerFilePathsInformation(powerRelatedFilePathsRData, powerEstimate,
                                   folderName, performAdditionalKMeansStep,
                                   tomFilePath, fullDiseaseName,
                                   bodyRegion, outputAddOn, allObjectsFilePath,
                                   genePhenoAssociationFilePath)

load(powerRelatedFilePathsRData)
genePhenoList = findingGenePhenotypeCorrelationsAndPValues(datTraits, datExpr,
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
                                                              dynamicColorsKmeans)


genePhenoList
updatedGeneTraitPhenotypeDF = genePhenoList$updatedGeneTraitPhenotypeDF
geneTraitSignificanceDF = genePhenoList$geneTraitSignificanceDF
GSPvalueDF = genePhenoList$GSPvalueDF
significantGeneTraitsDF = genePhenoList$significantGeneTraitsDF

###########################################################
###############################################
#### Module Trait Correlations:

# Finding Module-Trait Correlations:

modulePhenoList = findingModulePhenotypeCorrelationsAndPValues(MEs, datTraits,
                                                               modTraitCorrelationOutputNameCSV,
                                                               datExpr, modTraitCorrPValOutputNameCSV,
                                                                modTraitCorrAndCorrPValueOutputNameRData,
                                                                modulePhenoAssociationFilePath,
                                                               performAdditionalKMeansStep,
                                                                pValueCutOffForSignificance, powerEstimate,
                                                                signifModTraitCorrelationOutputNameCSV,
                                                                positiveCorrelationsOnly,
                                                                posCorsOnlyModuleTraitsOutputNameCSV,
                                                               numberOfRoundingDigits, dynamicColorsKmeans, fullDiseaseName,                                                                           tissueName,                                                                           bodyRegion,                                                                           dataScaling,
                                                               updatedModuleTraitPhenotypeDFFilePath)


significantModTraits = modulePhenoList$significantModTraits
updatedModuleTraitPhenotypeDF = modulePhenoList$updatedModuleTraitPhenotypeDF
moduleTraitCor = modulePhenoList$moduleTraitCor
moduleTraitPvalue = modulePhenoList$moduleTraitPvalue
info = modulePhenoList$info
##############################################################################
#  Saving all of the data objects from WGCNA
keyModules = unique(updatedModuleTraitPhenotypeDF$moduleName)

# getting metrics on the gene, assigned module, significant phenotypes for gene, and significant phenotypes of gene's assigned module:
source_python(pathForPythonCode)
geneModuleAssignmentsPhenosDF = pleaseGetGenesModulesPhenotypesAssignmentsOrganized(powerEstimate, geneAssignmentsWithKmeansOutputNameCSV, updatedModuleTraitPhenotypeDFFilePath, updatedGeneTraitPhenotypeDFFilePath, geneAssignmentsModAssignmentsPhenosWithKmeansOutputNameCSV)
geneModuleAssignmentsPhenosDF
View(geneModuleAssignmentsPhenosDF)
currentDate = as.character(now())
save(significantModTraits,
     updatedModuleTraitPhenotypeDF,
     updatedGeneTraitPhenotypeDF,
     significantGeneTraitsDF,
     bodyRegion,
     currentDate,
     folderName,
     disease,
     fullDiseaseName,
     tissueName,
     dataScaling,
     datExpr,
     geneModuleAssignmentsPhenosDF,
     geneAndEntrezIDMappingDF,
     tfsDF,
     datTraits,
     moduleTraitCor,
     moduleTraitPvalue,
     geneModuleMembership,
     MMPvalue,
     clusterInfo,
     dynamicColorsKmeans,
     MEList, MEs, hubs,
     allNetworkModulesTomSummaryStatisticsDF,
     powerEstimate, info,
     geneTraitSignificanceDF,
     GSPvalueDF, keyModules,
     file = wgcnaWithKmeansAllObjectsExceptTOM)


print(paste0(":) please note that there are ", length(keyModules), " gene modules that are associated with at least 1 of the phenotypes"))
clusterInfoTemp = clusterInfo # please hold onto the original clusterInfo for all modules
clusterInfo = clusterInfo[which(clusterInfo[,grep("module", colnames(clusterInfo))] %in% keyModules),]
write.csv(clusterInfo, keyGeneAssignmentsWithKmeansOutputNameCSV) #paste0(outputPath, "wgcna_1029_with_kmeans_pow", powerEstimate, "_moduleAssignments__newRNASeqAlzh.csv"))
save(clusterInfo, file = keyGeneAssignmentsWithKmeansOutputNameRData) #"D:\\DLPFC_updatedGeneClusterInfo.RData")
clusterInfo = clusterInfoTemp # please reset clusterInfo back to all the modules after we saved keyModules
rm(clusterInfoTemp) # please remove temp variable

#################################################################################

##########################################################

# Step 3:  Annotate modular functions by enrichment analysis
load(powerRelatedFilePathsRData)
load(wgcnaWithKmeansAllObjectsExceptTOM)

region = folderName
regionName = folderName

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
colnames(clusterInfo) #  "??.."      "entrezID" "module24" "module26" "module44" "module30" "module" "module37"
dim(clusterInfo)
entrezIDCol = grep("entrez", colnames(clusterInfo))
print(paste0(":) Please note the entrezIDCol is: ", entrezIDCol))
modColNum = grep("module", colnames(clusterInfo)) # module column number

print(paste0(":) please note that there are ", length(keyModules), " gene modules that are associated with at least 1 of the phenotypes"))


if (performMESHEnrichment){
  print(paste0(":) Please note that since performMESHEnrichment is", performMESHEnrichment, " we are performing Medical Subject Enrichment Analysis! :)"))
  #meshEnrichmentsFilePath = paste0(moduleEnrichmentsFilePath, "//MeSH")
  #dir.create(meshEnrichmentsFilePath, showWarnings = FALSE)
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
  performMolSigDBEnrichment = TRUE

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



###############################################################################################################################################################################
######## SNPs breaking TF Binding Sites:
motifBreakRSNPsDF = findingSNPsBreakingTFBindingSitesFromGWASDataFunction(gwasDataFile,
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
                                                                          snpDFPathRData)


#####################################################################################
# 5.  Please predict gene regulatory networks for genes and modules
# Step 5:
################################################################################
# Chromatin Interaction Part:
# Please note that since the chromosomes are usually smaller in size the higher the number of the chromosome,
# we start with the later chromosomes to ensure those work :)
# # Step2: please get TFs for each promoter and enhancer
# # https://github.com/daifengwanglab/scGRN

##########################################################################################
# PLEASE NOTE THIS DEMO FOR 1 CHROMOSOME:
# user can have the option to create an interactionsDF or to run scgrnom to get it
# caseInfo

scgrnGetTFsForEachChromosomesChromatinNetworkVec = pleaseBuildChromatinGeneRegulatoryNetwork(numOfCoresToUseForScgrnom,
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
                                                                                             enhancerAndPromoterInteraction)


#########################################################################
## GENE EXPRESSION GENE REGULATORY NETWORKS:
##################################################################
load(powerRelatedFilePathsRData)
load(wgcnaWithKmeansAllObjectsExceptTOM)


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
combinedGeneExpressGRNDF


parentName = paste0(transcriptionalRegulatoryNetworkOutputPath, "//")# #'F://organizedAlzheimers//FinalGeneRegulatoryNetworks//Alzheimers_Brain_Hippocampus//GeneRegulatNetwork//'



scgrnGetTFsForEachChromosomesChromatinNetworkVec = paste0(chromatinInteraction_Step2_FilePath_CleanedCSVObjects, list.files(chromatinInteraction_Step2_FilePath_CleanedCSVObjects, pattern = "scgrnomGetTFs_"))
scgrnGetTFsForEachChromosomesChromatinNetworkVec
SNPOutputNameStem
if (haveChromatinInteractionRegulatoryNetwork){
  hasFinalFullNetworkForChromosome = "TRUE"
} else {
  hasFinalFullNetworkForChromosome = "FALSE"
}


# chromatinRegNetFileNameList = paste0(chromatinInteraction_Step2_FilePath_CleanedCSVObjects,
#                                      scgrnGetTFsForEachChromosomesChromatinNetworkVec)
chromatinRegNetFileNameList = scgrnGetTFsForEachChromosomesChromatinNetworkVec
source_python(pathForPythonCode)

motifbreakRFinalResultsFilePath = snpDFPath
#combinedSubsForAllChromsDF = gettingSNPSubNetwork(chromatinRegNetFileNameList, grnResultsCombinedFileNameCSV, motifbreakRFinalResultsFilePath, parentName, disease, bodyRegion, dataScalingOutputMini, geneToEntrezIDMappingPath, minNumSourcesGeneGRN, enhancer_buffer_kbp, promoter_buffer_kbp, includeTheIndividualTFsInGroup_Python, SNPOutputNameStem, geneAssignmentsModAssignmentsPhenosWithKmeansOutputNameCSV,powerEstimate,  progressBarPythonIterations = 1000, filePathForFullNetwork = "", hasFinalFullNetworkForChromosome = "FALSE")
#View(combinedSubsForAllChromsDF)
#combinedSubsForAllChromsDF






pathForPythonCode
motifbreakRFinalResultsFilePath = snpDFPath

source_python(pathForPythonCode)
snpDFPath
combinedSubsForAllChromsDF = gettingSNPSubNetwork(chromatinRegNetFileNameList,
                                                  grnResultsCombinedFileNameCSV,
                                                  snpDFPath,
                                                  parentName, disease, bodyRegion,
                                                  dataScalingOutputMini,
                                                  geneToEntrezIDMappingPath,
                                                  minNumSourcesGeneGRN,
                                                  enhancer_buffer_kbp,
                                                  promoter_buffer_kbp,
                                                  snpDFPath,
                                                  SNPOutputNameStem,
                                                  geneAssignmentsModAssignmentsPhenosWithKmeansOutputNameCSV,
                                                  powerEstimate,
                                                  progressBarPythonIterations, filePathForFullNetwork, hasFinalFullNetworkForChromosome)#, filePathForFullNetwork = NULL, hasFinalFullNetworkForChromosome = "False")



finalGRN = pleaseGetOrganizedFinalGeneRegulatoryNetwork(geneAssignmentsModAssignmentsPhenosWithKmeansOutputNameCSV,
                                                        finalFullNetworkChromatinAndGeneExpressCSV,
                                                        parentName, disease, bodyRegion,
                                                        dataScalingOutputMini, minNumSourcesGeneGRN,
                                                        chromatinRegNetFileNameList, grnResultsCombinedFileNameCSV,
                                                        geneToEntrezIDMappingPath,
                                                        enhancer_buffer_kbp, promoter_buffer_kbp,
                                                        powerEstimate,
                                                        progressBarPythonIterations,
                                                        hasFinalFullNetworkForChromosome,
                                                        filePathForFullNetwork)  #

finalGRN
combinedSubsForAllChromsDF




#############################################################################
print(paste("chromatinRegNetFileNameList:",chromatinRegNetFileNameList))
print()
print(paste("grnResultsCombinedFileNameCSV:",grnResultsCombinedFileNameCSV))
print()
print(paste("motifbreakRFinalResultsFilePath:",motifbreakRFinalResultsFilePath))
print()
print(paste("parentName:",parentName))
print()
print(paste("disease:",disease))
print()
print(paste("bodyRegion:",bodyRegion))
print()
print(paste("dataScalingOutputMini:",dataScalingOutputMini))
print()
print(paste("geneToEntrezIDMappingPath:",geneToEntrezIDMappingPath))
print()
print(paste("minNumSourcesGeneGRN:",minNumSourcesGeneGRN))
print()
print(paste("enhancer_buffer_kbp:",enhancer_buffer_kbp))
print()
print(paste("promoter_buffer_kbp:",promoter_buffer_kbp))
print()
print(paste("includeTheIndividualTFsInGroup_Python:",includeTheIndividualTFsInGroup_Python))
print()
print(paste("SNPOutputNameStem:",SNPOutputNameStem))
print()
print(paste("filePathForFullNetwork:", filePathForFullNetwork))
print()
print(paste("hasFinalFullNetworkForChromosome", hasFinalFullNetworkForChromosome))
print()


print(paste0(":) YAY!! WE ARE ALL DONE!! :)"))

