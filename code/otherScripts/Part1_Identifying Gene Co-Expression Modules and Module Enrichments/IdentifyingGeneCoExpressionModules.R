#########################################################################################################################################
### INFORMATION:  Please note that the code below was written by Saniya Khullar and Daifeng Wang, Ph.D.
#########################################################################################################################################
### USER INPUTS:

adsnphenoDirectory = "F://organizedAlzheimers//ADSNPheno//"  # please supply full directory path
setwd(adsnphenoDirectory) # please set this folder as working directory for code and other resources
codeFolder = paste0(adsnphenoDirectory, "code//") # the main folder where code would be within ADSNPheno

dataFolder = paste0(adsnphenoDirectory, "data//") # the main folder where data would be within ADSNPheno
# We use a combined list by Lambert and Jaspar for our TFs
transcriptionFactorsFilePath = "LambertAndJasparTFs.csv"

tfsUsed = "LambertAndJaspar"

# please note the location of the Gene Name to Entrez ID mapping
geneToEntrezIDMappingPath = "F://organizedAlzheimers//saniyaGithubCodeJune2021//data//dataDefaults//entrezMappingInfo//infoDFMappingsGeneSymbolAndID.csv"

# our background gene regulatory data came from TRRUST2, so we input:
# please note that this data source should have columns like this: "TF", "RegulatedGene", "Source", "Info", "CombinedName"
additionalGeneRegulatoryNetworkFilePath = "Trrust2.csv"

tfSource = "TRRUST2"

numberOfRoundingDigits = 3
includeTimeStamp = "TRUE"

outputPathNameADSNPhenoOutputs = paste0("F://organizedAlzheimers//ADSNPhenoOutputs_NEW//")



#getwd() # please note this is our current working directory
#setwd(codeFolder) # please set this folder as working directory


disease = "Alzheimers"
diseasePetName = "AD" # or diseaseNickName

tissueName = "Brain"
bodyRegion = "MiniFinalDemo_LTL1" 
inputGeneExpressionDataFilePath = paste0(dataFolder, "//originalGeneExpressionDataSets//LateralTemporalLobeRegionGeneExpressionData_miniDemo_200genes.csv") #"F://organizedAlzheimers//Setup//InputDataFiles//originalGeneExpressionDataSets//LateralTemporalLobeRegionGeneExpressionData_miniDemo.csv" #"F://organizedAlzheimers//Setup//InputDataFiles//originalGeneExpressionDataSets//originalLateralTemporalLobeRegionGeneExpressionData.csv"

useTFsFoundByADSNPheno = TRUE # please note that by default we use  combined list of Lambert and Jaspar TFs.
filePathOfTFsToUse = NULL # please note that if this is NULL, we will use combined list of Lambert and Jaspar TFs.  Else, please specify csv filepath of TFs, where the TFs are in 1 column with the title called "tfName", with names that are found exactly in gene expression dataset

log2transformInputData = "TRUE" # (should we apply a log2(x+1) transform on data
scaleInputData = "FALSE"  # should we apply a scale() on data
phenotypesFilePath = paste0(dataFolder, "phenotypes//AlzheimersLateralTemporalLobePhenotypesUpdated.csv") 
########################################################################################################################
### Weighted Gene Co-Expression Network Analysis (WGCNA) with K-Means:
# Step 2:  How to construct co-expression network and modules (WGCNA and K-means too).  Define parameters

saveOutputs = "TRUE"
netType = "signed"  # # network type:but please check WGCNA for additional parameters.  We retained the defaults

# Soft-Thresholding Power to be used:
maxNumPowers = 50 # when selecting a Soft-Thresholding Power, what is the maximum Power, beta, that you want to consider?

ourPower = NULL # if we have a recommended power in mind please put that there (this can come from analyzing the Soft-Thresholding Analysis results)
useRecommendedPower = "TRUE" #FALSE" # otherwise, please use the default that WGCNA selects


# module trait and gene-trait associations:
positiveCorrelationsOnly = TRUE #   for module-trait associations and gene-trait associations
pValueCutOffForSignificance = 0.05 # please note that associations with a p-value < pValueCutOffForSignificance will be considered statistically significant

# co-expression network part:
percentilesForCoexpressionNetworkVec = c(0, 0.25, 0.50, 0.75, 0.85, 0.90, 0.95, 1)
minTOMThresholdValue = 0.0001

# when WGCNA is performing initial step of creating modules from the TOM, what is the minimum module size?
minModuleSize = 30;
# tree cut height ??  we kept it at 0, but individuals could adjust as needed to help merge modules?
treeCutHeightForMergingCloseModules = 0  # we usually do not merge close modules

# K-Means step performed after WGCNA
performAdditionalKMeansStep = TRUE # or can be FALSE if you want original WGCNA results

# please note source of this code:
# https://github.com/juanbot/km2gcn
min.genes.for.grey=0 # originally, grey module corresponds a group of unclustered genes.  How many minimum genes do we want in grey module after k-means?
# stopping conditions
nIterations = 201 # what is the maximum # of iterations of k-means that we are willing to perform?
meg = 0 # Minimum number of genes exchanged: if before the number of iterations, the number of genes exchanged in the current iteration is less than this parameter, k-means algorithm stops.




########################################################################################################################################
### PLEASE RUN: 

########################################################### running other things:
generalDataForPipeline = paste0(dataFolder, "generalDataForPipeline.RData")
load(generalDataForPipeline)
pathForPythonCode = paste0(codeFolder, "setupScripts//BuildingFullAndSubNetworksPythonCode.py")
#filePathDerivationsSourceCode = paste0(codeFolder, "setupScripts//filePathDerivations.R")
packagesNeededSourceCode = paste0(codeFolder, "setupScripts//packagesNeeded.R")
functionsNeededSourceCode = paste0(codeFolder, "setupScripts//functionsNeeded.R")
#userInputsSourceCode = paste0(codeFolder, "pythonCode//BuildingFullAndSubNetworksPythonCode.py")
#generalDataForPipeline = paste0(dataFolder, "generalDataForPipeline.RData")

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



######################################################################################


initialRDataFilePath_FilePathDerivationsAndInfoObjectsForWGCNA = pleaseGetInitialFilePathDerivationsAndInfoObjectsForWGCNA(outputPathNameADSNPhenoOutputs,
disease, tissueName, bodyRegion, tfsUsed, performAdditionalKMeansStep, phenotypesFilePath, log2transformInputData, scaleInputData, includeTimeStamp, numberOfRoundingDigits, netType)

if (performAdditionalKMeansStep){
print("hi yes")
  } else {
    print("nooo")
  #allNetworkModulesCombinedTomSummaryStatisticsFileName = paste0(wgcnaAndKMeansOutputPath, "//", bodyRegion, "_wgcnaOnly_Power", powerEstimate, "finalCombo_CoExpressNetSummaryStats_forAllModules",outputAddOn,".csv")
}
if (performAdditionalKMeansStep){
  allNetworkModulesCombinedTomSummaryStatisticsFileName = paste0("//", bodyRegion, "_wgcnaWith_kmeans_Power", powerEstimate, "finalCombo_CoExpressNetSummaryStats_forAllModules",outputAddOn,".csv")
} else {
  allNetworkModulesCombinedTomSummaryStatisticsFileName = paste0("//", bodyRegion, "_wgcnaOnly_Power", powerEstimate, "finalCombo_CoExpressNetSummaryStats_forAllModules",outputAddOn,".csv")
}


load(initialRDataFilePath_FilePathDerivationsAndInfoObjectsForWGCNA)
#load(wgcnaWithKmeansAllObjectsExceptTOM_SimpleFileName) # please note that after we load initialRDataFilePath_FilePathDerivationsAndInfoObjectsForWGCNA, we have names for many file paths and can load this to get WGCNA results :) YAY!


print(paste0(":) please note body region = ", bodyRegion))

# Step 2:  How to construct co-expression network and modules (WGCNA and K-means too).  Define parameters

newDate = str_replace_all(Sys.Date(), "-", "_")

#numOfRoundingDigits = numberOfRoundingDigits



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
# powerRelatedFilePathsRData = saveWGCNAPowerFilePathsInformation(powerRelatedFilePathsRData, powerEstimate,
#                                    folderName, performAdditionalKMeansStep,
#                                    tomFilePath, fullDiseaseName,
#                                    bodyRegion, outputAddOn, allObjectsFilePath,
#                                    genePhenoAssociationFilePath)

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
     powerRelatedFilePathsRData,
     file = wgcnaWithKmeansAllObjectsExceptTOM)

# we also want to save these same objects in this field so we can call them from other methods without needing info on the power :)
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
     powerRelatedFilePathsRData,
     file = wgcnaWithKmeansAllObjectsExceptTOM_SimpleFileName)


#wgcnaWithKmeansAllObjectsExceptTOM_SimpleFileName


print(paste0(":) please note that there are ", length(keyModules), " gene modules that are associated with at least 1 of the phenotypes"))
clusterInfoTemp = clusterInfo # please hold onto the original clusterInfo for all modules
clusterInfo = clusterInfo[which(clusterInfo[,grep("module", colnames(clusterInfo))] %in% keyModules),]
write.csv(clusterInfo, keyGeneAssignmentsWithKmeansOutputNameCSV) #paste0(outputPath, "wgcna_1029_with_kmeans_pow", powerEstimate, "_moduleAssignments__newRNASeqAlzh.csv"))
save(clusterInfo, file = keyGeneAssignmentsWithKmeansOutputNameRData) #"D:\\DLPFC_updatedGeneClusterInfo.RData")
clusterInfo = clusterInfoTemp # please reset clusterInfo back to all the modules after we saved keyModules
rm(clusterInfoTemp) # please remove temp variable




