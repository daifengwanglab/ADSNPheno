powerEstimate = 32

# a check here to see if we wanted to use recommendations or our own power
adjacencyPower = adjacency(datExpr, power = powerEstimate, type = netType)
# "Please note that the Soft-Threshold Power Selected is: 11"
# --> NEW "Please note that the Soft-Threshold Power Selected is: 10"
# Please Turn adjacency into topological overlap :) matrix (TOM)
TOMPower = TOMsimilarity(adjacencyPower);
dissTOMPower = 1 - TOMPower


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






wgcnaFileInfoDF = data.frame(c("fileCreated:", as.character(now())))
wgcnaFileInfoDF = rbind(wgcnaFileInfoDF, data.frame(c("bodyRegion:", bodyRegion)))
wgcnaFileInfoDF = rbind(wgcnaFileInfoDF, data.frame(c("tissueName:", tissueName)))
wgcnaFileInfoDF = rbind(wgcnaFileInfoDF, data.frame(c("disease:", disease)))
wgcnaFileInfoDF = rbind(wgcnaFileInfoDF, data.frame(c("fullDiseaseName:", fullDiseaseName)))
wgcnaFileInfoDF = rbind(wgcnaFileInfoDF, data.frame(c("minModuleSize:", minModuleSize)))
wgcnaFileInfoDF = rbind(wgcnaFileInfoDF, data.frame(c("maxNumPowers", maxNumPowers)))
wgcnaFileInfoDF = rbind(wgcnaFileInfoDF, data.frame(c("min.genes.for.grey:", min.genes.for.grey)))
wgcnaFileInfoDF = rbind(wgcnaFileInfoDF, data.frame(c("max # of kMeans Iterations", nIterations)))
wgcnaFileInfoDF = rbind(wgcnaFileInfoDF, data.frame(c("meg", meg)))
wgcnaFileInfoDF = rbind(wgcnaFileInfoDF, data.frame(c("network type (netType)", netType)))
wgcnaFileInfoDF = rbind(wgcnaFileInfoDF, data.frame(c("numRoundingDigits:", numRoundingDigits)))
wgcnaFileInfoDF = rbind(wgcnaFileInfoDF, data.frame(c("pValueCutOffForSignificance:", pValueCutOffForSignificance)))
wgcnaFileInfoDF = rbind(wgcnaFileInfoDF, data.frame(c("use Recommended Power by WGCNA", useRecommendedPower)))
wgcnaFileInfoDF = rbind(wgcnaFileInfoDF, data.frame(c("ourPower:", ourPower)))
wgcnaFileInfoDF = rbind(wgcnaFileInfoDF, data.frame(c("positiveCorrelationsOnly?", positiveCorrelationsOnly)))
wgcnaFileInfoDF = rbind(wgcnaFileInfoDF, data.frame(c("inputGeneExpressionDataFilePath:", inputGeneExpressionDataFilePath)))
wgcnaFileInfoDF = rbind(wgcnaFileInfoDF, data.frame(c("wgcnaAndKMeansOutputPath:", wgcnaAndKMeansOutputPath)))
wgcnaFileInfoDF = rbind(wgcnaFileInfoDF, data.frame(c("geneToEntrezIDMappingPath:", geneToEntrezIDMappingPath)))
wgcnaFileInfoDF = rbind(wgcnaFileInfoDF, data.frame(c("phenotypesFilePath:", phenotypesFilePath)))
wgcnaFileInfoDF = rbind(wgcnaFileInfoDF, data.frame(c("log2transformInputData:", log2transformInputData)))
wgcnaFileInfoDF = rbind(wgcnaFileInfoDF, data.frame(c("scaleInputData:", scaleInputData)))
wgcnaFileInfoDF = rbind(wgcnaFileInfoDF, data.frame(c("dataScaling:", dataScaling)))
wgcnaFileInfoDF = rbind(wgcnaFileInfoDF, data.frame(c("geneToEntrezIDMapping RData File:", geneToEntrezIDPathForGeneExpressionDataRData)))
wgcnaFileInfoDF = rbind(wgcnaFileInfoDF, data.frame(c("geneToEntrezIDMapping CSV File:", geneToEntrezIDPathForGeneExpressionDataCSV)))
wgcnaFileInfoDF = rbind(wgcnaFileInfoDF, data.frame(c("original gene expression RData File (before filtering):", originalGeneFilteredRData)))
wgcnaFileInfoDF = rbind(wgcnaFileInfoDF, data.frame(c("original gene expression CSV File (before filtering):", originalGeneFilteredCSV)))
wgcnaFileInfoDF = rbind(wgcnaFileInfoDF, data.frame(c(paste0("final gene expression CSV File (after filtering) and ", dataScaling, " with date:"), geneExpressionOutputNameCSV_withDate)))
wgcnaFileInfoDF = rbind(wgcnaFileInfoDF, data.frame(c(paste0("final gene expression CSV File (after filtering) and ", dataScaling, " :"), updatedOutputAddOn)))
wgcnaFileInfoDF = rbind(wgcnaFileInfoDF, data.frame(c(paste0("final gene expression RData File (after filtering) and ", dataScaling, ":"), geneExpressionOutputNameRData_withDate)))
wgcnaFileInfoDF = rbind(wgcnaFileInfoDF, data.frame(c(paste0("final DatExpr CSV File (after filtering) and ", dataScaling, " with date:"), updatedOutputAddOnAndDatExpr_withDate)))
wgcnaFileInfoDF = rbind(wgcnaFileInfoDF, data.frame(c(paste0("final DatExpr CSV File (after filtering) and ", dataScaling, " :"), updatedOutputAddOnAndDatExpr)))




sftThresholdingOutputNameCSV = paste0(wgcnaAndKMeansPathsForOutputs, "//softThresholdPowerDF_",netType, "_network_recommendedPower_",fullDiseaseName, "_", bodyRegion, outputAddOn, ".csv")
sftThresholdingOutputNameRData = paste0(wgcnaAndKMeansPathsForOutputs, "//softThresholdPowerDF_",netType, "_network_recommendedPower_",fullDiseaseName, "_", bodyRegion, outputAddOn, ".RData")


dissTOMPowerOutputNameCSV = paste0(wgcnaAndKMeansPathsForOutputs, "//dissTOMPower_",netType, "_power",powerEstimate, "_", fullDiseaseName, "_", bodyRegion, outputAddOn, ".csv")
dissTOMPowerOutputNameRData = paste0(wgcnaAndKMeansPathsForOutputs, "//dissTOMPower_",netType, "_power",powerEstimate, "_", fullDiseaseName, "_", bodyRegion, outputAddOn, ".RData")

tomPowerOutputNameCSV = paste0(wgcnaAndKMeansPathsForOutputs, "//TOMPower_",netType, "_power",powerEstimate, "_", fullDiseaseName, "_", bodyRegion, outputAddOn, ".csv")
tomPowerOutputNameRData = paste0(wgcnaAndKMeansPathsForOutputs, "//TOMPower_",netType, "_power",powerEstimate, "_", fullDiseaseName, "_", bodyRegion, outputAddOn, ".RData")

wgcnaPowerOutputNameRData = paste0(wgcnaAndKMeansPathsForOutputs, "//wgcnaInitialResultsObjects_",netType, "_power",powerEstimate, "_", fullDiseaseName, "_", bodyRegion, outputAddOn, ".RData")

moduleCountsOutputName = paste0(wgcnaAndKMeansPathsForOutputs, "//WGCNAInitialModuleCounts_",netType, "_power",powerEstimate, outputAddOn, ".csv")

wgcnaOutputName = paste0(wgcnaAndKMeansPathsForOutputs, "//WGCNAInitialResultsObjects_",netType, "_power",powerEstimate, outputAddOn, ".RData")



MEsWithKmeansOutputName = paste0(wgcnaAndKMeansPathsForOutputs, "//WGCNA_withKMeans_",netType, "_power",powerEstimate, "_ModuleEigengenesMEs", outputAddOn, ".csv")
geneAssignmentsWithKmeansOutputNameCSV = paste0(wgcnaAndKMeansPathsForOutputs, "//WGCNA_withKMeans_",netType, "_power",powerEstimate, "_ModuleAssignments", outputAddOn, ".csv")
geneAssignmentsWithKmeansOutputNameRData = paste0(wgcnaAndKMeansPathsForOutputs, "//WGCNA_withKMeans_",netType, "_power",powerEstimate, "_ModuleAssignments", outputAddOn, ".RData")

moduleCountsWithKmeansOutputNameCSV = paste0(wgcnaAndKMeansPathsForOutputs, "//WGCNA_withKMeans_",netType, "_power",powerEstimate, "_ModuleCounts", outputAddOn, ".csv")


#print(paste0("Please note that the Soft Threshold Power Data Frame has been written here: ", sftThresholdingOutputNameCSV))



geneModMembershipOutputNameCSV = paste0(wgcnaAndKMeansPathsForOutputs, "//WGCNA_withKMeans_",netType, "_power",powerEstimate, "_geneModuleMembership", outputAddOn, ".csv")
mmpOutputNameCSV = paste0(wgcnaAndKMeansPathsForOutputs, "//WGCNA_withKMeans_",netType, "_power",powerEstimate, "_MMPvalues", outputAddOn, ".csv")

modTraitCorrelationOutputNameCSV = paste0(wgcnaAndKMeansPathsForOutputs, "//WGCNA_withKMeans_",netType, "_power",powerEstimate, "_moduleTraitCorrelation", outputAddOn, ".csv")
modTraitCorrPValOutputNameCSV = paste0(wgcnaAndKMeansPathsForOutputs, "//WGCNA_withKMeans_",netType, "_power",powerEstimate, "_moduleTraitCorrPValue", outputAddOn, ".csv")


modTraitCorrAndCorrPValueOutputNameRData = paste0(wgcnaAndKMeansPathsForOutputs, "//WGCNA_withKMeans_",netType, "_power",powerEstimate, "_moduleTraitCorrelationAndCorrPValue", outputAddOn, ".RData")
signifModTraitCorrelationOutputNameCSV = paste0(wgcnaAndKMeansPathsForOutputs, "//WGCNA_withKMeans_",netType, "_power",powerEstimate, "_StatisticallySignifModuleTraitCorrelations_pValThresh", pValueCutOffForSignificance, outputAddOn, ".csv")


updatedModuleTraitPhenotypeDFFilePath = paste0(wgcnaAndKMeansPathsForOutputs, "//", bodyRegion, "_wgcna_with_kmeans_Power", powerEstimate, "_StatisticallySignificantModuleTraitCorrelations_CombinedInformation",outputAddOn,".csv")


wgcnaWithKmeansAllObjectsExceptTOM = paste0(wgcnaAndKMeansPathsForOutputs, "//", bodyRegion, "_wgcna_with_kmeans_Power", powerEstimate, "_ALLKeyObjectsExceptTOM",outputAddOn,".RData")
wgcnaWithKmeansAllObjectsIncludingTOM = paste0(wgcnaAndKMeansPathsForOutputs, "//", bodyRegion, "_wgcna_with_kmeans_Power", powerEstimate, "_ALLKeyObjectsIncludingTOM",outputAddOn,".RData")


wgcnaWithKmeansAllObjectsExceptTOMWithDate = paste0(wgcnaAndKMeansPathsForOutputs, "//", bodyRegion, "_wgcna_with_kmeans_Power", powerEstimate, "_ALLKeyObjectsExceptTOM",outputAddOnWithDate,".RData")
wgcnaWithKmeansAllObjectsIncludingTOMWithDate = paste0(wgcnaAndKMeansPathsForOutputs, "//", bodyRegion, "_wgcna_with_kmeans_Power", powerEstimate, "_ALLKeyObjectsIncludingTOM",outputAddOnWithDate,".RData")


wgcnaWithKmeansGeneTraitSignificancePath = paste0(wgcnaAndKMeansPathsForOutputs, "//", bodyRegion, "_wgcna_with_kmeans_Power", powerEstimate, "_GeneTraitSignificanceCorrelation",outputAddOn,".csv")
wgcnaWithKmeansGeneTraitSigPValuePath = paste0(wgcnaAndKMeansPathsForOutputs, "//", bodyRegion, "_wgcna_with_kmeans_Power", powerEstimate, "_GeneTraitCorrPValue",outputAddOn,".csv")


wgcnaWithKmeansGeneModuleMembershipPath  = paste0(wgcnaAndKMeansPathsForOutputs, "//", bodyRegion, "_wgcna_with_kmeans_Power", powerEstimate, "_GeneModuleMembershipCorrelation",outputAddOn,".csv")
wgcnaWithKmeansGeneModMemberPValuePath= paste0(wgcnaAndKMeansPathsForOutputs, "//", bodyRegion, "_wgcna_with_kmeans_Power", powerEstimate, "_GeneModMembershipCorrPValue",outputAddOn,".csv")


wgcnaWithKmeansGeneTraitsEdgeTablePathCSV  = paste0(wgcnaAndKMeansPathsForOutputs, "//", bodyRegion, "_wgcna_with_kmeans_Power", powerEstimate, "_GeneTraitCorrelationAndPValueTable",outputAddOn,".csv")
wgcnaWithKmeansGeneTraitsEdgeTablePathRData  = paste0(wgcnaAndKMeansPathsForOutputs, "//", bodyRegion, "_wgcna_with_kmeans_Power", powerEstimate, "_GeneTraitCorrelationAndPValueTable",outputAddOn,".RData")
signifGeneTraitCorrelationOutputNameCSV = paste0(wgcnaAndKMeansPathsForOutputs, "//", bodyRegion, "_wgcna_with_kmeans_Power", powerEstimate, "_SignificantGeneTraitCorrelationAndPValueTable_forPValThresh", pValueCutOffForSignificance,outputAddOn,".csv")
signifPositiveGeneTraitCorrelationOutputNameCSV = paste0(wgcnaAndKMeansPathsForOutputs, "//", bodyRegion, "_wgcna_with_kmeans_Power", powerEstimate, "_Significant_Positive_GeneTraitCorrelationAndPValueTable_forPValThresh", pValueCutOffForSignificance,outputAddOn,".csv")
updatedGeneTraitPhenotypeDFFilePath = paste0(wgcnaAndKMeansPathsForOutputs, "//", bodyRegion, "_wgcna_with_kmeans_Power", powerEstimate, "finalCombined_Significant_GeneTraitCorrelationAndPValueTable_forPValThresh", pValueCutOffForSignificance,outputAddOn,".csv")

wgcnaFileInfoDF = rbind(wgcnaFileInfoDF, c("bodyRegion:", bodyRegion))
wgcnaFileInfoDF = rbind(wgcnaFileInfoDF, c("tissueName:", tissueName))
wgcnaFileInfoDF = rbind(wgcnaFileInfoDF, c("disease:", disease))
wgcnaFileInfoDF = rbind(wgcnaFileInfoDF, c("fullDiseaseName:", fullDiseaseName))
wgcnaFileInfoDF = rbind(wgcnaFileInfoDF, c("minModuleSize:", minModuleSize))
wgcnaFileInfoDF = rbind(wgcnaFileInfoDF, c("maxNumPowers", maxNumPowers))
wgcnaFileInfoDF = rbind(wgcnaFileInfoDF, c("min.genes.for.grey:", min.genes.for.grey))
wgcnaFileInfoDF = rbind(wgcnaFileInfoDF, c("max # of kMeans Iterations", nIterations))
wgcnaFileInfoDF = rbind(wgcnaFileInfoDF, c("meg", meg))
wgcnaFileInfoDF = rbind(wgcnaFileInfoDF, c("network type (netType)", netType))
wgcnaFileInfoDF = rbind(wgcnaFileInfoDF, c("numRoundingDigits:", numRoundingDigits))
wgcnaFileInfoDF = rbind(wgcnaFileInfoDF, c("pValueCutOffForSignificance:", pValueCutOffForSignificance))
wgcnaFileInfoDF = rbind(wgcnaFileInfoDF, c("use Recommended Power by WGCNA", useRecommendedPower))
wgcnaFileInfoDF = rbind(wgcnaFileInfoDF, c("ourPower:", ourPower) )
wgcnaFileInfoDF = rbind(wgcnaFileInfoDF, c("positiveCorrelationsOnly?", positiveCorrelationsOnly))
wgcnaFileInfoDF = rbind(wgcnaFileInfoDF, c("inputGeneExpressionDataFilePath:", inputGeneExpressionDataFilePath))
wgcnaFileInfoDF = rbind(wgcnaFileInfoDF, c("wgcnaAndKMeansOutputPath:", wgcnaAndKMeansOutputPath))
wgcnaFileInfoDF = rbind(wgcnaFileInfoDF, c("geneToEntrezIDMappingPath:", geneToEntrezIDMappingPath))
wgcnaFileInfoDF = rbind(wgcnaFileInfoDF, c("phenotypesFilePath:", phenotypesFilePath))
wgcnaFileInfoDF = rbind(wgcnaFileInfoDF, c("log2transformInputData:", log2transformInputData))
wgcnaFileInfoDF = rbind(wgcnaFileInfoDF, c("scaleInputData:", scaleInputData))
wgcnaFileInfoDF = rbind(wgcnaFileInfoDF, c("dataScaling:", dataScaling))
wgcnaFileInfoDF = rbind(wgcnaFileInfoDF, c("geneToEntrezIDMapping RData File:", geneToEntrezIDPathForGeneExpressionDataRData))
wgcnaFileInfoDF = rbind(wgcnaFileInfoDF, c("geneToEntrezIDMapping CSV File:", geneToEntrezIDPathForGeneExpressionDataCSV))
wgcnaFileInfoDF = rbind(wgcnaFileInfoDF, c("original gene expression RData File (before filtering):", originalGeneFilteredRData))
wgcnaFileInfoDF = rbind(wgcnaFileInfoDF, c("original gene expression CSV File (before filtering):", originalGeneFilteredCSV))
wgcnaFileInfoDF = rbind(wgcnaFileInfoDF, c(paste0("final gene expression CSV File (after filtering) and ", dataScaling, " with date:"), geneExpressionOutputNameCSV_withDate))
wgcnaFileInfoDF = rbind(wgcnaFileInfoDF, c(paste0("final gene expression CSV File (after filtering) and ", dataScaling, " :"), updatedOutputAddOn))
wgcnaFileInfoDF = rbind(wgcnaFileInfoDF, c(paste0("final gene expression RData File (after filtering) and ", dataScaling, ":"), geneExpressionOutputNameRData_withDate))
wgcnaFileInfoDF = rbind(wgcnaFileInfoDF, c(paste0("final DatExpr CSV File (after filtering) and ", dataScaling, " with date:"), updatedOutputAddOnAndDatExpr_withDate))
wgcnaFileInfoDF = rbind(wgcnaFileInfoDF, c(paste0("final DatExpr CSV File (after filtering) and ", dataScaling, " :"), updatedOutputAddOnAndDatExpr))

wgcnaFileInfoDF

#source("D:\\organizedAlzheimers\\Setup\\parametersDerivedFromUserInput.R")



write.csv(dissTOMPower, dissTOMPowerOutputNameCSV)
save(dissTOMPower,  file = dissTOMPowerOutputNameRData)
save(TOMPower,  file = tomPowerOutputNameRData)


#dissTOMPowerOutputName = paste0(outputPath, "//dissTOMPower_",netType, "_power",powerEstimate, outputAddRData)
#wgcnaOutputName = paste0(outputPath, "//WGCNAInitialResultsObjects_",netType, "_power",powerEstimate, outputAddRData)
write.csv(TOMPower, file = tomPowerOutputName) #"D://dlpfc_dissTOMPower_newRNASeq_Alzh_1029NEW.RData")

#write.csv(dissTOMPower, file = dissTOMPowerOutputName) #"D://dlpfc_dissTOMPower_newRNASeq_Alzh_1029NEW.RData")
#rm(TOMPower)
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

#outputFileName = paste0(parentPath, "wgcna_initialResultsPower_newRNASeqData_1029_dlpfc_power_", powerEstimate, ".RData")
save(dissTOMPower, MEList, MEs, MEDiss, geneTreePower, powerEstimate,  dynamicModsPower, dynamicColorsPower, 
     file = wgcnaPowerOutputNameRData) #outputFileName)



#data.frame(table(dynamicColorsPower))

#moduleCountsOutputName = paste0(outputPath, "//WGCNAInitialModuleCounts_",netType, "_power",powerEstimate, outputAddOn)

write.csv(data.frame(table(dynamicColorsPower)), moduleCountsOutputName) #"D:\\initialDLPFC_modules.csv")

#####################
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = paste("Initial Clustering of Module Eigengenes (", diseaseName, "; ", bodyRegion,")"),
     xlab = paste0(netType, " Network with Power: ", powerEstimate), sub = "")
#head(datExpr)

initialModuleCounts = data.frame(table(dynamicColorsPower))
print(initialModuleCounts)



# Please calculate initial Module Eigengenes

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

print("Getting initial eigengenes")
if(sum(partition.in.colors == "grey") < min.genes.for.grey){
  eigengenes = moduleEigengenes(datExpr,partition.in.colors, excludeGrey=TRUE)
} else {
  eigengenes = moduleEigengenes(datExpr,partition.in.colors, excludeGrey=F)
}
cat("We got",length(eigengenes$eigengenes)," eigengene vectors\n")
print(head(eigengenes$eigengenes))

#This variable is fixed and used as a reference to indicate the
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
  print(paste0("Number of centroids before getting new partition ",ncol(centroids)))
  
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

data.frame(table(dynamicColorsKmeans))
dim(datExpr)
MEList = moduleEigengenes(datExpr, colors = dynamicColorsKmeans)
MEs = MEList$eigengenes


write.csv(MEs, MEsWithKmeansOutputName) #paste0(outputPath, "dlpfc_1029_wgcna_with_kmeans_Power", powerEstimate, "ModuleEigengenesMEs_newRNASeqAlzh.csv"))

gene_names = colnames(datExpr)
clusterInfo = data.frame(module = dynamicColorsKmeans)
rownames(clusterInfo)=gene_names
dim(clusterInfo)
clusterInfo = merge(clusterInfo, geneAndEntrezIDMappingDF, by = 0)
rownames(clusterInfo)=gene_names
clusterInfo = clusterInfo[,-1]
colnames(clusterInfo) = c("module", "entrezID")
head(clusterInfo)


write.csv(clusterInfo, geneAssignmentsWithKmeansOutputNameCSV) #paste0(outputPath, "wgcna_1029_with_kmeans_Power", powerEstimate, "_moduleAssignments__newRNASeqAlzh.csv"))
save(clusterInfo, file = geneAssignmentsWithKmeansOutputNameRData) #"D:\\DLPFC_updatedGeneClusterInfo.RData")
data.frame(table(clusterInfo))

write.csv(data.frame(table(clusterInfo$module)), moduleCountsWithKmeansOutputNameCSV) #"D:\\DLPFCModules.csv")

dim(datExpr)
nSamples = nrow(datExpr)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
#geneModMembershipOutputNameCSV = paste0(outputPath, "//WGCNA_withKMeans_",netType, "_power",powerEstimate, "_geneModuleMembership", outputAddOn)
#mmpOutputNameCSV = paste0(outputPath, "//WGCNA_withKMeans_",netType, "_power",powerEstimate, "_MMPvalues", outputAddOn)

geneModuleMembership = merge(geneAndEntrezIDMappingDF, geneModuleMembership, by = 0)
rownames(geneModuleMembership)=gene_names
geneModuleMembership = geneModuleMembership[,-1]
colnames(geneModuleMembership)[1] = "entrezID"
MMPvalue = merge(geneAndEntrezIDMappingDF, MMPvalue, by = 0)
colnames(MMPvalue)[1] = "entrezID"
rownames(MMPvalue)=gene_names
MMPvalue = MMPvalue[,-1]
head(MMPvalue)
#colnames(clusterInfo) = c("module", "entrezID")
head(clusterInfo)


write.csv(geneModuleMembership, geneModMembershipOutputNameCSV)

write.csv(MMPvalue, mmpOutputNameCSV)

head(diseaseTraits)
colnames(diseaseTraits)
str(diseaseTraits)

# chr

class(diseaseTraits[,2])

typesToRemoveVec = c(1)
traitsToKeepVec = c()
for (k in 3:ncol(diseaseTraits)){
  if (class(diseaseTraits[,k]) == "character"){
    typesToRemoveVec = c(typesToRemoveVec, k)
  } else {
    traitsToKeepVec = c(traitsToKeepVec, k)
  }
  
}
typesToRemoveVec
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
traitRows = match(alzhSamples, allTraits[,1]);


#names(binaryTraits)
#[1] "Patient"       "controlStage"  "initialStage"  "moderateStage" "severeStage" 

datTraits = allTraits[traitRows, -1];

rownames(datTraits) = allTraits[traitRows, 1];
datTraits


moduleTraitCor = cor(MEs, datTraits, use = "p") #, robustY = FALSE, maxPOutliers = 0.1);

#modTraitCorrelationOutputNameCSV = paste0(outputPath, "//WGCNA_withKMeans_",netType, "_power",powerEstimate, "_moduleTraitCorrelation", outputAddOn)
#modTraitCorrPValOutputNameCSV = paste0(outputPath, "//WGCNA_withKMeans_",netType, "_power",powerEstimate, "_moduleTraitCorrPValue", outputAddOn)

write.csv(data.frame(moduleTraitCor), modTraitCorrelationOutputNameCSV)

#names(datTraits[,otherTraits]) # [1] "Age"  "nft"  "mmse" "pmi" 
nSamples = nrow(datExpr)
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
#write.csv(data.frame(moduleTraitPvalue), "C://Users//sk792//Documents//Summer 2020//Alzheimers Project//mostUpdatedAlzheimers//mergedGroupsPow44_line0.10_moduleTraitPvalue.csv")
write.csv(data.frame(moduleTraitPvalue), modTraitCorrPValOutputNameCSV)
currentDate = as.character(now()) # current date and time str_replace_all(Sys.Date(), "-", "_")

save(moduleTraitCor, moduleTraitPvalue, currentDate, file = modTraitCorrAndCorrPValueOutputNameRData)
















library(reshape2)
info = "WGCNA with kMeans Applied"
modTraitDF1 = melt(moduleTraitCor)
modTraitDF1$newName = paste0(modTraitDF1[,1], "_", modTraitDF1[,2])

#modTraitDF1 = modTraitDF1[,c(4,3)]
#colnames(modTraitDF1) = c("comboName", "correlation_r")
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
  write.csv(significantModTraits, paste0(outputPath, "wgcna_with_kmeans_Power", powerEstimate, "_StatisticallySignificantModuleTraitCorrelations_POSITIVEonly_.csv"))
  
}
#numRoundingDigits = 3
significantModTraits$roundedCorr = round(significantModTraits$correlation_r, numRoundingDigits)
#miniDF$roundedCorr = round(significantModTraits$correlation_r, numRoundingDigits)

head(significantModTraits)
uniqueModules = as.vector(unique(significantModTraits[,2])) #
moduleNamesVec = unique(dynamicColorsKmeans) 
#as.vector(unique(significantModTraits[,2]))
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





save(significantModTraits,
     updatedModuleTraitPhenotypeDF,
     bodyRegion,
     currentDate,
     disease,
     dataScaling,
     datExpr,
     moduleTraitCor,
     moduleTraitPvalue,
     geneModuleMembership,
     MMPvalue,
     clusterInfo,
     dynamicColorsKmeans,
     MEList,
     powerEstimate,
     info,
     geneTraitSignificanceDF,
     GSPvalueDF,
     file = wgcnaWithKmeansAllObjectsExceptTOM)

save(significantModTraits,
     updatedModuleTraitPhenotypeDF,
     bodyRegion,
     tissueName,
     dataScaling,
     currentDate,
     disease,
     datExpr,
     moduleTraitCor,
     moduleTraitPvalue,
     geneModuleMembership,
     MMPvalue,
     clusterInfo,
     dynamicColorsKmeans,
     MEList,
     TOMPower,
     powerEstimate,
     info,
     GSPvalueDF,
     geneTraitSignificanceDF,
     file = wgcnaWithKmeansAllObjectsIncludingTOM)


uniqueModules = unique(dynamicColorsKmeans)


allTOM = TOM#[inModule, inModule];
allTOM[upper.tri(allTOM)] <- NA
probes = names(datExpr)
dimnames(allTOM) = list(probes, probes)

allTOMInteractionsDF = reshape2::melt(allTOM, varnames = c('row', 'col'), na.rm = TRUE)
mainDiagIndices = which(allTOMInteractionsDF[,1] == allTOMInteractionsDF[,2])
moduleInteractionsDF = moduleInteractionsDF[-mainDiagIndices,]


allTOMInteractions = melt(TOMPower)
head(allTOMInteractions)
inModule = (moduleColors==module);
modProbes = probes#[inModule];
# Select the corresponding Topological Overlap
modTOM = TOMPower#[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
modTOM[upper.tri(modTOM)] <- NA
cbind(which(!is.na(m),arr.ind = TRUE),na.omit(as.vector(m)))


allTOMInteractions = melt(modTOM)

interactions = melt(modTOM)

hubs  = chooseTopHubInEachModule(datExpr, dynamicColorsKmeans)

#   quantile(moduleInteractionsDF$TOMvalue, probs = percentilesForCoexpressionNetworkVec)
#moduleToUse = uniqueModules[1]

allModulesCoExpressionNetworkModularTOMResultsDF = data.frame() # only for modular edges 
allNetworkModulesTomSummaryStatisticsDF = data.frame()
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
#moduleCoExpressionResultsList1 = pleaseGetCoexpressionNetworkForAWGCNAModuleAfterKmeans(module = moduleToUse)




geneTraitSignificanceResultsList = pleaseGetTheGeneTraitSignificanceValuesForTheTraits(datTraits,datExpr,
                                                                                       wgcnaWithKmeansGeneTraitSignificancePath,
                                                                                       wgcnaWithKmeansGeneTraitSigPValuePath,
                                                                                       bodyRegion, tissueName, fullDiseaseName, dataScaling,
                                                                                       geneAndEntrezIDMappingDF)
geneTraitSignificanceDF = data.frame(geneTraitSignificanceResultsList$geneTraitSignificanceDF)
GSPvalueDF = data.frame(geneTraitSignificanceResultsList$GSPvalueDF )

dim(geneTraitSignificanceDF)
dim(GSPvalueDF)




pleaseGetTheGeneModuleMembershipForTheGenes <- function(MEs, powerEstimate,
                                                        datExpr,
                                                        wgcnaWithKmeansGeneTraitSignificancePath,
                                                        wgcnaWithKmeansGeneTraitSigPValuePath,
                                                        bodyRegion, tissueName, fullDiseaseName, dataScaling, geneAndEntrezIDMappingDF) {
  
  nSamples = nrow(datExpr)
  gene_names = colnames(datExpr)
  
  geneModuleMembershipDF = as.data.frame(cor(datExpr, MEs, use = "p"));
  MMPvalueDF = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembershipDF), nSamples));
  modNames = substring(names(MEs), 3)
  names(geneModuleMembershipDF) = paste("GeneModuleMembership_", modNames, sep="");
  names(MMPvalueDF) = paste("pValue_GeneModMembership_", modNames, sep="");
  
  geneModuleMembershipDF$explanation = paste0("Module Membership: Correlation (r) between gene ", 
                                              colnames(datExpr), " and the Module Eigengene of Each of the ", ncol(MEs), " gene modules for ",
                                              disease, " trait observed in the ", bodyRegion, " of the ",
                                              tissueName)
  geneModuleMembershipDF$dataRelatedExplanation = rep(paste0("Please note the gene expression data had a ", dataScaling, 
                                                             " transformation and WGCNA with kMeans was applied to a ",
                                                             netType, " network with soft threshold power: ", powerEstimate, "."),
                                                      nrow(geneModuleMembershipDF))
  
  MMPvalueDF$explanation = paste0("Module Membership P-Value: Associated P-value for the Correlation between gene ", 
                                  colnames(datExpr), " and the Module Eigengene of Each of the ", ncol(MEs), " gene modules for ",
                                  disease, " trait observed in the ", bodyRegion, " of the ",
                                  tissueName)
  MMPvalueDF$dataRelatedExplanation = rep(paste0("Please note the gene expression data had a ", dataScaling, 
                                                 " transformation and WGCNA with kMeans was applied to a ",
                                                 netType, " network with soft threshold power: ", powerEstimate, "."),
                                          nrow(geneModuleMembershipDF))
  
  
  
  geneModuleMembershipDF = merge(geneAndEntrezIDMappingDF, geneModuleMembershipDF, by = 0)
  rownames(geneModuleMembershipDF)=gene_names
  geneModuleMembershipDF = geneModuleMembershipDF[,-1]
  colnames(geneModuleMembershipDF)[1] = "entrezID"
  
  MMPvalueDF = merge(geneAndEntrezIDMappingDF, MMPvalueDF, by = 0)
  rownames(MMPvalueDF)=gene_names
  MMPvalueDF = MMPvalueDF[,-1]
  colnames(MMPvalueDF)[1] = "entrezID"
  
  head(MMPvalueDF)
  
  
  write.csv(geneModuleMembershipDF, wgcnaWithKmeansGeneModuleMembershipPath)
  
  write.csv(MMPvalueDF, wgcnaWithKmeansGeneModMemberPValuePath)
  
  geneModuleMembershipResultsList = list()
  geneModuleMembershipResultsList$geneModuleMembershipDF = geneModuleMembershipDF
  geneModuleMembershipResultsList$MMPvalueDF = MMPvalueDF
  print(paste0(":) Please note we return the results for gene trait significance correlations and associated correlation p-values here: ",
               "geneModuleMembershipResultsList"))
  
  print(paste0("Please note :) geneModuleMembershipResultsList$geneModuleMembershipDF has the geneModuleMembership correlations"))
  print(paste0("Please note :) geneModuleMembershipResultsList$MMPvalueDF has the p-values for the geneModuleMembership correlations"))
  print(str(geneModuleMembershipResultsList))
  return(geneModuleMembershipResultsList)
}

geneModuleMembershipResultsList = pleaseGetTheGeneModuleMembershipForTheGenes(MEs, powerEstimate,datExpr,
                                                                              wgcnaWithKmeansGeneTraitSignificancePath,
                                                                              wgcnaWithKmeansGeneTraitSigPValuePath,
                                                                              bodyRegion, tissueName, fullDiseaseName, dataScaling, geneAndEntrezIDMappingDF)

geneModuleMembershipDF = data.frame(geneModuleMembershipResultsList$geneModuleMembershipDF)
MMPvalueDF = data.frame(geneModuleMembershipResultsList$MMPvalueDF)
# gene mod membership...

dim(geneModuleMembershipDF)
dim(MMPvalueDF)