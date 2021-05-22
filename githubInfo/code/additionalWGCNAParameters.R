
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


tfsUsed = "LambertAndJaspar"



