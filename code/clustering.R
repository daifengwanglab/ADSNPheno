source("D:\\organizedAlzheimers\\Setup\\parametersDerivedFromUserInput.R")

geneAndEntrezIDMappingDF = geneExpressionDataSetToEntrezIDMapping(inputGeneExpressionDataFilePath,
                                                   geneToEntrezIDPathForGeneExpressionDataCSV, 
                                                   geneToEntrezIDPathForGeneExpressionDataRData)

#geneAndEntrezIDMappingDF$geneName = row.names(geneAndEntrezIDMappingDF)

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



powerEstimate = selectingWGCNAPower(datExpr = datExpr, maxNumPowers = maxNumPowers,
                                netType = netType, ourPower = ourPower, useRecommendedPower = useRecommendedPower,
                                log2transformInputData = log2transformInputData, scaleInputData = scaleInputData) 


require(parallel)
require(doParallel)

cpucores <- makeCluster(detectCores(), type='PSOCK') ; 

registerDoParallel(cpucores) ; 

Sys.setenv("MC_CORES"=cpucores) ;
allowWGCNAThreads(nThreads = 20)

# a check here to see if we wanted to use recommendations or our own power
adjacencyPower = adjacency(datExpr, power = powerEstimate, type = netType)
# "Please note that the Soft-Threshold Power Selected is: 11"
# --> NEW "Please note that the Soft-Threshold Power Selected is: 10"
# Please Turn adjacency into topological overlap :) matrix (TOM)
TOMPower = TOMsimilarity(adjacencyPower);
#save(TOMPower,  file = tomPowerOutputNameRData)
rm(adjacencyPower)
dissTOMPower = 1 - TOMPower


#source("D:\\organizedAlzheimers\\Setup\\parametersDerivedFromUserInput.R")
source("D:\\organizedAlzheimers\\Setup\\additionalWGCNAParameters.R")


save(dissTOMPower,  file = dissTOMPowerOutputNameRData)



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
# save(MEList, MEs, MEDiss, geneTreePower, powerEstimate,  dynamicModsPower, dynamicColorsPower, 
#      file = wgcnaPowerOutputNameRData)


#data.frame(table(dynamicColorsPower))

#moduleCountsOutputName = paste0(outputPath, "//WGCNAInitialModuleCounts_",netType, "_power",powerEstimate, outputAddOn)

write.csv(data.frame(table(dynamicColorsPower)), moduleCountsOutputName) #"D:\\initialDLPFC_modules.csv")

#####################
# Cluster module eigengenes
diseaseName = disease #'Alzheimers'
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

#write.csv(MEs,"D://keyMEsLTL_March28th.csv")
write.csv(MEs, MEsWithKmeansOutputName) #paste0(outputPath, "dlpfc_1029_wgcna_with_kmeans_Power", powerEstimate, "ModuleEigengenesMEs_newRNASeqAlzh.csv"))

gene_names = colnames(datExpr)
clusterInfo = data.frame(module = dynamicColorsKmeans)
rownames(clusterInfo)=gene_names
dim(clusterInfo)
clusterInfo = merge(clusterInfo, geneAndEntrezIDMappingDF, by = 0)
gene_names = clusterInfo[,1] # newer Chanu :)
rownames(clusterInfo)= clusterInfo[,1]#gene_names
head(clusterInfo)
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

#gene_names = 
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
#traitRows = match(alzhSamples, allTraits[,1]);
traitRows = match(alzhSamples, diseaseTraits[,2]);
datTraits = allTraits[traitRows, -1];
length(traitRows)
traitRows = traitRows[which(!isNA(traitRows))]
length(traitRows)

#names(binaryTraits)
#[1] "Patient"       "controlStage"  "initialStage"  "moderateStage" "severeStage" 
#rownames(datTraits) = allTraits[traitRows, 1];
#datTraits = allTraits[traitRows, -1];

rownames(datTraits) = diseaseTraits[traitRows, 1];

#

datTraits
#MEs = read.csv("D://organizedAlzheimers//GeneExpressionPreparation//outputs//neuroscience//WGCNA_withKMeans_signed_power30_ModuleEigengenesMEs_OriginalInputData.csv", header = TRUE, row.names = 1)

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






head(clusterInfo)









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
numRoundingDigits = numOfRoundingDigits
#numRoundingDigits = 3
significantModTraits$roundedCorr = round(significantModTraits$correlation_r, numRoundingDigits)
#miniDF$roundedCorr = round(significantModTraits$correlation_r, numRoundingDigits)

head(significantModTraits)
uniqueModules = as.vector(unique(significantModTraits[,2])) #
# head(clusterInfo)
# dynamicColorsKmeans = clusterInfo[,1]
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



# 
# 
# save(significantModTraits,
#      updatedModuleTraitPhenotypeDF,
#      bodyRegion,
#      currentDate,
#      disease,
#      dataScaling,
#      datExpr,
#      moduleTraitCor,
#      moduleTraitPvalue,
#      geneModuleMembership,
#      MMPvalue,
#      clusterInfo,
#      dynamicColorsKmeans,
#      MEList,
#      powerEstimate,
#      info,
#      geneTraitSignificanceDF,
#      GSPvalueDF,
#      file = wgcnaWithKmeansAllObjectsExceptTOM)
# 
# save(significantModTraits,
#      updatedModuleTraitPhenotypeDF,
#      bodyRegion,
#      tissueName,
#      dataScaling,
#      currentDate,
#      disease,
#      datExpr,
#      moduleTraitCor,
#      moduleTraitPvalue,
#      geneModuleMembership,
#      MMPvalue,
#      clusterInfo,
#      dynamicColorsKmeans,
#      MEList,
#      TOMPower,
#      powerEstimate,
#      info,
#      GSPvalueDF,
#      geneTraitSignificanceDF,
#      file = wgcnaWithKmeansAllObjectsIncludingTOM)


uniqueModules = unique(dynamicColorsKmeans)
print(paste0(":) Please note that there are: ", length(uniqueModules), " total unique gene co-expression modules here!"))
#modDF= read.csv("D://organizedAlzheimers//GeneExpressionPreparation//outputs//neuroscience//WGCNA_withKMeans_signed_power30_ModuleAssignments_OriginalInputData.csv", header = TRUE)
#uniqueModules = unique(modDF[,2])
#dynamicColorsKmeans = modDF[,2]
allTOM = 1 - dissTOMPower # allTOM = TOMPower#[inModule, inModule];
allTOM[upper.tri(allTOM)] <- NA
probes = names(datExpr)
dimnames(allTOM) = list(probes, probes)

allTOMInteractionsDF = reshape2::melt(allTOM, varnames = c('row', 'col'), na.rm = TRUE)
mainDiagIndices = which(allTOMInteractionsDF[,1] == allTOMInteractionsDF[,2])
moduleInteractionsDF = moduleInteractionsDF[-mainDiagIndices,]

###################################
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
#######################################
#dynamicColorsKmeans = clusterInfo[,1]
#dynamicColorsKmeans = str_replace_all(vectorInfo, "_47", "")
#dfInfo = read.csv("D:\\organizedAlzheimers\\GeneExpressionPreparation\\outputs\\neuroscience\\WGCNA_withKMeans_signed_power30_ModuleAssignments_OriginalInputData.csv", header = TRUE)
#dynamicColorsKmeans= dfInfo[,2]
hubs  = chooseTopHubInEachModule(datExpr, dynamicColorsKmeans)
# hubsDF = data.frame(hubs)
# hubsDF$power = rep(powerEstimate,nrow(hubsDF))
# hubsDF$region = rep(bodyRegion,nrow(hubsDF))
# hubsDF$date = rep("May 13, 2021",nrow(hubsDF))
# write.csv(hubsDF, "D://hubsNeuroscience_Power30_updatedMay13.csv")

#write.csv(data.frame(colnames(datExpr)), "D://neuroOrder.csv")

######colorh = labels2colors(dynamicColorsKmeans)
###### hubs    = chooseTopHubInEachModule(datExpr, colorh)
######head(colorh)
head(dynamicColorsKmeans)

#write.csv(colnames(datExpr), "D://neuroscienceHubs.csv")
hubsDF = data.frame(hubs)
hubsDF$power = rep(powerEstimate,nrow(hubsDF))
hubsDF$region = rep(bodyRegion,nrow(hubsDF))
head(hubsDF)
dim(hubsDF)
write.csv(hubsDF, "D://hubsNeuroscience_Power30_May13thNewest.csv")


#   quantile(moduleInteractionsDF$TOMvalue, probs = percentilesForCoexpressionNetworkVec)
#moduleToUse = uniqueModules[1]

allModulesCoExpressionNetworkModularTOMResultsDF = data.frame() # only for modular edges 
allNetworkModulesTomSummaryStatisticsDF = data.frame()
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
#moduleCoExpressionResultsList1 = pleaseGetCoexpressionNetworkForAWGCNAModuleAfterKmeans(module = moduleToUse)

allNetworkModulesCombinedTomSummaryStatisticsFileName = paste0(wgcnaAndKMeansPathsForOutputs, "//", bodyRegion, "_wgcna_with_kmeans_Power", powerEstimate, "finalCombined_CoExpressionNetwork_forAll_", length(uniqueModules), "_Modules",outputAddOn,".csv")
print(paste0(":) please note we are writing out the Co-Expression Network for ALL ", length(uniqueModules), " here:", allNetworkModulesCombinedTomSummaryStatisticsFileName))
write.csv(allNetworkModulesTomSummaryStatisticsDF, allNetworkModulesCombinedTomSummaryStatisticsFileName)



geneTraitSignificanceResultsList = pleaseGetTheGeneTraitSignificanceValuesForTheTraits(datTraits,datExpr,
                                                                                                   wgcnaWithKmeansGeneTraitSignificancePath,
                                                                                                   wgcnaWithKmeansGeneTraitSigPValuePath,
                                                                                       bodyRegion, tissueName, fullDiseaseName, dataScaling,
                                                                                       geneAndEntrezIDMappingDF)
geneTraitSignificanceDF = data.frame(geneTraitSignificanceResultsList$geneTraitSignificanceDF)
GSPvalueDF = data.frame(geneTraitSignificanceResultsList$GSPvalueDF )

dim(geneTraitSignificanceDF)
dim(GSPvalueDF)

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
     hubs,
     allNetworkModulesTomSummaryStatisticsDF,
     powerEstimate,
     info,
     geneTraitSignificanceDF,
     GSPvalueDF,
     file = wgcnaWithKmeansAllObjectsExceptTOM)


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
gene_names = geneModuleMembershipDF[,1]
rownames(geneModuleMembershipDF)= gene_names
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
#write.csv(geneModuleMembershipDF, wgcnaWithKmeansGeneModuleMembershipPath)

#write.csv(MMPvalueDF, wgcnaWithKmeansGeneModMemberPValuePath)
dim(geneModuleMembershipDF)
dim(MMPvalueDF)


save(significantModTraits,
     updatedModuleTraitPhenotypeDF,
     bodyRegion,
     currentDate,
     datTraits,
     geneAndEntrezIDMappingDF,
     disease,
     dataScaling,
     datExpr,
     moduleTraitCor,
     moduleTraitPvalue,
     geneModuleMembership,
     MMPvalue,
     clusterInfo,
     uniqueModules,
     dynamicColorsKmeans,
     MEList,
     powerEstimate,
     info,
     geneTraitSignificanceDF,
     GSPvalueDF,
     file = wgcnaWithKmeansAllObjectsExceptTOM)

save(significantModTraits,
     updatedModuleTraitPhenotypeDF,
     geneAndEntrezIDMappingDF,
     bodyRegion,
     tissueName,
     dataScaling,
     currentDate,
     datTraits,
     disease,
     datExpr,
     moduleTraitCor,
     moduleTraitPvalue,
     geneModuleMembership,
     MMPvalue,
     clusterInfo,
     dynamicColorsKmeans,
     uniqueModules,
     MEList,
     TOMPower,
     powerEstimate,
     info,
     GSPvalueDF,
     geneTraitSignificanceDF,
     file = wgcnaWithKmeansAllObjectsIncludingTOM)

#save(results_401_to_500, file= "D://results_401_to_500.RData")
#install.packages("tictoc")
#library("tictoc")
# 
# tic("starting")
# results_501_to_600 <- motifbreakR(snpList = snps.mb_10001_to_20000[501:600], 
#                                   filterp = TRUE,
#                                    pwmList = motifbreakR_motif,
#                                     threshold = motifBreakrThreshold, #1e-4,
#                                     method = "default", #ic",
#                                      bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
#                                      BPPARAM = BiocParallel::SnowParam(17))
# save(results_501_to_600, file = "D://results_501_to_600.RData") # 9432.72
# toc("ending")
# 
# 
# tic("starting")
# results_601_to_700 <- motifbreakR(snpList = snps.mb_10001_to_20000[601:700], 
#                                   filterp = TRUE,
#                                   pwmList = motifbreakR_motif,
#                                   threshold = motifBreakrThreshold, #1e-4,
#                                   method = "default", #ic",
#                                   bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
#                                   BPPARAM = BiocParallel::SnowParam(17))
# save(results_601_to_700, file = "D://results_601_to_700.RData") # 9432.72
# toc("ending") # 10,406.37 
# 
# 
# 
# 
# 
# 
# tic("starting")
# results_701_to_800 <- motifbreakR(snpList = snps.mb_10001_to_20000[701:800], 
#                                   filterp = TRUE,
#                                   pwmList = motifbreakR_motif,
#                                   threshold = motifBreakrThreshold, #1e-4,
#                                   method = "default", #ic",
#                                   bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
#                                   BPPARAM = BiocParallel::SnowParam(17))
# save(results_701_to_800, file = "D://results_701_to_800.RData") # 9432.72
# toc("ending") 
# 
# 
# 
# rm(results_701_to_800)
# 
# tic("starting")
# results_801_to_900 <- motifbreakR(snpList = snps.mb_10001_to_20000[801:900], 
#                                   filterp = TRUE,
#                                   pwmList = motifbreakR_motif,
#                                   threshold = motifBreakrThreshold, #1e-4,
#                                   method = "default", #ic",
#                                   bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
#                                   BPPARAM = BiocParallel::SnowParam(17))
# save(results_801_to_900, file = "D://results_801_to_900.RData") # 9432.72
# toc("ending") 
# 
# 
# rm(results_801_to_900)
# 
# tic("starting")
# results_901_to_1000 <- motifbreakR(snpList = snps.mb_10001_to_20000[901:1000], 
#                                   filterp = TRUE,
#                                   pwmList = motifbreakR_motif,
#                                   threshold = motifBreakrThreshold, #1e-4,
#                                   method = "default", #ic",
#                                   bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
#                                   BPPARAM = BiocParallel::SnowParam(17))
# save(results_901_to_1000, file = "D://results_901_to_1000.RData") # 9432.72
# toc("ending") 
# 




# to be safe, please save this:
save(significantModTraits,
     updatedModuleTraitPhenotypeDF,
     bodyRegion,
     currentDate,
     disease,
     dataScaling,
     datExpr,
     datTraits,
     moduleTraitCor,
     moduleTraitPvalue,
     geneModuleMembership,
     MMPvalue,
     clusterInfo,
     dynamicColorsKmeans,
     uniqueModules,
     MEList,
     geneAndEntrezIDMappingDF,
     powerEstimate,
     info,
     geneTraitSignificanceDF,
     GSPvalueDF,
     file = wgcnaWithKmeansAllObjectsExceptTOMWithDate)

save(significantModTraits,
     updatedModuleTraitPhenotypeDF,
     bodyRegion,
     tissueName,
     geneAndEntrezIDMappingDF,
     dataScaling,
     currentDate,
     disease,
     datExpr,
     datTraits,
     moduleTraitCor,
     moduleTraitPvalue,
     geneModuleMembership,
     MMPvalue,
     clusterInfo,
     dynamicColorsKmeans,
     uniqueModules,
     MEList,
     TOMPower,
     powerEstimate,
     info,
     GSPvalueDF,
     geneTraitSignificanceDF,
     file = wgcnaWithKmeansAllObjectsIncludingTOMWithDate)




