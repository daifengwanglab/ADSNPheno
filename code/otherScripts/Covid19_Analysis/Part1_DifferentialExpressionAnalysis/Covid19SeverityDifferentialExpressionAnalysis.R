#########################################################################################################################################
### INFORMATION:  Please note that the code below was written by Saniya Khullar and Daifeng Wang, Ph.D.
#########################################################################################################################################
### USER INPUTS:

#!/usr/bin/env Rscript

# Reading in all of the proper data files:
adsnphenoDirectory = "F://organizedAlzheimers//ADSNPheno//"  # please supply full directory path

# Please write the output locations for your files (where you would like them to be stored)
originalRawGeneExpressionDataBoxplotFilePathPNG = "F://covidAD_ICU_vs_NonICU_InfoRawDataBoxplots.png"
originalRawGeneExpressionDataBoxplotFilePathPDF = "F://covidAD_ICU_vs_NonICU_InfoRawDataBoxplots.pdf"
nomalizedGeneExpressionDataBoxplotFilePathPNG = "F://covidAD_ICU_vs_NonICU_InfoNormalizedDataBoxplots.png"
nomalizedGeneExpressionDataBoxplotFilePathPDF = "F://covidAD_ICU_vs_NonICU_InfoNormalizedDataBoxplots.pdf"
dataSetOutputPaths =  "F://covid19SeverityOutputs_"
volcanoPlotOutputFileName = "F://organizedAlzheimers//differentialExpressionAnalysisCovidAD//volcanoPlotADCovidICUvsNon-ICU_",


# Please note you can adjust these other parameters for the plots:
boxplotsPlotWidthPNG = 2501
boxplotsPlotHeightPNG = 1501
boxplotsPlotWidthPDF = 26
boxplotsPlotHeightPDF = 15
volcanoPlotWidth = 2001
volcanoPlotHeight = 1001

#################################################################################################################################
########### PLEASE RUN:
setwd(adsnphenoDirectory) # please set this folder as working directory for code and other resources
codeFolder = paste0(adsnphenoDirectory, "code//") # the main folder where code would be within ADSNPheno
dataFolder = paste0(adsnphenoDirectory, "data//") # the main folder where code would be within ADSNPheno

virus_data_raw = read.csv(paste0(dataFolder, "covidDataSets//covidForADProjectInfo.csv"), header = TRUE,check.names = FALSE)[,-1]

adCovidGRNGenesListFilePath = read.csv(paste0(dataFolder, "covidDataSets//regionsGeneTF_CovidADRelated_AtLeast2Sources_FullGRN.csv"), header = TRUE,check.names = FALSE)[,-1]

#adCovidGRNGenesListFilePath = "F://organizedAlzheimers//regionsGeneTF_CovidADRelated_AtLeast2Sources_FullGRN.csv"
library(EBSeq)
library(DESeq2)
library(edgeR)
library(plyr)

dim(virus_data_raw) # [1] 19967 genes  79 samples  # [1] before: 19874 genes   68 samples
rows = as.character(virus_data_raw[,1])
virus_data_raw = virus_data_raw[,-1]
rownames(virus_data_raw) = make.names(rows, unique=TRUE)
colnames(virus_data_raw)
virus_data_rawOLD = virus_data_raw
colnames(virus_data_raw)

head(row.names(virus_data_raw))
head(row.names(virus_data_raw))

virus_data_raw[1:5, 1:5]
colnames(virus_data_raw)

s1_rowsCOVID_ICU = c(1:50) # [1] 1st 50 rows are ICU
s1_rowsCOVID_ICU # 1 2 3

s2_rowsCOVID_nonICU = c(51:100) # next 50 rows are non-ICU
s2_rowsCOVID_nonICU # [1] 4 5 6
print(paste0(":) s1_rowsCOVID_ICU: ", length(s1_rowsCOVID_ICU)))
print(paste0(":) s2_rowsCOVID_nonICU: ", length(s2_rowsCOVID_nonICU)))


#####################################
# NORMALIZATION OF THE RAW DATA: 

RoundedInteger_virus_ALL <- apply(virus_data_raw, c(1, 2), function (x) {
  (round(x))
})

RoundedInteger_virus_ALL = RoundedInteger_virus_ALL


dim(RoundedInteger_virus_ALL)
virus_data_raw[1:5, 1:5]
RoundedInteger_virus_ALL[1:5, 1:5]



################################
# the main groups: 

groupsViruses = factor(c(rep("CovidICUSevere", length(s1_rowsCOVID_ICU)),
                         rep("CovidNonICU", length(s2_rowsCOVID_nonICU))))
groupsViruses

library("stringr")
allLabels = str_replace_all(colnames(virus_data_raw), "C", "")
allCovid = str_replace_all(colnames(virus_data_raw), "C", "")


covidColors  = c(rep("tomato", length(s1_rowsCOVID_ICU)),
                 rep("yellow", length(s2_rowsCOVID_nonICU)))#,

allColors = c(rep("tomato", length(s1_rowsCOVID_ICU)),
              rep("yellow", length(s2_rowsCOVID_nonICU)))

# PNG of the raw gene expression data comparative boxplots:

png(originalRawGeneExpressionDataBoxplotFilePathPNG, width = boxplotsPlotWidthPNG, height = boxplotsPlotHeightPNG) #original_selected_logQnorm.csv", header = TRUE, row.names = 1) #Advanced Bioinformatics\\Covid19\\April18\\medianNormalizedBoxplotCellLines.png", width = 2101, height = 1501)
par()
par(mfrow=c(2,1))
boxplot(virus_data_raw,  cex.main = 3, cex.lab = 3, ylim = c(0, 30), col = allColors, col.main = "red", main = "Raw Total Gene Read Counts for 100 Covid-19 Infected Human Samples (GSM4753022)",
        xaxt = NULL, xlab = "Different Covid-19 Infected Human Plasma and Leukocyte Samples", #ylab = "Total Gene Read Counts (GSM4753022)",
        cex.col = "blue", names = allLabels, las = 2)


plot(0,type='n',axes=FALSE,ann=FALSE) #, main = "Cell Line Samples Taken from a 79-Year Old Caucasian Woman", col.main = "blue")
legend("top", legend=c("Intensive Care Unit (ICU) Covid-19 Patients: Severe", 
                       "Non-ICU  Covid-19 Patients:  Not Severe"),
       
       col=unique(allColors), pch = 15, cex=3, pt.cex = 3.5)
dev.off()




# PDF of the raw gene expression data comparative boxplots:
pdf(originalRawGeneExpressionDataBoxplotFilePathPDF, width = boxplotsPlotWidthPDF, height = boxplotsPlotHeightPDF) #original_selected_logQnorm.csv", header = TRUE, row.names = 1) #Advanced Bioinformatics\\Covid19\\April18\\medianNormalizedBoxplotCellLines.png", width = 2101, height = 1501)
par()
par(mfrow=c(2,1))

boxplot(virus_data_raw,  cex.main = 3, cex.lab = 3, ylim = c(0, 30), col = allColors, col.main = "red", main = "Raw Total Gene Read Counts for 100 Covid-19 Infected Human Samples (GSM4753022)",
        xaxt = NULL, xlab = "Different Covid-19 Infected Human Plasma and Leukocyte Samples", #ylab = "Total Gene Read Counts (GSM4753022)",
        cex.col = "blue", names = allLabels, las = 2)


plot(0,type='n',axes=FALSE,ann=FALSE) #, main = "Cell Line Samples Taken from a 79-Year Old Caucasian Woman", col.main = "blue")
legend("top", legend=c("Intensive Care Unit (ICU) Covid-19 Patients: Severe", 
                       "Non-ICU  Covid-19 Patients:  Not Severe"),
       
       col=unique(allColors), pch = 15, cex=3, pt.cex = 3.5)
dev.off()












# Median Normalization of Gene Expression Data:
dds_ALL <- DESeqDataSetFromMatrix(RoundedInteger_virus_ALL, DataFrame(groupsViruses), ~ groupsViruses)
dds_ALL <- estimateSizeFactors(dds_ALL)
dds_ALL <- DESeq(dds_ALL)
rleNormalized_virus_ALL = counts(dds_ALL, normalized = TRUE)



# PNG of the median normalized gene expression data comparative boxplots:
png(nomalizedGeneExpressionDataBoxplotFilePathPNG, width = boxplotsPlotWidthPNG, height = boxplotsPlotHeightPNG) #original_selected_logQnorm.csv", header = TRUE, row.names = 1) #Advanced Bioinformatics\\Covid19\\April18\\medianNormalizedBoxplotCellLines.png", width = 2101, height = 1501)
par()
par(mfrow=c(2,1))
boxplot(rleNormalized_virus_ALL,  cex.main = 3, cex.lab = 3, ylim = c(0, 30), col = allColors, col.main = "red", main = "Median Normalized Total Gene Read Counts for 100 Covid-19 Infected Human Samples (GSM4753022)",
        xaxt = NULL, xlab = "Different Covid-19 Infected Human Plasma and Leukocyte Samples", #ylab = "Total Gene Read Counts (GSM4753022)",
        cex.col = "blue", names = allLabels, las = 2)


plot(0,type='n',axes=FALSE,ann=FALSE) #, main = "Cell Line Samples Taken from a 79-Year Old Caucasian Woman", col.main = "blue")
legend("top", legend=c("Intensive Care Unit (ICU) Covid-19 Patients: Severe", 
                       "Non-ICU  Covid-19 Patients:  Not Severe"),
       
       col=unique(allColors), pch = 15, cex=3, pt.cex = 3.5)

dev.off()


# PDF of the median normalized gene expression data comparative boxplots:
pdf(nomalizedGeneExpressionDataBoxplotFilePathPDF, width = boxplotsPlotWidthPDF, height = boxplotsPlotHeightPDF) #original_selected_logQnorm.csv", header = TRUE, row.names = 1) #Advanced Bioinformatics\\Covid19\\April18\\medianNormalizedBoxplotCellLines.png", width = 2101, height = 1501)
par()
par(mfrow=c(2,1))
boxplot(rleNormalized_virus_ALL,  cex.main = 3, cex.lab = 3, ylim = c(0, 30), col = allColors, col.main = "red", main = "Median Normalized Total Gene Read Counts for 100 Covid-19 Infected Human Samples (GSM4753022)",
        xaxt = NULL, xlab = "Different Covid-19 Infected Human Plasma and Leukocyte Samples", #ylab = "Total Gene Read Counts (GSM4753022)",
        cex.col = "blue", names = allLabels, las = 2)


plot(0,type='n',axes=FALSE,ann=FALSE) #, main = "Cell Line Samples Taken from a 79-Year Old Caucasian Woman", col.main = "blue")
legend("top", legend=c("Intensive Care Unit (ICU) Covid-19 Patients: Severe", 
                       "Non-ICU  Covid-19 Patients:  Not Severe"),
       
       col=unique(allColors), pch = 15, cex=3, pt.cex = 3.5)
dev.off()

################################################
severeDF  = rleNormalized_virus_ALL[,s1_rowsCOVID_ICU]
nonSevereDF  = rleNormalized_virus_ALL[, s2_rowsCOVID_nonICU]



################### DESEQ2 FUNCTIONS #####################
# this function calculates the group mean for each row in a dataSet
calculatingGroupMean <- function(dataSet){
  groupMeanVec = c()
  for (i in 1:nrow(dataSet)){
    row = as.numeric(dataSet[i,])
    rowMean = mean(row)
    groupMeanVec = c(groupMeanVec, rowMean)
  }
  return(groupMeanVec)
}


specificGroupMean <- function(groupName){
  groupMeanVec = c()
  if (groupName == "CovidICUSevere"){
    groupMeanVec = calculatingGroupMean(severeDF)
  }
  else if (groupName == "CovidNonICU"){
    groupMeanVec = calculatingGroupMean(nonSevereDF)
  }
  print(groupMeanVec)
  return(groupMeanVec)
}


# this function gets the column numbers
groupColumns <- function(groupName){
  columnNumbers = c()
  if (groupName == "CovidICUSevere"){
    columnNumbers = s1_rowsCOVID_ICU
  }
  else if (groupName == "CovidNonICU"){
    columnNumbers = s2_rowsCOVID_nonICU
  }
  #print(columnNumbers)
  return(columnNumbers)
  
  
}



dataFrameToAttach <- function(groupName){
  DFToAttach = NULL
  if (groupName == "CovidICUSevere"){
    dataFrameToAttach = severeDF
  }
  else {#if (groupName == "CovidNonICU"){
    dataFrameToAttach = nonSevereDF
  }
  #print(columnNumbers)
  return(dataFrameToAttach)
}

# this function performs deseq 2 and organizes and writes out (as a csv) the result of the differential gene analysis for the 2 health groups
twoGeneDifferentialAnalysis <- function(ddsDataFrame, res_2groups, grouping, baselineGroup, otherGroup, logFoldTermToAdd = 0.01){
  #twoGeneDifferentialAnalysis(rleNormalized_virus_ALL, res_SARSCOV2_Covid_NonIcu_vs_Icu, combiningGroups, "CovidNonICU", "CovidICUSevere"
  mcols(res_2groups, use.names=TRUE)
  resDataFrame_2Groups = as.data.frame(res_2groups)
  
  
  empiricalLFC = log2((rowMeans(ddsDataFrame[,groupColumns(otherGroup)]) + logFoldTermToAdd)/(rowMeans(ddsDataFrame[,groupColumns(baselineGroup)]) + logFoldTermToAdd))  # logFoldTermToAdd = 0.01 by default
  
  names(empiricalLFC) = rownames(res_2groups)
  
  resDataFrame_2Groups = cbind(specificGroupMean(baselineGroup), specificGroupMean(otherGroup), resDataFrame_2Groups, empiricalLFC)
  
  
  
  resDataFrame_2Groups = resDataFrame_2Groups[,-c(4:5)]  # removing the unneccessary columns 4 and 5
  
  colName1 = paste0(baselineGroup, "_ColumnsMean")
  colName2 = paste0(otherGroup, "_ColumnsMean")
  
  
  numGenes = nrow(resDataFrame_2Groups) # number of Genes :)
  starterVec = rep("Not Significant", numGenes)
  upRegulatedGeneIndexes = which(resDataFrame_2Groups$padj < 0.05 & resDataFrame_2Groups$empiricalLFC > 0)
  downRegulatedGeneIndexes = which(resDataFrame_2Groups$padj < 0.05 & resDataFrame_2Groups$empiricalLFC < 0)
  starterVec[upRegulatedGeneIndexes] = paste("UP-Regulated from: ", baselineGroup, " to ", otherGroup)
  starterVec[downRegulatedGeneIndexes] = paste("DOWN-Regulated from: ", baselineGroup, " to ", otherGroup)
  
  
  resDataFrame_2Groups = cbind(resDataFrame_2Groups, starterVec)
  colnames(resDataFrame_2Groups) <- c(colName1, colName2, "baseMean", "stat", "pvalue", "padj", "empiricalLFC", "regulationStatus")
  
  resDataFrame_2Groups <- cbind(resDataFrame_2Groups, dataFrameToAttach(baselineGroup))  # adding the baseline group columns
  resDataFrame_2Groups <- cbind(resDataFrame_2Groups, dataFrameToAttach(otherGroup))  # adding the other group columns
  
  fileName = paste0("DE_", baselineGroup, "baseLine_vs_", otherGroup, ".csv")
  
  # writing the output data frame to a csv file
  FilePath = paste0(dataSetOutputPaths, fileName)
  print(FilePath)
  FilePath  # [1] "C:\\Users\\sk792\\Documents\\Spring 2020```````.csv"
  write.csv(resDataFrame_2Groups, FilePath)
  return(resDataFrame_2Groups)
}

################### ################### ################### 

res_SARSCOV2_Covid_NonIcu_vs_Icu <- results(dds_ALL, contrast=c("groupsViruses","CovidNonICU","CovidICUSevere"))
length(which(res_SARSCOV2_Covid_NonIcu_vs_Icu$padj < 0.05)) # 
combiningGroups = groupsViruses
GenesComparison_SARSCOV2_Covid_NonIcu_vs_Icu = twoGeneDifferentialAnalysis(rleNormalized_virus_ALL, res_SARSCOV2_Covid_NonIcu_vs_Icu, combiningGroups, "CovidNonICU", "CovidICUSevere")
#"C://Users//sk792//Documents//Spring 2020//Advanced Bioinformatics//Covid19//April15//GroupDEComparison_s1_COVID_NHBELungCells_NormalbaseLine_vs_s1_COVID_NHBELungCells_Infected_MOI2.csv"







#######################################################################
########### 
region = c("LTL", "Hippocampus", "DLPFC")

###################################################################################
# LTL:
num = 1

myList = read.csv(adCovidGRNGenesListFilePath header = TRUE)[,-1]

myList = myList[which(myList$region == region[num]),]
dim(myList)
volcanoPlotOutputFileName
#outputPath = paste0("F://organizedAlzheimers//differentialExpressionAnalysisCovidAD//volcanoPlotADCovidICUvsNon-ICU_", region[num], "_Updated.png")

outputPath = paste0(volcanoPlotOutputFileName, region[num], ".png")

head(myList)
myList = as.vector(myList[,1])
myList
require(ggplot2)
library("ggrepel")
library("ggplot2")

resDataFrame = GenesComparison_SARSCOV2_Covid_NonIcu_vs_Icu #GenesComparison_SARSCOV2_Covid_NonIcu_vs_Icu #res_SARSCOV2_Covid_NonIcu_vs_Icu
bodyRegion = region[num]
print(paste0(":) please note the bodyRegion is: ", bodyRegion))
if(region[num] == "Hippocampus"){
  groupComparisonString = "Not Severe (Non-ICU) Versus Severe (ICU) Covid-19 Human Samples in the Hippocampus Brain Region"
} else if (region[num] == "LTL"){
  groupComparisonString = "Not Severe (Non-ICU) Versus Severe (ICU) Covid-19 Human Samples in the LTL Brain Region"
  
}else if (region[num] == "DLPFC"){
  groupComparisonString = "Not Severe (Non-ICU) Versus Severe (ICU) Covid-19 Human Samplesin the DLPFC Brain Region"
  
}

####################################

downRegColor = "green4"
noSigChangeColor = "khaki" #yellow" #"grey"#rgb(0, 0, 255, max = 255, alpha = 125, names = "grey") # "grey"
upRegColor = "red"
log2FoldChangeCutOff = 0
numOfGenesToLabel = 60
# Please note that this function by Saniya creates a Volcano Plot object for a given pairwise comparison of Adenocarcinoma Patients :)

####################################

resDataFrame = resDataFrame[-which(is.na(resDataFrame$padj)),] # removing NAs for padj
dim(resDataFrame)
log2FoldChangeCutOff = 0  # please note that this is the default we use
adjustedPvalueCutoff = 0.05 # please note that this is the default we use

log2col  =  grep("empirical", colnames(resDataFrame))  
log2col # 7 is the column number with the log_2_fold changes

pAdjcol  =  grep("padj", colnames(resDataFrame))  
pAdjcol

pAdj_values =  resDataFrame[,pAdjcol]

log2col_values =  resDataFrame[,log2col]
diff_df = resDataFrame

intersect(myList, row.names(diff_df))#$external_gene_name)
matches = match(myList, row.names(diff_df))#diff_df$external_gene_name)
matches = matches[!is.na(matches)]
matches
length(matches) # 126
up = which(log2col_values > log2FoldChangeCutOff & pAdj_values < adjustedPvalueCutoff)
down = which(log2col_values < log2FoldChangeCutOff & pAdj_values < adjustedPvalueCutoff)
print(paste0("up genes: ", length(up))) # 2,505
print(paste0("down genes: ", length(down))) # 2,580
print(paste0("Differentially Expressed Gene Totals: ", 
             length(unique(c(up, down)))))

notSignificant = which(pAdj_values >= adjustedPvalueCutoff)
notSignificant
notSignificantMatches = intersect(matches, notSignificant)

othersNotSignificant = setdiff(notSignificant, matches)
print(paste0(":) othersNotSignificant: ", length(othersNotSignificant)))

print(paste0(":) notSignificantMatches: ", length(notSignificantMatches)))

upReg_Matches = intersect(matches, up)
downReg_Matches = intersect(matches, down)

upRegMatchesLTL = row.names(diff_df)[upReg_Matches]
downRegMatchesLTL = row.names(diff_df)[downReg_Matches]
print(paste0("Differentially Expressed Gene Matches ", bodyRegion, ": ", 
             length(unique(c(upReg_Matches, downReg_Matches)))))
print(paste0("upReg_Matches AD/Covid genes: ", length(upReg_Matches))) # 45
print(paste0("downReg_Matches AD/Covid genes: ", length(downReg_Matches))) # 19

newGenes = c(up, down)
print(length(newGenes))
otherSignificant = setdiff(newGenes, matches)
otherSignificant
print(paste0("otherSignificant genes: ", length(otherSignificant))) # 19

length(upReg_Matches) + length(downReg_Matches)  + length(otherSignificant) + length(notSignificantMatches) + length(othersNotSignificant)
dim(diff_df)# 12927  

regulationStatusVec = c()
#colorVector = c()
groupVector = c()
for (i in 1:length(log2col_values)){
  if (i %in%  upReg_Matches)
  {
    regulationStatusVec = c(regulationStatusVec, "Up Regulated Matches")
    groupVector = c(groupVector, "Up_Regulated_Match")
    
  }
  else if (i %in% downReg_Matches)
  {
    regulationStatusVec = c(regulationStatusVec, "Down Regulated Matches")
    groupVector = c(groupVector, "Down_Regulated_Match")
    
  }  else if (i %in% otherSignificant)
  {
    regulationStatusVec = c(regulationStatusVec, "Other Significant Genes")
    groupVector = c(groupVector, "Other_Significant_Gene")
    
  } 
  
  else if (i %in% notSignificantMatches)
  {
    regulationStatusVec = c(regulationStatusVec, "Not Significant Matches")
    groupVector = c(groupVector, "Not_Significant_Match")
    
  } else if (i %in% othersNotSignificant) {
    regulationStatusVec = c(regulationStatusVec, "Not Significant Genes")
    groupVector = c(groupVector, "Not_Significant_Gene")
    
  }
  else
  {
    regulationStatusVec = c(regulationStatusVec, "No Significant Changes")
    
  }
  
}

print(paste0("Please note that this is the Breakdown Table of Down-Regulated, No Significant Changes, and Up-Regulated Genes for : ", groupComparisonString))
print(table(regulationStatusVec))

pAdjcol  =  grep("padj", colnames(resDataFrame))  
pAdjcol # 6 is the column number with the adjusted p-values
pAdj_values =  resDataFrame[,pAdjcol]


diff_df["group"] <- "No_Significant_Changes"

# for our plot, we want to highlight 
# FDR < 0.05 (significance level)
# Fold Change > 1.5

# change the grouping for the entries with significance but not a large enough Fold change
diff_df$group = groupVector

diff_df["external_gene_name"] = row.names(diff_df)

head(diff_df)
# Add gene labels for all of the top genes we found
# here we are creating an empty list, and filling it with entries for each row in the dataframe
# each list entry is another list with named items that will be used by Plot.ly

diff_df["Significance"] <- diff_df["group"] 
levels(diff_df$Significance)[levels(diff_df$Significance) == "Not_Significant_Gene"] <- "Not Significant Gene"
levels(diff_df$Significance)[levels(diff_df$Significance) == "Other_Significant_Gene"] <- "Other Significant Gene"
levels(diff_df$Significance)[levels(diff_df$Significance) == "Down_Regulated_Match"] <- "Down Regulated Match"
levels(diff_df$Significance)[levels(diff_df$Significance) == "Up_Regulated_Match"] <- "Up Regulated Match"
levels(diff_df$Significance)[levels(diff_df$Significance) == "Not_Significant_Match"] <- "Not Significant Match"


diff_df$Significance[which(diff_df$Significance == "Not_Significant_Gene")] <- "No Significant Change"
diff_df$Significance[which(diff_df$Significance == "Other_Significant_Gene")] <- "Other Significant Gene"
diff_df$Significance[which(diff_df$Significance == "Down_Regulated_Match")] <- "Down Regulated Match"
diff_df$Significance[which(diff_df$Significance == "Up_Regulated_Match")] <- "Up Regulated Match"
diff_df$Significance[which(diff_df$Significance == "Not_Significant_Match")] <- "Not Significant Match"  

head(diff_df)

diff_df["AlphaValues"] = rep(1, nrow(diff_df))
diff_df$AlphaValues[which(diff_df$group == "Down_Regulated_Match")] <- 1
diff_df$AlphaValues[which(diff_df$group == "Up_Regulated_Match")] <- 1
diff_df$AlphaValues[which(diff_df$group == "Not_Significant_Match")] <- 1
diff_df$AlphaValues[which(diff_df$group == "Not_Significant_Gene")] <- 1


table(diff_df$AlphaValues)
# 0.01     1 
# 12801   126
# Find and label the top peaks..
numOfGenesToLabel = 60
group_numGenes = round(numOfGenesToLabel/3)  # there are 3 groups below that we label:
print(paste0("Please note the actual total number of labeled genes for this ", groupComparisonString, " Volcano Plot is: ", (3*group_numGenes)))


#these are the 3 groups below that we will label:
diff_df = diff_df[with(diff_df, order(-diff_df$empiricalLFC, diff_df$padj)),]

genesToSelectFromEachCategory = 20

upMatchDF = diff_df[which(diff_df$group == "Up_Regulated_Match"),]
upMatchDF = upMatchDF[order(upMatchDF$padj),]


upMatchDFpart2 = diff_df[which(diff_df$group == "Up_Regulated_Match"),]

upMatchDFpart2 = upMatchDF[order(upMatchDF$empiricalLFC, decreasing = TRUE),]

upMatchDF = upMatchDF[1:genesToSelectFromEachCategory,]


upMatchDFpart2 = upMatchDFpart2[1:genesToSelectFromEachCategory,]


downMatchDF = diff_df[which(diff_df$group == "Down_Regulated_Match"),]
downMatchDF = downMatchDF[order(downMatchDF$padj),]#1:genesToSelectFromEachCategory]    
downMatchDF = downMatchDF[1:genesToSelectFromEachCategory,]


nonSigMatchDF = diff_df[which(diff_df$group == "Not_Significant_Match"),]
nonSigMatchDF_part1 = nonSigMatchDF[order(nonSigMatchDF$empiricalLFC),] #1:(genesToSelectFromEachCategory/2)]    
nonSigMatchDF_part2 = nonSigMatchDF[order(nonSigMatchDF$empiricalLFC, decreasing = TRUE), ] #decreasing = TRUE),1:(genesToSelectFromEachCategory/2)]    
nonSigMatchDF_part1 = nonSigMatchDF_part1[1:genesToSelectFromEachCategory,]
nonSigMatchDF_part2 = nonSigMatchDF_part2[1:genesToSelectFromEachCategory,]

top_peaks <- upMatchDF #diff_df[which(diff_df$Significance == "Down Regulated Match"),][1:newGenes,] #diff_df[with(diff_df, order(diff_df$empiricalLFC, diff_df$padj)),][1:group_numGenes,]  # positive logs
top_peaks <- rbind(top_peaks, upMatchDFpart2) #diff_df[which(diff_df$Significance != "Up Regulated Match"),][1:newGenes,]) #diff_df[with(diff_df, order(diff_df$empiricalLFC, diff_df$padj)),][1:group_numGenes,]  # positive logs
top_peaks <- rbind(top_peaks, downMatchDF) #diff_df[which(diff_df$Significance != "Up Regulated Match"),][1:newGenes,]) #diff_df[with(diff_df, order(diff_df$empiricalLFC, diff_df$padj)),][1:group_numGenes,]  # positive logs
top_peaks <- rbind(top_peaks, nonSigMatchDF_part1)
top_peaks <- rbind(top_peaks, nonSigMatchDF_part2)
top_peaks = unique(top_peaks)

top_peaks = diff_df[grep("Up_Regulated_Match", diff_df$group),]
top_peaks = rbind(top_peaks, diff_df[grep("Down_Regulated_Match", diff_df$group),])
top_peaks = unique(top_peaks)

unique(diff_df$group)
newGenes =group_numGenes/2

nonMatches = diff_df[-grep("Match", diff_df$group),]
matches = diff_df[grep("Match", diff_df$group),]
dim(nonMatches) # 12801   112
dim(matches) # 126   112

diff_df = nonMatches
# so matches are plotted after and come on top more
diff_df = rbind(diff_df, matches)
head(top_peaks)

volcanoTitle = paste0("Volcano Plot: ", groupComparisonString) 
new_demo = ggplot(diff_df, aes(x = empiricalLFC, y = -log10(padj))) + geom_point(aes(color = Significance, alpha = AlphaValues, size = 2)) + scale_color_manual(values = c(downRegColor, noSigChangeColor, "khaki", "khaki", upRegColor)) + theme_bw(base_size = 12) + theme(legend.position = "bottom") + geom_text_repel(data = top_peaks, aes(label = external_gene_name), size = 5, box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"))
new_demo = new_demo + ylim(0, 10) + xlim(-3, 3.5)

mynamestheme <- theme(plot.title = element_text(family = "Helvetica", face = "bold", size = (25), hjust = 0.5), legend.title = element_text(colour = "steelblue", face = "bold.italic", family = "Helvetica", size = (15)), 
                      legend.text = element_text(face = "italic", colour = "steelblue4", family = "Helvetica", size = (17)),
                      axis.title = element_text(family = "Helvetica", size = (21), colour = "steelblue4"),
                      axis.text = element_text(family = "Courier", colour = "cornflowerblue", size = (19)))
final_plot = new_demo + mynamestheme + labs(title = volcanoTitle, x = expression(Log[2]~(Fold~Change)), y = expression(-Log[10]~(Adjusted~P-Value))) + guides(colour = guide_legend(override.aes = list(size = 10)))

final_plot = final_plot + geom_hline(yintercept=1, linetype="dotted", size = 1, colour = "steelblue4")


final_plot


png(outputPath, width = volcanoPlotWidth, height = volcanoPlotHeight)
print(final_plot)
dev.off()


#######################################################################################

###################################################################################
# Hippocampus:
num = 2
# labeling only a few points since many are significant for Hippocampus

myList = read.csv(adCovidGRNGenesListFilePath header = TRUE)[,-1]

myList = myList[which(myList$region == region[num]),]
dim(myList)

outputPath = paste0(volcanoPlotOutputFileName, region[num], ".png")head(myList)
myList = as.vector(myList[,1])
myList


resDataFrame = GenesComparison_SARSCOV2_Covid_NonIcu_vs_Icu #GenesComparison_SARSCOV2_Covid_NonIcu_vs_Icu #res_SARSCOV2_Covid_NonIcu_vs_Icu
bodyRegion = region[num]
print(paste0(":) please note the bodyRegion is: ", bodyRegion))
if(region[num] == "Hippocampus"){
  groupComparisonString = "Not Severe (Non-ICU) Versus Severe (ICU) Covid-19 Human Samples in the Hippocampus Brain Region"
} else if (region[num] == "LTL"){
  groupComparisonString = "Not Severe (Non-ICU) Versus Severe (ICU) Covid-19 Human Samples in the LTL Brain Region"
  
}else if (region[num] == "DLPFC"){
  groupComparisonString = "Not Severe (Non-ICU) Versus Severe (ICU) Covid-19 Human Samplesin the DLPFC Brain Region"
  
}

####################################

downRegColor = "green4"
noSigChangeColor = "khaki"
upRegColor = "red"
log2FoldChangeCutOff = 0
numOfGenesToLabel = 60
# Please note that this function by Saniya creates a Volcano Plot object for a given pairwise comparison of Adenocarcinoma Patients :)

####################################

resDataFrame = resDataFrame[-which(is.na(resDataFrame$padj)),] # removing NAs for padj
dim(resDataFrame)
log2FoldChangeCutOff = 0  # please note that this is the default we use
adjustedPvalueCutoff = 0.05 # please note that this is the default we use

log2col  =  grep("empirical", colnames(resDataFrame))  
log2col # 7 is the column number with the log_2_fold changes

pAdjcol  =  grep("padj", colnames(resDataFrame))  
pAdjcol

pAdj_values =  resDataFrame[,pAdjcol]

log2col_values =  resDataFrame[,log2col]
diff_df = resDataFrame

intersect(myList, row.names(diff_df))#$external_gene_name)
matches = match(myList, row.names(diff_df))#diff_df$external_gene_name)
matches = matches[!is.na(matches)]
matches
length(matches) # 126
up = which(log2col_values > log2FoldChangeCutOff & pAdj_values < adjustedPvalueCutoff)
down = which(log2col_values < log2FoldChangeCutOff & pAdj_values < adjustedPvalueCutoff)
print(paste0("up genes: ", length(up))) # 2,505
print(paste0("down genes: ", length(down))) # 2,580
print(paste0("Differentially Expressed Gene Totals: ", 
             length(unique(c(up, down)))))

notSignificant = which(pAdj_values >= adjustedPvalueCutoff)
notSignificant
notSignificantMatches = intersect(matches, notSignificant)

othersNotSignificant = setdiff(notSignificant, matches)
print(paste0(":) othersNotSignificant: ", length(othersNotSignificant)))

print(paste0(":) notSignificantMatches: ", length(notSignificantMatches)))

upReg_Matches = intersect(matches, up)
downReg_Matches = intersect(matches, down)


upRegMatchesHipp= row.names(diff_df)[upReg_Matches]
downRegMatchesHipp = row.names(diff_df)[downReg_Matches]

print(paste0("Differentially Expressed Gene Matches ", bodyRegion, ": ", 
             length(unique(c(upReg_Matches, downReg_Matches)))))
print(paste0("upReg_Matches AD/Covid genes: ", length(upReg_Matches))) # 45
print(paste0("downReg_Matches AD/Covid genes: ", length(downReg_Matches))) # 19

newGenes = c(up, down)
print(length(newGenes))
otherSignificant = setdiff(newGenes, matches)
otherSignificant
print(paste0("otherSignificant genes: ", length(otherSignificant))) # 19

length(upReg_Matches) + length(downReg_Matches)  + length(otherSignificant) + length(notSignificantMatches) + length(othersNotSignificant)
dim(diff_df)# 12927  

regulationStatusVec = c()
#colorVector = c()
groupVector = c()
for (i in 1:length(log2col_values)){
  if (i %in%  upReg_Matches)
  {
    regulationStatusVec = c(regulationStatusVec, "Up Regulated Matches")
    #colorVector = c(colorVector, "red")
    groupVector = c(groupVector, "Up_Regulated_Match")
    
  }
  else if (i %in% downReg_Matches)
  {
    regulationStatusVec = c(regulationStatusVec, "Down Regulated Matches")
    #colorVector = c(colorVector, "dark green")
    groupVector = c(groupVector, "Down_Regulated_Match")
    
  }  else if (i %in% otherSignificant)
  {
    regulationStatusVec = c(regulationStatusVec, "Other Significant Genes")
    groupVector = c(groupVector, "Other_Significant_Gene")
    
  } 
  
  else if (i %in% notSignificantMatches)
  {
    regulationStatusVec = c(regulationStatusVec, "Not Significant Matches")
    #colorVector = c(colorVector, "dark green")
    groupVector = c(groupVector, "Not_Significant_Match")
    
  } else if (i %in% othersNotSignificant) {
    regulationStatusVec = c(regulationStatusVec, "Not Significant Genes")
    groupVector = c(groupVector, "Not_Significant_Gene")
    
  }
  else
  {
    regulationStatusVec = c(regulationStatusVec, "No Significant Changes")
    
    # colorVector = c(colorVector, "light gray")
  }
  
}

print(paste0("Please note that this is the Breakdown Table of Down-Regulated, No Significant Changes, and Up-Regulated Genes for : ", groupComparisonString))
print(table(regulationStatusVec))

pAdjcol  =  grep("padj", colnames(resDataFrame))  
pAdjcol # 6 is the column number with the adjusted p-values
pAdj_values =  resDataFrame[,pAdjcol]


diff_df["group"] <- "No_Significant_Changes"

# for our plot, we want to highlight 
# FDR < 0.05 (significance level)
# Fold Change > 1.5

# change the grouping for the entries with significance but not a large enough Fold change
diff_df$group = groupVector

diff_df["external_gene_name"] = row.names(diff_df)

head(diff_df)
# Add gene labels for all of the top genes we found
# here we are creating an empty list, and filling it with entries for each row in the dataframe
# each list entry is another list with named items that will be used by Plot.ly

diff_df["Significance"] <- diff_df["group"] 
levels(diff_df$Significance)[levels(diff_df$Significance) == "Not_Significant_Gene"] <- "Not Significant Gene"
levels(diff_df$Significance)[levels(diff_df$Significance) == "Other_Significant_Gene"] <- "Other Significant Gene"
levels(diff_df$Significance)[levels(diff_df$Significance) == "Down_Regulated_Match"] <- "Down Regulated Match"
levels(diff_df$Significance)[levels(diff_df$Significance) == "Up_Regulated_Match"] <- "Up Regulated Match"
levels(diff_df$Significance)[levels(diff_df$Significance) == "Not_Significant_Match"] <- "Not Significant Match"


diff_df$Significance[which(diff_df$Significance == "Not_Significant_Gene")] <- "No Significant Change"
diff_df$Significance[which(diff_df$Significance == "Other_Significant_Gene")] <- "Other Significant Gene"
diff_df$Significance[which(diff_df$Significance == "Down_Regulated_Match")] <- "Down Regulated Match"
diff_df$Significance[which(diff_df$Significance == "Up_Regulated_Match")] <- "Up Regulated Match"
diff_df$Significance[which(diff_df$Significance == "Not_Significant_Match")] <- "Not Significant Match"  

head(diff_df)

diff_df["AlphaValues"] = rep(1, nrow(diff_df))
diff_df$AlphaValues[which(diff_df$group == "Down_Regulated_Match")] <- 1
diff_df$AlphaValues[which(diff_df$group == "Up_Regulated_Match")] <- 1
diff_df$AlphaValues[which(diff_df$group == "Not_Significant_Match")] <- 1
diff_df$AlphaValues[which(diff_df$group == "Not_Significant_Gene")] <- 1


table(diff_df$AlphaValues)
# 0.01     1 
# 12801   126
# Find and label the top peaks..
numOfGenesToLabel = 60
group_numGenes = round(numOfGenesToLabel/3)  # there are 3 groups below that we label:
print(paste0("Please note the actual total number of labeled genes for this ", groupComparisonString, " Volcano Plot is: ", (3*group_numGenes)))


#these are the 3 groups below that we will label:
diff_df = diff_df[with(diff_df, order(-diff_df$empiricalLFC, diff_df$padj)),]

genesToSelectFromEachCategory = 72

upMatchDF = diff_df[which(diff_df$group == "Up_Regulated_Match"),]
upMatchDF = upMatchDF[order(upMatchDF$padj),]


upMatchDFpart2 = diff_df[which(diff_df$group == "Up_Regulated_Match"),]

upMatchDFpart2 = upMatchDF[order(upMatchDF$empiricalLFC, decreasing = TRUE),]

upMatchDF = upMatchDF[1:genesToSelectFromEachCategory,]


upMatchDFpart2 = upMatchDFpart2[1:genesToSelectFromEachCategory,]


downMatchDF = diff_df[which(diff_df$group == "Down_Regulated_Match"),]
downMatchDF = downMatchDF[order(downMatchDF$padj),]#1:genesToSelectFromEachCategory]    
downMatchDF = downMatchDF[1:genesToSelectFromEachCategory,]


nonSigMatchDF = diff_df[which(diff_df$group == "Not_Significant_Match"),]
nonSigMatchDF_part1 = nonSigMatchDF[order(nonSigMatchDF$empiricalLFC),] #1:(genesToSelectFromEachCategory/2)]    
nonSigMatchDF_part2 = nonSigMatchDF[order(nonSigMatchDF$empiricalLFC, decreasing = TRUE), ] #decreasing = TRUE),1:(genesToSelectFromEachCategory/2)]    
nonSigMatchDF_part1 = nonSigMatchDF_part1[1:genesToSelectFromEachCategory,]
nonSigMatchDF_part2 = nonSigMatchDF_part2[1:genesToSelectFromEachCategory,]

top_peaks <- upMatchDF #diff_df[which(diff_df$Significance == "Down Regulated Match"),][1:newGenes,] #diff_df[with(diff_df, order(diff_df$empiricalLFC, diff_df$padj)),][1:group_numGenes,]  # positive logs
top_peaks <- rbind(top_peaks, upMatchDFpart2) #diff_df[which(diff_df$Significance != "Up Regulated Match"),][1:newGenes,]) #diff_df[with(diff_df, order(diff_df$empiricalLFC, diff_df$padj)),][1:group_numGenes,]  # positive logs
top_peaks <- rbind(top_peaks, downMatchDF) #diff_df[which(diff_df$Significance != "Up Regulated Match"),][1:newGenes,]) #diff_df[with(diff_df, order(diff_df$empiricalLFC, diff_df$padj)),][1:group_numGenes,]  # positive logs
top_peaks <- rbind(top_peaks, nonSigMatchDF_part1)
top_peaks <- rbind(top_peaks, nonSigMatchDF_part2)
top_peaks = unique(top_peaks)
dim(top_peaks)
otherTop_peaks = diff_df[grep("Up_Regulated_Match", diff_df$group),]
otherTop_peaks = rbind(otherTop_peaks, diff_df[grep("Down_Regulated_Match", diff_df$group),])
otherTop_peaks = unique(otherTop_peaks)

top_peaks = top_peaks[which(row.names(top_peaks) %in% row.names(otherTop_peaks)),]
dim(top_peaks)

unique(diff_df$group)
newGenes =group_numGenes/2

nonMatches = diff_df[-grep("Match", diff_df$group),]
matches = diff_df[grep("Match", diff_df$group),]
dim(nonMatches) # 12801   112
dim(matches) # 126   112

diff_df = nonMatches
# so matches are plotted after and come on top more
diff_df = rbind(diff_df, matches)
head(top_peaks)

volcanoTitle = paste0("Volcano Plot: ", groupComparisonString) 
new_demo = ggplot(diff_df, aes(x = empiricalLFC, y = -log10(padj))) + geom_point(aes(color = Significance, alpha = AlphaValues, size = 2)) + scale_color_manual(values = c(downRegColor, noSigChangeColor, "khaki", "khaki", upRegColor)) + theme_bw(base_size = 12) + theme(legend.position = "bottom") + geom_text_repel(data = top_peaks, aes(label = external_gene_name), size = 5, box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"))
new_demo = new_demo + ylim(0, 10) + xlim(-3, 3.5)

mynamestheme <- theme(plot.title = element_text(family = "Helvetica", face = "bold", size = (25), hjust = 0.5), legend.title = element_text(colour = "steelblue", face = "bold.italic", family = "Helvetica", size = (15)), 
                      legend.text = element_text(face = "italic", colour = "steelblue4", family = "Helvetica", size = (17)),
                      axis.title = element_text(family = "Helvetica", size = (21), colour = "steelblue4"),
                      axis.text = element_text(family = "Courier", colour = "cornflowerblue", size = (19)))
final_plot = new_demo + mynamestheme + labs(title = volcanoTitle, x = expression(Log[2]~(Fold~Change)), y = expression(-Log[10]~(Adjusted~P-Value))) + guides(colour = guide_legend(override.aes = list(size = 10)))

final_plot = final_plot + geom_hline(yintercept=1, linetype="dotted", size = 1, colour = "steelblue4")


final_plot


png(outputPath, width = volcanoPlotWidth, height = volcanoPlotHeight)
print(final_plot)
dev.off()




###################################################################################
# DLPFC:
num = 3

myList = read.csv(adCovidGRNGenesListFilePath header = TRUE)[,-1]

myList = myList[which(myList$region == region[num]),]
dim(myList)
outputPath = paste0(volcanoPlotOutputFileName, region[num], ".png")
head(myList)
myList = as.vector(myList[,1])
myList


resDataFrame = GenesComparison_SARSCOV2_Covid_NonIcu_vs_Icu 
bodyRegion = region[num]
print(paste0(":) please note the bodyRegion is: ", bodyRegion))
if(region[num] == "Hippocampus"){
  groupComparisonString = "Not Severe (Non-ICU) Versus Severe (ICU) Covid-19 Human Samples in the Hippocampus Brain Region"
} else if (region[num] == "LTL"){
  groupComparisonString = "Not Severe (Non-ICU) Versus Severe (ICU) Covid-19 Human Samples in the LTL Brain Region"
  
}else if (region[num] == "DLPFC"){
  groupComparisonString = "Not Severe (Non-ICU) Versus Severe (ICU) Covid-19 Human Samplesin the DLPFC Brain Region"
  
}

####################################

downRegColor = "green4"
noSigChangeColor = "khaki" #yellow"
upRegColor = "red"
log2FoldChangeCutOff = 0
numOfGenesToLabel = 60
# Please note that this function by Saniya creates a Volcano Plot object for a given pairwise comparison of Adenocarcinoma Patients :)

####################################

resDataFrame = resDataFrame[-which(is.na(resDataFrame$padj)),] # removing NAs for padj
dim(resDataFrame)
log2FoldChangeCutOff = 0  # please note that this is the default we use
adjustedPvalueCutoff = 0.05 # please note that this is the default we use

log2col  =  grep("empirical", colnames(resDataFrame))  
log2col # 7 is the column number with the log_2_fold changes

pAdjcol  =  grep("padj", colnames(resDataFrame))  
pAdjcol

pAdj_values =  resDataFrame[,pAdjcol]

log2col_values =  resDataFrame[,log2col]
diff_df = resDataFrame

intersect(myList, row.names(diff_df))#$external_gene_name)
matches = match(myList, row.names(diff_df))#diff_df$external_gene_name)
matches = matches[!is.na(matches)]
matches
length(matches) # 126
up = which(log2col_values > log2FoldChangeCutOff & pAdj_values < adjustedPvalueCutoff)
down = which(log2col_values < log2FoldChangeCutOff & pAdj_values < adjustedPvalueCutoff)
print(paste0("up genes: ", length(up))) # 2,505
print(paste0("down genes: ", length(down))) # 2,580
print(paste0("Differentially Expressed Gene Totals: ", 
             length(unique(c(up, down)))))

notSignificant = which(pAdj_values >= adjustedPvalueCutoff)
notSignificant
notSignificantMatches = intersect(matches, notSignificant)

othersNotSignificant = setdiff(notSignificant, matches)
print(paste0(":) othersNotSignificant: ", length(othersNotSignificant)))

print(paste0(":) notSignificantMatches: ", length(notSignificantMatches)))

upReg_Matches = intersect(matches, up)
downReg_Matches = intersect(matches, down)


upRegMatchesDLPFC = row.names(diff_df)[upReg_Matches]
downRegMatchesDLPFC = row.names(diff_df)[downReg_Matches]


print(paste0("Differentially Expressed Gene Matches ", bodyRegion, ": ", 
             length(unique(c(upReg_Matches, downReg_Matches)))))
print(paste0("upReg_Matches AD/Covid genes: ", length(upReg_Matches))) # 45
print(paste0("downReg_Matches AD/Covid genes: ", length(downReg_Matches))) # 19

newGenes = c(up, down)
print(length(newGenes))
otherSignificant = setdiff(newGenes, matches)
otherSignificant
print(paste0("otherSignificant genes: ", length(otherSignificant))) # 19

length(upReg_Matches) + length(downReg_Matches)  + length(otherSignificant) + length(notSignificantMatches) + length(othersNotSignificant)
dim(diff_df)

regulationStatusVec = c()
groupVector = c()
for (i in 1:length(log2col_values)){
  if (i %in%  upReg_Matches)
  {
    regulationStatusVec = c(regulationStatusVec, "Up Regulated Matches")
    groupVector = c(groupVector, "Up_Regulated_Match")
    
  }
  else if (i %in% downReg_Matches)
  {
    regulationStatusVec = c(regulationStatusVec, "Down Regulated Matches")
    groupVector = c(groupVector, "Down_Regulated_Match")
    
  }  else if (i %in% otherSignificant)
  {
    regulationStatusVec = c(regulationStatusVec, "Other Significant Genes")
    groupVector = c(groupVector, "Other_Significant_Gene")
    
  } 
  
  else if (i %in% notSignificantMatches)
  {
    regulationStatusVec = c(regulationStatusVec, "Not Significant Matches")
    groupVector = c(groupVector, "Not_Significant_Match")
    
  } else if (i %in% othersNotSignificant) {
    regulationStatusVec = c(regulationStatusVec, "Not Significant Genes")
    groupVector = c(groupVector, "Not_Significant_Gene")
    
  }
  else
  {
    regulationStatusVec = c(regulationStatusVec, "No Significant Changes")
    
  }
  
}

print(paste0("Please note that this is the Breakdown Table of Down-Regulated, No Significant Changes, and Up-Regulated Genes for : ", groupComparisonString))
print(table(regulationStatusVec))

pAdjcol  =  grep("padj", colnames(resDataFrame))  
pAdjcol # 6 is the column number with the adjusted p-values
pAdj_values =  resDataFrame[,pAdjcol]


diff_df["group"] <- "No_Significant_Changes"

# for our plot, we want to highlight 
# FDR < 0.05 (significance level)
# Fold Change > 1.5

# change the grouping for the entries with significance but not a large enough Fold change
diff_df$group = groupVector

diff_df["external_gene_name"] = row.names(diff_df)

head(diff_df)
# Add gene labels for all of the top genes we found
# here we are creating an empty list, and filling it with entries for each row in the dataframe
# each list entry is another list with named items that will be used by Plot.ly

diff_df["Significance"] <- diff_df["group"] 
levels(diff_df$Significance)[levels(diff_df$Significance) == "Not_Significant_Gene"] <- "Not Significant Gene"
levels(diff_df$Significance)[levels(diff_df$Significance) == "Other_Significant_Gene"] <- "Other Significant Gene"
levels(diff_df$Significance)[levels(diff_df$Significance) == "Down_Regulated_Match"] <- "Down Regulated Match"
levels(diff_df$Significance)[levels(diff_df$Significance) == "Up_Regulated_Match"] <- "Up Regulated Match"
levels(diff_df$Significance)[levels(diff_df$Significance) == "Not_Significant_Match"] <- "Not Significant Match"


diff_df$Significance[which(diff_df$Significance == "Not_Significant_Gene")] <- "No Significant Change"
diff_df$Significance[which(diff_df$Significance == "Other_Significant_Gene")] <- "Other Significant Gene"
diff_df$Significance[which(diff_df$Significance == "Down_Regulated_Match")] <- "Down Regulated Match"
diff_df$Significance[which(diff_df$Significance == "Up_Regulated_Match")] <- "Up Regulated Match"
diff_df$Significance[which(diff_df$Significance == "Not_Significant_Match")] <- "Not Significant Match"  

head(diff_df)

diff_df["AlphaValues"] = rep(1, nrow(diff_df))
diff_df$AlphaValues[which(diff_df$group == "Down_Regulated_Match")] <- 1
diff_df$AlphaValues[which(diff_df$group == "Up_Regulated_Match")] <- 1
diff_df$AlphaValues[which(diff_df$group == "Not_Significant_Match")] <- 1
diff_df$AlphaValues[which(diff_df$group == "Not_Significant_Gene")] <- 1


table(diff_df$AlphaValues)
# 0.01     1 
# 12801   126
# Find and label the top peaks..
numOfGenesToLabel = 60
group_numGenes = round(numOfGenesToLabel/3)  # there are 3 groups below that we label:
print(paste0("Please note the actual total number of labeled genes for this ", groupComparisonString, " Volcano Plot is: ", (3*group_numGenes)))


#these are the 3 groups below that we will label:
diff_df = diff_df[with(diff_df, order(-diff_df$empiricalLFC, diff_df$padj)),]

genesToSelectFromEachCategory = 20

upMatchDF = diff_df[which(diff_df$group == "Up_Regulated_Match"),]
upMatchDF = upMatchDF[order(upMatchDF$padj),]


upMatchDFpart2 = diff_df[which(diff_df$group == "Up_Regulated_Match"),]

upMatchDFpart2 = upMatchDF[order(upMatchDF$empiricalLFC, decreasing = TRUE),]

upMatchDF = upMatchDF[1:genesToSelectFromEachCategory,]


upMatchDFpart2 = upMatchDFpart2[1:genesToSelectFromEachCategory,]


downMatchDF = diff_df[which(diff_df$group == "Down_Regulated_Match"),]
downMatchDF = downMatchDF[order(downMatchDF$padj),]#1:genesToSelectFromEachCategory]    
downMatchDF = downMatchDF[1:genesToSelectFromEachCategory,]


nonSigMatchDF = diff_df[which(diff_df$group == "Not_Significant_Match"),]
nonSigMatchDF_part1 = nonSigMatchDF[order(nonSigMatchDF$empiricalLFC),] #1:(genesToSelectFromEachCategory/2)]    
nonSigMatchDF_part2 = nonSigMatchDF[order(nonSigMatchDF$empiricalLFC, decreasing = TRUE), ] #decreasing = TRUE),1:(genesToSelectFromEachCategory/2)]    
nonSigMatchDF_part1 = nonSigMatchDF_part1[1:genesToSelectFromEachCategory,]
nonSigMatchDF_part2 = nonSigMatchDF_part2[1:genesToSelectFromEachCategory,]

top_peaks <- upMatchDF #diff_df[which(diff_df$Significance == "Down Regulated Match"),][1:newGenes,] #diff_df[with(diff_df, order(diff_df$empiricalLFC, diff_df$padj)),][1:group_numGenes,]  # positive logs
top_peaks <- rbind(top_peaks, upMatchDFpart2) #diff_df[which(diff_df$Significance != "Up Regulated Match"),][1:newGenes,]) #diff_df[with(diff_df, order(diff_df$empiricalLFC, diff_df$padj)),][1:group_numGenes,]  # positive logs
top_peaks <- rbind(top_peaks, downMatchDF) #diff_df[which(diff_df$Significance != "Up Regulated Match"),][1:newGenes,]) #diff_df[with(diff_df, order(diff_df$empiricalLFC, diff_df$padj)),][1:group_numGenes,]  # positive logs
top_peaks <- rbind(top_peaks, nonSigMatchDF_part1)
top_peaks <- rbind(top_peaks, nonSigMatchDF_part2)
top_peaks = unique(top_peaks)

top_peaks = diff_df[grep("Up_Regulated_Match", diff_df$group),]
top_peaks = rbind(top_peaks, diff_df[grep("Down_Regulated_Match", diff_df$group),])
top_peaks = unique(top_peaks)

unique(diff_df$group)
newGenes =group_numGenes/2

nonMatches = diff_df[-grep("Match", diff_df$group),]
matches = diff_df[grep("Match", diff_df$group),]
dim(nonMatches) # 12801   112
dim(matches) # 126   112

diff_df = nonMatches
# so matches are plotted after and come on top more
diff_df = rbind(diff_df, matches)
head(top_peaks)

volcanoTitle = paste0("Volcano Plot: ", groupComparisonString) 
new_demo = ggplot(diff_df, aes(x = empiricalLFC, y = -log10(padj))) + geom_point(aes(color = Significance, alpha = AlphaValues, size = 2)) + scale_color_manual(values = c(downRegColor, noSigChangeColor, "khaki", "khaki", upRegColor)) + theme_bw(base_size = 12) + theme(legend.position = "bottom") + geom_text_repel(data = top_peaks, aes(label = external_gene_name), size = 5, box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"))
new_demo = new_demo + ylim(0, 10) + xlim(-3, 3.5)

mynamestheme <- theme(plot.title = element_text(family = "Helvetica", face = "bold", size = (25), hjust = 0.5), legend.title = element_text(colour = "steelblue", face = "bold.italic", family = "Helvetica", size = (15)), 
                      legend.text = element_text(face = "italic", colour = "steelblue4", family = "Helvetica", size = (17)),
                      axis.title = element_text(family = "Helvetica", size = (21), colour = "steelblue4"),
                      axis.text = element_text(family = "Courier", colour = "cornflowerblue", size = (19)))
final_plot = new_demo + mynamestheme + labs(title = volcanoTitle, x = expression(Log[2]~(Fold~Change)), y = expression(-Log[10]~(Adjusted~P-Value))) + guides(colour = guide_legend(override.aes = list(size = 10)))

final_plot = final_plot + geom_hline(yintercept=1, linetype="dotted", size = 1, colour = "steelblue4")


final_plot


png(outputPath, width = volcanoPlotWidth, height = volcanoPlotHeight)
print(final_plot)
dev.off()



#######################################################################

## Please note that we organize the up-regulated and down-regulated genes here
# for each brain region:
# Lateral Temporal Lobe: 
ltlUpDF = cbind(upRegMatchesLTL, rep("up-regulated", length(upRegMatchesLTL)),
                                     rep("LTL", length(upRegMatchesLTL)))
colnames(ltlUpDF) = c("geneName", "regulationDirection", "brainRegion")
ltlDownDF = cbind(downRegMatchesLTL, rep("down-regulated", length(downRegMatchesLTL)),
                rep("LTL", length(downRegMatchesLTL)))                

colnames(ltlDownDF) = c("geneName", "regulationDirection", "brainRegion")

ltlDEDF = rbind(ltlUpDF, ltlDownDF)



# Hippocampus: 
hippUpDF = cbind(upRegMatchesHipp, rep("up-regulated", length(upRegMatchesHipp)),
                rep("Hippocampus", length(upRegMatchesHipp)))

hippDownDF = cbind(downRegMatchesHipp, rep("down-regulated", length(downRegMatchesHipp)),
                  rep("Hippocampus", length(downRegMatchesHipp)))                

colnames(hippUpDF) = c("geneName", "regulationDirection", "brainRegion")
colnames(hippDownDF) = c("geneName", "regulationDirection", "brainRegion")

hippDEDF = rbind(hippUpDF, hippDownDF)

# DLPFC: 
dlpfcUpDF = cbind(upRegMatchesDLPFC, rep("up-regulated", length(upRegMatchesDLPFC)),
                 rep("DLPFC", length(upRegMatchesDLPFC)))

dlpfcDownDF = cbind(downRegMatchesDLPFC, rep("down-regulated", length(downRegMatchesDLPFC)),
                   rep("DLPFC", length(downRegMatchesDLPFC)))                
colnames(dlpfcUpDF) = c("geneName", "regulationDirection", "brainRegion")
colnames(dlpfcDownDF) = c("geneName", "regulationDirection", "brainRegion")


dlpfcDEDF = rbind(dlpfcUpDF, dlpfcDownDF)



allDE3BrainRegionsResultsDF = rbind(ltlDEDF, hippDEDF)
allDE3BrainRegionsResultsDF = rbind(allDE3BrainRegionsResultsDF, dlpfcDEDF)
allDE3BrainRegionsResultsDF
allDE3BrainRegionsResultsDF = unique(allDE3BrainRegionsResultsDF)
deOutputPathName_All3RegionsADCovidGRN = paste0(dataSetOutputPaths, "all3BrainRegionsDifferentiallyExpressedADCovidGRNGenes.csv")
print(paste0(":) Please note that the output is here: ", ddeOutputPathName_All3RegionsADCovidGRN)
dataSetOutputPaths
#write.csv(allDE3BrainRegionsResultsDF, "F://organizedAlzheimers//allDE3BrainRegionsResultsDF.csv")
write.csv(allDE3BrainRegionsResultsDF, deOutputPathName_All3RegionsADCovidGRN)
