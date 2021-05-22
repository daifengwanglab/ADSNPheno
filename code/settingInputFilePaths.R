# please note that the user enters in the file paths for input files

colnamesVec = c("TF", "RegulatedGene", "Source", "Info", "CombinedName")
filePathOfGeneExpressionDataUsedByWGCNA = "D:\\organizedAlzheimers\\GeneExpressionPreparation\\outputs\\finalGeneExpressionData_AlzheimersDisease_HippocampusCa1_Brain_log2TransformedAndThenScaledInput_2020_11_20.csv"
geneInfoFilePath = "D:\\organizedAlzheimers\\entrezMappingInfo\\entrezMappingInfo.csv"

filePathOfRTN = "D:\\organizedAlzheimers\\TranscriptionalRegulator\\outputs\\finalRTN_"
filePathOfGenie3 = "D:\\organizedAlzheimers\\TranscriptionalRegulator\\outputs\\finalGenie3_"

infoDF = read.csv(geneInfoFilePath, header = TRUE)#[,c(1, 2)]#, row.names = 1)
colnames(infoDF) = c("geneSymbolName", "entrezID",          "status")

tfsDF <- read.csv("D:\\organizedAlzheimers\\TranscriptionalRegulator\\inputs\\LambertAndJasparTFs.csv", 
                  header = TRUE)[,c(1,2)]
tfs = tfsDF[,1]
#tfsDF <- read.csv("D:\\organizedAlzheimers\\TranscriptionalRegulator\\inputs\\LambertAndJasparTFs.csv", 
                  