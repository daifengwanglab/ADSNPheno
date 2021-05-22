library(biomaRt);

#sub("[*.]")
#unlist(strsplit("ENSG00000167578.11", "[.]"))[1]
geneExpressionDF = read.csv("D://organizedAlzheimers//dlpfcRosmapJan2021Data.csv", header = TRUE)




# protein coding only
ensembl_mart=useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl",host="grch37.ensembl.org",mirror = "useast")
ensembl_types=listAttributes(ensembl_mart);
#ensembl_types[grep('entrez',ensembl_types$name),]
hid=getBM(attributes=c('ensembl_gene_id','hgnc_symbol','entrezgene_id','gene_biotype','ensembl_transcript_id'), mart = ensembl_mart)
hid
#write.csv(hid, "D://hid_ensemblIDs.csv")

gettingEnsembleIDsForTheGeneExpressionDF <- function(geneExpressionDF, ensemblIDCol = 2,
                                                     outputName = "D://organizedAlzheimers//dlpfcGeneExpressionDataset.csv"){
print(":) please note this function by Saniya helps get the EntrezIDs and Gene Names for the gene expression dataset.")
  library("biomaRt")
ensemblGeneIDDF = geneExpressionDF[,ensemblIDCol]
allIndexesOfEnsembls = seq(1, length(ensemblGeneIDDF), by = 1)

indexesOfEnsemblsWithDecimals = grep("[.]", ensemblGeneIDDF)
print(paste0(":) Please note that ", length(indexesOfEnsemblsWithDecimals), " of ",
             length(allIndexesOfEnsembls), " have decimal places in Ensembl IDs, which we will remove."))
indexesOfEnsemblsWithNoDecimals = setdiff(allIndexesOfEnsembls, indexesOfEnsemblsWithDecimals)

ensemblIDsWithDecimals = ensemblGeneIDDF[indexesOfEnsemblsWithDecimals]
brokenIDs = unlist(strsplit(ensemblGeneIDDF, "[.]"))#[1]
idsToUse = seq(1, length(brokenIDs), by = 2)
ensemblIDs = brokenIDs[idsToUse]

if (length(indexesOfEnsemblsWithNoDecimals) > 0){
  print(paste0(":) Please note that ", length(indexesOfEnsemblsWithNoDecimals), " of ",
               length(allIndexesOfEnsembls), " do NOT have decimal places in Ensembl IDs."))
  ensemblsWithNoDecimals = ensemblGeneIDDF[indexesOfEnsemblsWithNoDecimals]
    
  ensemblIDs = c(ensemblIDs, ensemblsWithNoDecimals)
}

geneExpressionDF$simplifiedEnsemblID = ensemblIDs
ensembl_mart=useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl",host="grch37.ensembl.org",mirror = "useast")


genesALL <- getBM(filters = "ensembl_gene_id",
                  attributes = c("ensembl_gene_id","entrezgene_id", "hgnc_symbol"),
                  values = ensemblIDs, 
                  mart = ensembl_mart)
colnames(genesALL) = c("simplifiedEnsemblID", "entrezID", "geneName")
combinedDF = merge(genesALL, geneExpressionDF, by = c("simplifiedEnsemblID"))
#View(combinedDF)
head(combinedDF)
write.csv(combinedDF, outputName)
print(paste0(":) Please note that we wrote out the output here: ", outputName))
return(combinedDF)

}
#install.packages("gbutils")
library("gbutils")

#missingEntrezID = genesALL[intersect(which(isNA(combinedDF$entrezID)),]

missingDF = read.csv("D://organizedAlzheimers//dlpfcMissingInfo.csv", header = TRUE)
missingDF

library("org.Hs.eg.db")
newEntrezIDs = mapIds(org.Hs.eg.db, missingDF[,2], 'ENTREZID', 'SYMBOL')
#BiocManager::install("org.Hs.egENSEMBL")
library("org.Hs.egENSEMBL")
newDF = data.frame(newEntrezIDs)
newDF$geneName = missingDF[,2]
write.csv(newDF, "D://organizedAlzheimers//dlpfcMissing.csv")

datanewEntrezIDs

newEntrezIDs1 = mapIds(org.Hs.eg.db, ensemblIDs, 'ENSEMBL', 'SYMBOL')

library("EnsDb.Hsapiens.v75")
newEntrezIDs1 = mapIds(EnsDb.Hsapiens.v75, ensemblIDs, 'ENSEMBL', 'SYMBOL')



newEntrezIDs1 = mapIds(org.Hs.eg.db, ensemblIDs, 'ENSEMBL', 'SYMBOL')

newEntrezIDs = mapIds(ensemblGeneIDDF, newSymbols, 'ENTREZID', 'SYMBOL')


newEntrezIDs = mapIds(org.Hs.eg.db, ensemblIDs, 'ENSEMBL', 'ENTREZID')

#BiocManager::install("EnsDb.Hsapiens.v75")
library(EnsDb.Hsapiens.v75)



#BiocManager::install("EnrichmentBrowser")
library("EnrichmentBrowser")
map.ids( eset, org="hsa", from="ENSEMBL", to="ENTREZID" )

eset <- make.example.data("eset", nfeat=3, nsmpl=3)
featureNames(eset) <- paste0("ENSG00000000", c("003","005", "419"))
eset <- map.ids(eset, org="hsa")


genesALL <- getBM(filters = "ensembl_gene_id",
                  attributes = c("hgnc_symbol","entrezgene_id"),
                  values = missingDF[,2], 
                  mart = ensembl_mart)

updatedGeneInfoDF = gettingEnsembleIDsForTheGeneExpressionDF(geneExpressionDF, ensemblIDCol = 2,
                                         outputName = "D://organizedAlzheimers//dlpfcGeneExpressionDataset.csv")



#grep(ensemblGeneIDDF, "[.]", )

unlist(strsplit("ENSG00000167578.11", "[.]"))[1]



ensembl_mart=useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl",host="grch37.ensembl.org",mirror = "useast")

genes <- getBM(filters = "ensembl_gene_id",
               attributes = c("ensembl_gene_id","entrezgene_id"),
               values = ensemblIDs, 
               mart = ensembl_mart)

genes

write.csv(genes, "D://organizedAlzheimers//dlpfcEntrezIDs.csv")

ensembl_mart=useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl",host="grch37.ensembl.org",mirror = "useast")


genes <- getBM(filters = "ensembl_gene_id",
               attributes = c("ensembl_gene_id","entrezgene_id"),
               values = ensemblIDs, 
               mart = ensembl_mart)

genes


geneNames <- getBM(filters = "ensembl_gene_id",
                   attributes = c("ensembl_gene_id","hgnc_symbol"),
                   values = ensemblIDs, 
                   mart = ensembl_mart)

geneNames


genesALL <- getBM(filters = "ensembl_gene_id",
                  attributes = c("ensembl_gene_id","entrezgene_id", "hgnc_symbol"),
                  values = ensemblIDs, 
                  mart = ensembl_mart)





"hgnc_symbol"

write.csv(geneNames, "D://organizedAlzheimers//dlpfcEntrezIDsNEW.csv")

# consider log2(x+1)
