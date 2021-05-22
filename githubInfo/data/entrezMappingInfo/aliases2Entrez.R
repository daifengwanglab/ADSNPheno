#https://cran.r-project.org/web/packages/aliases2entrez/aliases2entrez.pdf

#BiocManager::install("aliases2entrez")
library("aliases2entrez")
# import the correspondence file
file <- system.file("extdata", "HGNC.txt", package = "aliases2entrez")
HGNC <- read.delim(file)
# alterntatively update a new one with update_symbols()
symbols <- c("BRCA1", "TP53")
# run the main function
ids <- convert_symbols(symbols, HGNC)

symbolsDf = read.csv("D:\\organizedAlzheimers\\entrezMappingInfo\\entrezMappingInfoNewest.csv")
symbols = symbolsDf[,1]
ids <- convert_symbols(symbols, HGNC)
ids
write.csv(ids, "D:\\aliases2Entrezold.csv")
