# Please note that this file has all the packages that are needed
#BiocManager::install("GenomicInteractions")
#BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")

library("aliases2entrez")
library("biomaRt")
library("doRNG")
library("GENIE3")
library("GenomicInteractions")
library("hash")
library("lubridate")
library("motifbreakR")
library("rentrez")
library("reshape2")
library("reticulate")
library("stringr")
library("trena")
library("WGCNA")
library(clusterProfiler)
library(DOSE)
library(dplyr)
library(lubridate)
library(meshes)
library(MotifDb)
library(msigdbr)
library(RedeR)
library(RTN)
library(rvest)
library(stringr)

library('TxDb.Hsapiens.UCSC.hg19.knownGene')
library(MeSH.Hsa.eg.db)
library("org.Hs.eg.db")
library(BSgenome.Hsapiens.UCSC.hg19)
library("SNPlocs.Hsapiens.dbSNP144.GRCh37")


###
# Please note that here we are creating a Python virtual environment:
############################
# #install.packages("reticulate")
# # https://rstudio.github.io/reticulate/
# # create a new environment
virtualenv_create("r-reticulate")
#
# # install SciPy
virtualenv_install("r-reticulate", "scipy", python_version  = "3.8.5")
virtualenv_install("r-reticulate", "pandas", python_version  = "3.8.5")
virtualenv_install("r-reticulate", "numpy", python_version  = "3.8.5")
virtualenv_install("r-reticulate", "datetime", python_version  = "3.8.5")
virtualenv_install("r-reticulate", "fsspec", python_version  = "3.8.5")

use_virtualenv("r-reticulate")

py_install("pandas", python_version  = "3.8.5")
py_install("numpy", python_version  = "3.8.5")
py_install("math", python_version  = "3.8.5")
py_install("datetime", python_version  = "3.8.5")
py_install("sklearn", python_version  = "3.8.5")
py_install("fsspec", python_version  = "3.8.5")


#################################
# use_condaenv("myenv")
# # create a new environment
# conda_create("r-reticulate")
#
# # install SciPy
# conda_install("r-reticulate", "scipy", python_version  = "3.8.5")
# conda_install("r-reticulate", "pandas", python_version  = '3.8.5')
# conda_install("r-reticulate", "numpy", python_version  = "3.8.5")
# conda_install("r-reticulate", "datetime", python_version  = "3.8.5")
# #################################################

#
# # please indicate that we want to use a specific virtualenv
#use_virtualenv("r-reticulate")

py_install("pandas")
py_install("numpy")
py_install("math")
py_install("datetime")
py_install("sklearn")





#print(paste0(":) please note that there are ", length(keyModules), " gene modules that are associated with at least 1 of the phenotypes"))



#save(filePathInfoDF, file = filePathInfoDF_FileNameRData)
#write.csv(filePathInfoDF, filePathInfoDF_FileNameCSV)
