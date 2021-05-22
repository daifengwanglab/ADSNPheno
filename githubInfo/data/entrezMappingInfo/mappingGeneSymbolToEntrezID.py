#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np

print(":) please note we are reading in function: geneNameToIDMappingFunctionForBigDataSets(infoDF, dfToMap, targetGeneColNameInDFToMap = 'target gene', tfColNameInDFToMap = 'TF', outputPath = '')"
# By Saniya Khullar :)
def geneNameToIDMappingFunctionForBigDataSets(infoDF, dfToMap, targetGeneColNameInDFToMap, tfColNameInDFToMap, outputPath = ""):
    mappingDF = infoDF.drop_duplicates()
    #mappingDF = pd.read_csv(mappingPath)


    # In[4]:


    #'Unnamed: 0'
    mappingDF.columns.tolist()


    # In[ ]:

    comboDF = dfToMap #pd.read_csv(dfInputFilePath)
    comboDF = comboDF.drop_duplicates()
    print("dfToMap:")
    print(comboDF.head())

    # In[ ]:


    newDF = comboDF.merge(mappingDF, left_on=str(targetGeneColNameInDFToMap), right_on='geneSymbolName')
    newDF = newDF.rename(columns={"entrezID": "regulatedGene_entrezID"})
    newDF = newDF.merge(mappingDF, left_on=str(tfColNameInDFToMap), right_on='geneSymbolName')
    newDF = newDF.rename(columns={"entrezID": "currentTF_entrezID"})


    # In[ ]:


    if 'Unnamed: 0' in newDF.columns.tolist():
        newDF = newDF.drop(columns=['Unnamed: 0'])
    newDF = newDF.drop(columns=['status_x', 'status_y', 'geneSymbolName_x', 'geneSymbolName_y', 'Unnamed: 0'])
    print(":) Please note that the newDF is:", newDF.shape[0], " rows and ", newDF.shape[1], "columns! :)")
    print(":) Please note the 1st 5 rows:")
    print(newDF.head())
    if outputPath != "":
        newDF.to_csv(outputPath)
        print(":) Please note we finished writing out the output :)")
    return(newDF)


print(":) please note we have finished reading in function: geneNameToIDMappingFunctionForBigDataSets(infoDF, dfToMap, targetGeneColNameInDFToMap = 'target gene', tfColNameInDFToMap = 'TF', outputPath = '')"
