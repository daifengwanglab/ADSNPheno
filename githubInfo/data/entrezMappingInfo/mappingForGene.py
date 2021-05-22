#!/usr/bin/env python
# coding: utf-8

# In[64]:


import pandas as pd
import numpy as np
#import math


print("Please note that this mapping function helps us get EntrezIDs for our data")
print(":) Please note we are reading in this function: mappingGenesToIDsTranscriptionalDF(transcriptionalDF, infoDF, outputNameTranscriptional)")
print("reading in function...")

def mappingGenesToIDsTranscriptionalDF(transcriptionalDF, infoDF, outputNameTranscriptional):
    infoDF = infoDF.drop_duplicates()
    geneSymbolList = infoDF["geneSymbolName"].tolist()
    entrezIDList = infoDF["entrezID"].tolist()
    geneSymbolToIDDict = {}
    for i in range(0, len(geneSymbolList)):
        if i % 10000 == 0:
            print(":) i = ", i)
        geneSymbol = geneSymbolList[i]
        geneId = entrezIDList[i]
        geneSymbolToIDDict[geneSymbol] = geneId
    print("geneSymbolToIDDict['C1orf35']", geneSymbolToIDDict['C1orf35'])

    geneSymbolToIDDict['MARCH1'] =	55016
    geneSymbolToIDDict['MARCH10'] =	162333
    geneSymbolToIDDict['MARCH11'] =	441061
    geneSymbolToIDDict['MARCH2'] =	51257
    geneSymbolToIDDict['MARCH3'] =	115123
    geneSymbolToIDDict['MARCH4'] =	57574
    geneSymbolToIDDict['MARCH5'] =	54708
    geneSymbolToIDDict['MARCH6'] =	10299
    geneSymbolToIDDict['MARCH7'] =	64844
    geneSymbolToIDDict['MARCH8'] =	220972
    geneSymbolToIDDict['MARCH9'] =	92979
    geneSymbolToIDDict['MARC1'] =	64757
    geneSymbolToIDDict['MARC2'] =	54996
    geneSymbolToIDDict['KMT2C'] =	58508
    geneSymbolToIDDict['MLL3'] =	58508
    transcriptionalDFList_TF = transcriptionalDF['currentTF'].tolist()
    transcriptionalDFList_RegulatedGene = transcriptionalDF['regulatedGene'].tolist()

    print(len(transcriptionalDFList_TF))
    print(len(transcriptionalDFList_RegulatedGene))

    transcriptionalDFListIDs_TF = []
    transcriptionalDFListIDs_RegulatedGene = []
    transcriptionalDFListIDs_Combined = []


    for i in range(0, len(transcriptionalDFList_TF)):
        if i % 1000 == 0:
            print(":) i = ", i)
        geneSymbolTF = transcriptionalDFList_TF[i]
        

        if geneSymbolTF == "NAN":
            geneIDTF = 11280
        else:
            geneIDTF = geneSymbolToIDDict[geneSymbolTF]
        transcriptionalDFListIDs_TF.append(geneIDTF)

        geneSymbolRegGene = transcriptionalDFList_RegulatedGene[i]

        if geneSymbolRegGene == "NAN":
            geneIDRegGene = 11280
        else:
            geneIDRegGene = geneSymbolToIDDict[geneSymbolRegGene]
        transcriptionalDFListIDs_RegulatedGene.append(geneIDRegGene)
        comboName = str(geneIDTF) + " || " + str(geneIDRegGene)
        transcriptionalDFListIDs_Combined.append(comboName)

    transcriptionalDF['currentTF_entrezID'] = transcriptionalDFListIDs_TF
    transcriptionalDF['regulatedGene_entrezID'] = transcriptionalDFListIDs_RegulatedGene
    transcriptionalDF['CombinedEntrezName'] = transcriptionalDFListIDs_Combined
    print(":) Please note transcriptionalDF:", transcriptionalDF.shape[0])
    print(transcriptionalDF.head())
    transcriptionalDF.to_csv(outputNameTranscriptional)
    print(":) please note our output path is: ", outputNameTranscriptional)
    return transcriptionalDF

print(":) Please note we are reading in this function: mappingGenesToIDsScGRNDF(scgrnOutputDF, infoDF, outputNameScgrn)")
print("reading in function...")
def mappingGenesToIDsScGRNDF(scgrnOutputDF, infoDF, outputNameScgrn):
    infoDF = infoDF.drop_duplicates()

    geneSymbolList = infoDF["geneSymbolName"].tolist()
    entrezIDList = infoDF["entrezID"].tolist()
    geneSymbolToIDDict = {}
    for i in range(0, len(geneSymbolList)):
        if i % 10000 == 0:
            print(":) i = ", i)
        geneSymbol = geneSymbolList[i]
        geneId = entrezIDList[i]
        geneSymbolToIDDict[geneSymbol] = geneId
    print(geneSymbolToIDDict)
    geneSymbolToIDDict['MARCH1'] =	55016
    geneSymbolToIDDict['MARCH10'] =	162333
    geneSymbolToIDDict['MARCH11'] =	441061
    geneSymbolToIDDict['MARCH2'] =	51257
    geneSymbolToIDDict['MARCH3'] =	115123
    geneSymbolToIDDict['MARCH4'] =	57574
    geneSymbolToIDDict['MARCH5'] =	54708
    geneSymbolToIDDict['MARCH6'] =	10299
    geneSymbolToIDDict['MARCH7'] =	64844
    geneSymbolToIDDict['MARCH8'] =	220972
    geneSymbolToIDDict['MARCH9'] =	92979
    geneSymbolToIDDict['MARC1'] =	64757
    geneSymbolToIDDict['MARC2'] =	54996
    geneSymbolToIDDict['KMT2C'] =	58508
    geneSymbolToIDDict['MLL3'] =	58508
    scgrnOutputDFList_TF = scgrnOutputDF['currentTF'].tolist()
    scgrnOutputDFList_RegulatedGene = scgrnOutputDF['regulatedGene'].tolist()
    print(len(scgrnOutputDFList_TF))
    print(len(scgrnOutputDFList_RegulatedGene))

    scgrnOutputDFListIDs_TF = []
    scgrnOutputDFListIDs_RegulatedGene = []
    scgrnOutputDFListIDs_Combined = []


    for i in range(0, len(scgrnOutputDFList_TF)):
        if i % 1000 == 0:
            print(":) i = ", i)
        #print(i)
        geneSymbolTF = scgrnOutputDFList_TF[i]

        if geneSymbolTF == "NAN":
            geneIDTF = 11280
        else:
            geneIDTF = geneSymbolToIDDict[geneSymbolTF]
        scgrnOutputDFListIDs_TF.append(geneIDTF)

        if geneSymbolRegGene == "NAN":
            geneIDRegGene = 11280
        else:
            geneIDRegGene = geneSymbolToIDDict[geneSymbolRegGene]

        scgrnOutputDFListIDs_RegulatedGene.append(geneIDRegGene)  
        comboName = str(geneIDTF) + " || " + str(geneIDRegGene)

        scgrnOutputDFListIDs_Combined.append(comboName)
    scgrnOutputDF['currentTF_entrezID'] = scgrnOutputDFListIDs_TF
    scgrnOutputDF['regulatedGene_entrezID'] = scgrnOutputDFListIDs_RegulatedGene
    scgrnOutputDF['CombinedEntrezName'] = scgrnOutputDFListIDs_Combined
    print(":) Please note scgrnOutputDF:", scgrnOutputDF.shape[0])
    print(scgrnOutputDF.head())
    scgrnOutputDF.to_csv(outputNameScgrn)
    print(":) please note our output path is: ", outputNameScgrn)
    return scgrnOutputDF


print(":) Please note we are reading in both functions: ")
print(":) mappingGenesToIDsTranscriptionalDF(transcriptionalDF, infoDF, outputNameTranscriptional)")
print("mappingGenesToIDsScGRNDF(scgrnOutputDF, infoDF, outputNameScgrn)")
print("reading in function...")
