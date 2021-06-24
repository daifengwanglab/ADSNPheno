print(":) Please note we are reading in this function: organizingOutputFileForTFToGeneRelationship(grnResultsCombinedFileNameCSV, parentName, disease, bodyRegion, dataScalingOutputMini, geneToEntrezIDMappingPath, minNumSourcesGeneGRN, progressBarPythonIterations)")
#def organizingOutputFileForTFToGeneRelationship(fileName, parentName, disease, bodyRegion, dataScaling, powerVal, newDate):


##import numpy as np
import pandas as pd
import datetime
from datetime import datetime
#from sklearn.preprocessing import LabelEncoder
#import matplotlib.pyplot as plt
#from sklearn.linear_model import LogisticRegression
#from sklearn.metrics import classification_report, confusion_matrix
##from numpy import mean
##from numpy import std
# Classification and ROC analysis
from itertools import cycle
#from sklearn.feature_selection import RFE
#from sklearn import svm, datasets
#from sklearn.metrics import auc
#from sklearn.metrics import plot_roc_curve
#from sklearn.model_selection import StratifiedKFold
#from sklearn import svm, datasets
#from sklearn.metrics import roc_curve, auc
#from sklearn.model_selection import train_test_split
#from sklearn.preprocessing import label_binarize
#from sklearn.multiclass import OneVsRestClassifier
#from scipy import interp
#from sklearn.metrics import roc_auc_score
#from sklearn.svm import SVC
#from sklearn.model_selection import StratifiedKFold
#from sklearn.feature_selection import RFECV
#from sklearn.datasets import make_classification
#from sklearn.model_selection import LeaveOneOut
#from sklearn.pipeline import Pipeline
#from sklearn.model_selection import cross_val_score, GridSearchCV, train_test_split


#seedVal = 123
#np.random.seed(seedVal)


def sourcesToRemoveListFunction(minNumSourcesGeneGRN):
    sourcesToRemoveList = []
    if minNumSourcesGeneGRN == 2:
        sourcesToRemoveList = ['Trust2', 'Trrust2 || Genie3', 'Trrust2 || TReNa', 'Trrust2', 'Genie3 || Trrust2', 'TReNa || Trrust2', 'Trrust2 || RTN', 'RTN || Trrust2']
    elif minNumSourcesGeneGRN == 1:
        sourcesToRemoveList = ['Trrust2']
    elif minNumSourcesGeneGRN == 3:
        sourcesToRemoveList = ['Trrust2', 'Trrust2 || Genie3', 'Trrust2 || TReNa', 'Trrust2', 'Genie3 || Trrust2', 'TReNa || Trrust2', 'Trrust2 || RTN', 'RTN || Trrust2',
        'Trrust2 || Genie3 || RTN', 'Trrust2 || RTN || Genie3', 'RTN || Trrust2 || Genie3',  'Genie3 || Trrust2 || RTN',  'RTN || Genie3 || Trrust2',  'Genie3 || RTN || Trrust2', 
        'Trrust2 || Genie3 || TReNA', 'Trrust2 || TReNA || Genie3', 'TReNA || Trrust2 || Genie3',  'Genie3 || Trrust2 || TReNA',  'TReNA || Genie3 || Trrust2',  'Genie3 || TReNA || Trrust2', 
        'Trrust2 || TReNA || RTN', 'Trrust2 || RTN || TReNA', 'RTN || Trrust2 || TReNA',  'TReNA || Trrust2 || RTN',  'RTN || TReNA || Trrust2',  'TReNA || RTN || Trrust2']
        
    print(":) minNumSourcesGeneGRN = ", minNumSourcesGeneGRN, " so sourcesToRemoveList = ", sourcesToRemoveList)
    return sourcesToRemoveList

def repeatVal(torep = "promoter", nrep = 4):
    repetitions=[]
    for i in range(nrep):
        i=torep
        repetitions.append(i)
    return repetitions



def pleaseGetCurrentDateInfo():
    from datetime import datetime
    # https://stackoverflow.com/questions/13650214/how-to-convert-a-python-datetime-datetime-into-a-list
    dateInfo = datetime.now()
    month = dateInfo.month
    if month == 1:
        monthName = "Jan_"
    elif month == 2:
        monthName = "Feb_"   
    elif month == 3:
        monthName = "Mar"
    elif month == 4:
        monthName = "Apr"   
    elif month == 5:
        monthName = "May"        
    elif month == 6:
        monthName = "Jun"   
    elif month == 7:
        monthName = "Jul"
    elif month == 8:
        monthName = "Aug"   
    elif month == 9:
        monthName = "Sep"   
    elif month == 10:
        monthName = "Oct"   
    elif month == 11:
        monthName = "Nov"       
    elif month == 12:
        monthName = "Dec"       
               
    dayInfo = monthName + str(dateInfo.day) + "_" + str(dateInfo.year)
    return(dayInfo)


def addingAdditionalGeneInfo(geneToEntrezIDMappingPath):
    
    #geneToEntrezIDMappingPath = "F://organizedAlzheimers//entrezMappingInfo//infoDFMappingsGeneSymbolAndID.csv"
    infoDF = pd.read_csv(geneToEntrezIDMappingPath)
    infoDF
    print(":) Please note that before our infoDF (geneName, entrezID mapping dataframe) had: ", infoDF.shape[0], " rows")

    # please note that we are removing withdrawn and/or unmatched and/or unknown genes from our list here
    infoDF = infoDF[infoDF['entrezID'] != "Withdrawn"]
    infoDF = infoDF[infoDF['entrezID'] != "unmatched"]
    infoDF = infoDF[infoDF['entrezID'] != "unknown"]

    infoDF
    # please note that this function by Saniya will add some additional genes info to a dataframe just to be safe :)
    geneSymbolToIDDict = {}
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

    geneSymbolToIDDict['SEPTIN10'] = 151011
    geneSymbolToIDDict['SEPTIN11'] = 55752
    geneSymbolToIDDict['SEPTIN8'] = 23176
    geneSymbolToIDDict['SEPTIN7'] = 989
    geneSymbolToIDDict['SEPTIN6'] = 23157
    geneSymbolToIDDict['SEPTIN1'] = 1731
    geneSymbolToIDDict['SEPTIN4'] = 5414
    geneSymbolToIDDict['SEPTIN3'] = 55964
    geneSymbolToIDDict['SEPTIN12'] = 12404
    geneSymbolToIDDict['DEC1'] = 50514

    geneSymbolToIDDict['NR1A4'] = 3164


    geneSymbolToIDDict['1-Mar'] =	64757
    geneSymbolToIDDict['10-Mar'] =	162333
    geneSymbolToIDDict['11-Mar'] =	441061
    geneSymbolToIDDict['2-Mar'] =	54996
    geneSymbolToIDDict['3-Mar'] =	115123
    geneSymbolToIDDict['4-Mar'] =	57574
    geneSymbolToIDDict['5-Mar'] =	54708
    geneSymbolToIDDict['6-Mar'] =	10299
    geneSymbolToIDDict['7-Mar'] =	64844
    geneSymbolToIDDict['8-Mar'] =	220972
    geneSymbolToIDDict['9-Mar'] =	92979


    geneSymbolToIDDict['12-Sep'] = 12404
    geneSymbolToIDDict['10-Sep'] = 151011
    geneSymbolToIDDict['11-Sep'] = 55752
    geneSymbolToIDDict['8-Sep'] = 23176
    geneSymbolToIDDict['7-Sep'] = 989
    geneSymbolToIDDict['6-Sep'] = 23157
    geneSymbolToIDDict['1-Sep'] = 1731
    geneSymbolToIDDict['4-Sep'] = 5414
    geneSymbolToIDDict['3-Sep'] = 55964
    geneSymbolToIDDict['1-Dec'] = 50514

    counter = 0
    newDataFrame = pd.DataFrame()
    for geneSymbol in geneSymbolToIDDict.keys():
        idVal = str(geneSymbolToIDDict[geneSymbol])
        row = [geneSymbol, idVal] #, geneSymbol, idVal, geneSymbol, idVal]
        rowToAdd = pd.DataFrame(row).transpose()
        #rowToAdd.columns = ['Symbols', 'entrezID'] #, 'targetGeneName', 'targetGeneEntrezId', 'TFName', 'TFEntrezId']
        newDataFrame = newDataFrame.append(rowToAdd)
        #geneMappingInfoDF.append(rowToAdd)
        counter = counter + 1
        if counter % 1000 == 0:
            print(":) counter = ", counter)
    newDataFrame.columns = ['Symbols', 'entrezID']
    newDataFrame
    
    updatedDF = infoDF
    updatedDF = updatedDF.append(newDataFrame)
    print("Please note that after removing unknown/unmatched/withdrawn genes, our infoDF had: ", infoDF.shape[0], " rows")
    print("Please note that after adding in some special date-related genes and their respective gene entrez IDs, our infoDF had: ", updatedDF.shape[0], " rows")

    updatedDF["targetGeneName"] = updatedDF["Symbols"]
    updatedDF["targetGeneEntrezId"] = updatedDF["entrezID"]
    updatedDF["TFName"] = updatedDF["Symbols"]
    updatedDF["tFEntrezId"] = updatedDF["entrezID"]

    updatedDF = updatedDF.set_index(['targetGeneName'])

    return updatedDF


def retrieveScgrnomOutputWithMappedEntrezIDs(scgrnomOutputFileName, geneToEntrezIDMappingPath):
    # organize scgrn output
    
    scgrnomDF = pd.read_csv(scgrnomOutputFileName)

    if 'Unnamed: 0' in scgrnomDF.columns:
        scgrnomDF = scgrnomDF.drop(['Unnamed: 0'], axis = 1)


    scgrnomDF

    infoDF = addingAdditionalGeneInfo(geneToEntrezIDMappingPath)
    infoDF

    regGeneMappingInfoDF = infoDF
    regGeneMappingInfoDF = regGeneMappingInfoDF[['targetGeneEntrezId']]
    regGeneMappingInfoDF


    tfMappingInfoDF = infoDF
    tfMappingInfoDF['TFName'] = tfMappingInfoDF['Symbols']
    tfMappingInfoDF['TFEntrezId'] = tfMappingInfoDF['entrezID']
    tfMappingInfoDF = tfMappingInfoDF[['TFName',  'TFEntrezId']]
    tfMappingInfoDF = tfMappingInfoDF.set_index(['TFName'])

    tfMappingInfoDF
    
    if 'target gene' in scgrnomDF.columns.tolist():
        scgrnomDF["targetGeneName"] = scgrnomDF["target gene"].tolist()
    elif 'targetGene' in scgrnomDF.columns.tolist():
        scgrnomDF["targetGeneName"] = scgrnomDF["targetGene"].tolist()
        
    scgrnomDF["TFName"] = scgrnomDF["TF"].tolist()

    newmasterDF = scgrnomDF.merge(tfMappingInfoDF, left_on='TF', right_on='TFName') #on = "targetGeneName")

    if 'target gene' in newmasterDF.columns.tolist():
        newmasterDF = newmasterDF.merge(regGeneMappingInfoDF, left_on='target gene', right_on='targetGeneName') #on = "targetGeneName")
    else:
        newmasterDF = newmasterDF.merge(regGeneMappingInfoDF, left_on='targetGene', right_on='targetGeneName') #on = "targetGeneName")

    newmasterDF["CombinedEntrezName"] = newmasterDF['TFEntrezId'] + ' || ' + newmasterDF['targetGeneEntrezId'] #.agg(' || '.join, axis=1)


    newmasterDF[newmasterDF["TF"] == newmasterDF["TFName"]]

    newmasterDF = newmasterDF[newmasterDF["TF"] == newmasterDF["TFName"]]
    newmasterDF = newmasterDF[newmasterDF["target gene"] == newmasterDF["targetGeneName"]]
    newmasterDF = newmasterDF.drop(columns = ['targetGeneName', 'TFName'], axis = 1)
    newmasterDF = newmasterDF.drop_duplicates()
    newmasterDF # 1686046  to 1685344 
    
    return newmasterDF



#scgrnomOutputFileName = "F://scGRNPart2//updatedNetworkForChromosome18.csv"

def retrieveScgrnomOutputWithMappedEntrezIDs(scgrnomOutputFileName, geneToEntrezIDMappingPath):
    # organize scgrn output
    
    scgrnomDF = pd.read_csv(scgrnomOutputFileName)

    if 'Unnamed: 0' in scgrnomDF.columns:
        scgrnomDF = scgrnomDF.drop(['Unnamed: 0'], axis = 1)


    scgrnomDF

    infoDF = addingAdditionalGeneInfo(geneToEntrezIDMappingPath)
    infoDF

    regGeneMappingInfoDF = infoDF
    regGeneMappingInfoDF = regGeneMappingInfoDF[['targetGeneEntrezId']]
    regGeneMappingInfoDF


    tfMappingInfoDF = infoDF
    tfMappingInfoDF['TFName'] = tfMappingInfoDF['Symbols']
    tfMappingInfoDF['TFEntrezId'] = tfMappingInfoDF['entrezID']
    tfMappingInfoDF = tfMappingInfoDF[['TFName',  'TFEntrezId']]
    tfMappingInfoDF = tfMappingInfoDF.set_index(['TFName'])

    tfMappingInfoDF
    
    if "targetGene" in scgrnomDF.columns.tolist():
        scgrnomDF["targetGeneName"] = scgrnomDF["targetGene"].tolist()
    else:# "target gene" in scgrnomDF.columns.tolist():
        scgrnomDF["targetGeneName"] = scgrnomDF["target gene"].tolist()
    scgrnomDF["TFName"] = scgrnomDF["TF"].tolist()

    newmasterDF = scgrnomDF.merge(tfMappingInfoDF, left_on='TF', right_on='TFName') #on = "targetGeneName")
    if "target gene" in newmasterDF.columns.tolist():
        newmasterDF = newmasterDF.merge(regGeneMappingInfoDF, left_on='target gene', right_on='targetGeneName') #on = "targetGeneName")
    else:
        newmasterDF = newmasterDF.merge(regGeneMappingInfoDF, left_on='targetGene', right_on='targetGeneName') #on = "targetGeneName")

    newmasterDF["CombinedEntrezName"] = newmasterDF['TFEntrezId'] + ' || ' + newmasterDF['targetGeneEntrezId'] #.agg(' || '.join, axis=1)


    newmasterDF[newmasterDF["TF"] == newmasterDF["TFName"]]

    newmasterDF = newmasterDF[newmasterDF["TF"] == newmasterDF["TFName"]]
    
    
    if "targetGene" in newmasterDF.columns.tolist():
        newmasterDF = newmasterDF[newmasterDF["targetGene"] == newmasterDF["targetGeneName"]]
    elif "target gene" in newmasterDF.columns.tolist():
        newmasterDF = newmasterDF[newmasterDF["target gene"] == newmasterDF["targetGeneName"]]
 
    #newmasterDF = newmasterDF[newmasterDF["target gene"] == newmasterDF["targetGeneName"]]
    newmasterDF = newmasterDF.drop(columns = ['targetGeneName', 'TFName'], axis = 1)
    newmasterDF = newmasterDF.drop_duplicates()
    newmasterDF # 1686046  to 1685344 
    
    return newmasterDF
    

    
def pleaseGetOutputNameForFinalFullNetwork(chromName, parentName, disease, bodyRegion, dataScalingOutputMini, minNumSourcesGeneGRN):
    if minNumSourcesGeneGRN > 1:
        outputName = parentName + "fullGeneRegulatoryNetwork_GenieRTNTrenaTrrust2_Scgrnom_TFtoGenes_"  + str(chromName) + "_" + disease + "_" + bodyRegion + "_" + dataScalingOutputMini + "_" + str(minNumSourcesGeneGRN) + ".0_minSources.csv"

    else:
        outputName = parentName + "fullGeneRegulatoryNetwork_GenieRTNTrenaTrrust2_Scgrnom_TFtoGenes_"  + str(chromName) + "_" + disease + "_" + bodyRegion + "_" + dataScalingOutputMini + ".csv"
    
    return outputName

def fullNetworkDFForChromFunction(scgrnomOutputFileName, grnResultsCombinedFileNameCSV, parentName, disease, bodyRegion, dataScalingOutputMini, geneToEntrezIDMappingPath, minNumSourcesGeneGRN, enhancer_buffer_kbp, promoter_buffer_kbp, progressBarPythonIterations = 1000):#, sourcesToRemoveList):

#masterDF = organizingOutputFileForTFToGeneRelationship(grnResultsCombinedFileNameCSV, parentName, disease, bodyRegion, dataScalingOutputMini, geneToEntrezIDMappingPath, progressBarPythonIterations = 1000)


    
    try:
        outputName = pleaseGetOutputNameForFinalComboGRNGeneExpression(parentName, disease, bodyRegion, dataScalingOutputMini, minNumSourcesGeneGRN)

        masterDF = pd.read_csv(outputName)
        print(":) Yay, we could load it in here!")
    except:
        print(":( Please note we cannot find the file here: ", outputName, " so we are gathering info for this df")
        masterDF = organizingOutputFileForTFToGeneRelationship(grnResultsCombinedFileNameCSV, parentName, disease, bodyRegion, dataScalingOutputMini, geneToEntrezIDMappingPath, minNumSourcesGeneGRN, progressBarPythonIterations)
    if 'Unnamed: 0' in masterDF.columns:
        masterDF = masterDF.drop(['Unnamed: 0'], axis = 1)
        
    scgrnomDF = retrieveScgrnomOutputWithMappedEntrezIDs(scgrnomOutputFileName, geneToEntrezIDMappingPath)
    scgrnomDF

    fullNetworkDFForChrom = masterDF.merge(scgrnomDF, left_on='CombinedEntrezName', right_on='CombinedEntrezName')
    #fullNetworkDFForChrom = fullNetworkDFForChrom[fullNetworkDFForChrom["targetGeneEntrezId_x"] == fullNetworkDFForChrom["targetGeneEntrezId_y"]]
    #fullNetworkDFForChrom = fullNetworkDFForChrom[fullNetworkDFForChrom["TFEntrezId_x"] == fullNetworkDFForChrom["TFEntrezId_y"]]

    ##fullNetworkDFForChrom = fullNetworkDFForChrom[fullNetworkDFForChrom["target gene"] == fullNetworkDFForChrom["regulatedGene"]]
    ##fullNetworkDFForChrom = fullNetworkDFForChrom[fullNetworkDFForChrom["currentTF"] == fullNetworkDFForChrom["TF"]]
    fullNetworkDFForChrom = fullNetworkDFForChrom.rename(columns={"TF_x": "TF"})
    fullNetworkDFForChrom = fullNetworkDFForChrom.rename(columns={"targetGene": "target gene"})#, "targetGeneEntrezId_x": "targetGeneEntrezId"})

    fullNetworkDFForChrom = fullNetworkDFForChrom.drop(["targetGeneEntrezId_y", "TFEntrezId_y", "target gene", "TF"], axis = 1)
    fullNetworkDFForChrom = fullNetworkDFForChrom.rename(columns={"TFEntrezId_x": "TFEntrezId", "targetGeneEntrezId_x": "targetGeneEntrezId"})
    fullNetworkDFForChrom = fullNetworkDFForChrom.drop_duplicates()
#     fullNetworkDFForChrom # 96024 
#     if fullNetworkDFForChrom.shape[0] == 0:
#         return "None"

    try:
        chromName = str(fullNetworkDFForChrom["chromosome"].tolist()[0])
        fullNetworkDFForChrom[["new_start", "new_end"]] = fullNetworkDFForChrom.apply(lambda row: pd.Series(do_widen_intervalFullNet(row, enhancer_buffer_kbp, promoter_buffer_kbp)), 
                                           axis=1)
        fullNetworkDFForChrom["tfAndChrom"] = fullNetworkDFForChrom["chromosome"] + "_" + fullNetworkDFForChrom["TFEntrezId"].values.astype(str)
    
    except:
        chromName = ""
        fullNetworkDFForChrom[["new_start", "new_end"]] = ["",""]
        fullNetworkDFForChrom["tfAndChrom"] = fullNetworkDFForChrom["chromosome"] + "_" #+ fullNetworkDFForChrom["TFEntrezId"].values.astype(str)

    outputName = pleaseGetOutputNameForFinalFullNetwork(chromName, parentName, disease, bodyRegion, dataScalingOutputMini, minNumSourcesGeneGRN)
    #outputName = parentName + "fullGeneRegulatoryNetwork_GenieRTNTrenaTrrust2_Scgrnom_TFtoGenes_"  + str(chromName) + "_" + disease + "_" + bodyRegion + dataScalingOutputMini + ".csv"
    
    
    
    
    # please filter this 
    fullNetworkDFForChrom["region"] = repeatVal(bodyRegion, fullNetworkDFForChrom.shape[0]) 
    #fullNetDF[fullNetDF["sources"] != "Trrust2"]
    fullNetworkDFForChrom = fullNetworkDFForChrom[fullNetworkDFForChrom["numSources"] >= minNumSourcesGeneGRN]
    
    
    sourcesToRemoveList = sourcesToRemoveListFunction(minNumSourcesGeneGRN)
    
    if len(sourcesToRemoveList) > 0:
        print(":) please note we are removing these sources: ", sourcesToRemoveList)
        print("beforehand, please note the # of rows is: ", fullNetworkDFForChrom.shape[0])
        for sourceToRemove in sourcesToRemoveList:
            fullNetworkDFForChrom = fullNetworkDFForChrom[fullNetworkDFForChrom["sources"] != sourceToRemove]
    print("after removing those ", len(sourcesToRemoveList), " sources, please note the # of rows is: ", fullNetworkDFForChrom.shape[0])
    fullNetworkDFForChrom[~fullNetworkDFForChrom["sources"].isin(sourcesToRemoveList)]
    fullNetworkDFForChrom = fullNetworkDFForChrom.drop_duplicates()
    
    fullNetworkDFForChrom.to_csv(outputName)
    print(fullNetworkDFForChrom)
    
    currentDate = pleaseGetCurrentDateInfo()

    outputName2 = parentName + "fullGeneRegulatoryNetwork_GenieRTNTrenaTrrust2_Scgrnom_TFtoGenes_"  + str(chromName) + "_" + disease + "_" + bodyRegion + dataScalingOutputMini + "_" + currentDate + "_"+str(minNumSourcesGeneGRN)+"minSources.csv"
    fullNetworkDFForChrom.to_csv(outputName2)
    print("please note the output is: ", outputName2)
    return fullNetworkDFForChrom




######## PLEASE NOTE NEW ADDED:
def pleaseGetFinalFullNetworkDFForAllChromosomes(parentName, disease, bodyRegion, dataScalingOutputMini, minNumSourcesGeneGRN, chromatinRegNetFileNameList, grnResultsCombinedFileNameCSV, geneToEntrezIDMappingPath,
    enhancer_buffer_kbp, promoter_buffer_kbp, progressBarPythonIterations):
    sourcesToRemoveList = sourcesToRemoveListFunction(minNumSourcesGeneGRN)
    finalFullNetworkDF = pd.DataFrame()
    chromeNameList = gettingchromeNameList(chromatinRegNetFileNameList)
    #chromeNameList = ['chr1', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr2', 'chr20', 'chr21', 'chr22', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chrX']
    for i in range(0, len(chromeNameList)):
        chromName = chromeNameList[i]
        print(":) i = ", i, " and chromosome is: ", chromName, " and we are getting FinalFullNetworkForThisGene")
        outputName = pleaseGetOutputNameForFinalFullNetwork(chromName, parentName, disease, bodyRegion, dataScalingOutputMini, minNumSourcesGeneGRN)
        scgrnomOutputFileName = chromatinRegNetFileNameList[i]
        try:
            miniDF = pd.read_csv(outputName)
        except:
            print(":( File not found so we are calculating this fullNetworkDFForChrom")
            miniDF = fullNetworkDFForChromFunction(scgrnomOutputFileName, grnResultsCombinedFileNameCSV, parentName, disease, bodyRegion, dataScalingOutputMini, geneToEntrezIDMappingPath, minNumSourcesGeneGRN, enhancer_buffer_kbp, promoter_buffer_kbp, progressBarPythonIterations = 1000)
        if 'Unnamed: 0' in miniDF.columns:
            miniDF.drop(['Unnamed: 0'], axis = 1)
        finalFullNetworkDF = finalFullNetworkDF.append(miniDF)
    if 'Unnamed: 0' in finalFullNetworkDF.columns.tolist():
        finalFullNetworkDF = finalFullNetworkDF.drop(['Unnamed: 0'], axis = 1)
    finalFullNetworkDF = finalFullNetworkDF.drop_duplicates()
    print("finalFullNetworkDF: ", finalFullNetworkDF.shape[0])
    if len(sourcesToRemoveList) > 0:
        finalFullNetworkDF = finalFullNetworkDF[~finalFullNetworkDF["sources"].isin(sourcesToRemoveList)] # 11659 
    outputName = pleaseGetOutputNameForALLChromosomesFinalFullNetwork(parentName, disease, bodyRegion, dataScalingOutputMini, minNumSourcesGeneGRN)
    print(":) please note the final output for the Full Network GRN for ALL Chromosomes will be written here: ", outputName)
    finalFullNetworkDF.to_csv(outputName)
    return finalFullNetworkDF




def pleaseGetOutputNameForALLChromosomesFinalFullNetwork(parentName, disease, bodyRegion, dataScalingOutputMini, minNumSourcesGeneGRN):
    if minNumSourcesGeneGRN > 1:
        outputName = parentName + "finalFullGeneRegulatoryNetwork_GenieRTNTrenaTrrust2_Scgrnom_TFtoGenes_ALLChromosomes_" + disease + "_" + bodyRegion + "_" + dataScalingOutputMini + "_" + str(minNumSourcesGeneGRN) + "_minSources.csv"

    else:
        outputName = parentName + "finalFullGeneRegulatoryNetwork_GenieRTNTrenaTrrust2_Scgrnom_TFtoGenes_ALLChromosomes_" + disease + "_" + bodyRegion + "_" + dataScalingOutputMini + ".csv"
    
    return outputName

############################################
def pleaseGetMotifbreakRResultsFromSNPs(motifbreakRFinalResultsFilePath, includeTheIndividualTFsInGroup_Python, geneToEntrezIDMappingPath):
    # GWAS Motifbreakr results
    snpsDF = pd.read_csv(motifbreakRFinalResultsFilePath)

    snpsCols = snpsDF.columns.tolist()
    unnamedCols = []
    for colName in snpsCols:
        if "Unnamed" in colName:
            unnamedCols.append(colName)
    unnamedCols
    snpsDF = snpsDF.rename(columns = {"geneSymbol":"TF"})
    snpsDF = snpsDF.drop(unnamedCols, axis = 1)
    snpsDF


    infoDF = addingAdditionalGeneInfo(geneToEntrezIDMappingPath)
    infoDF
    tfMappingInfoDF = infoDF
    tfMappingInfoDF['TFName'] = tfMappingInfoDF['Symbols']
    tfMappingInfoDF['TFEntrezId'] = tfMappingInfoDF['entrezID']
    tfMappingInfoDF = tfMappingInfoDF[['TFName',  'TFEntrezId']]
    tfMappingInfoDF = tfMappingInfoDF.set_index(['TFName'])

    tfMappingInfoDF


    # https://datatofish.com/substring-pandas-dataframe/

    if includeTheIndividualTFsInGroup_Python == "Yes":
        print(":) Please note that since ye selected TRUE for includeTheIndividualTFsInGroup, we will also break apart multiple TFs into individual TFs AS WELL:")
        #extraTFsDF = snpsDF[snpsDF['TF'].str.contains('::')]
        #extraTFsDF = extraTFsDF.reset_index(drop=True)
        
        extraTFsDF = snpsDF = snpsDF[~snpsDF["TF"].isna()] 
        extraTFsDF = snpsDF[snpsDF['TF'].str.contains('::')]
        extraTFsDF

        multipleTFsList = list(set(extraTFsDF["TF"].tolist()))
        multipleTFsList
        print(":) Please note that these are the TF groupings:")
        print(multipleTFsList)
        multTFsToIndividualTFsDict = {}
        for tfGroup in multipleTFsList:
            multTFsToIndividualTFsDict[tfGroup] = tfGroup.split("::")
        multTFsToIndividualTFsDict

        anotherDF = pd.DataFrame()
        for i in range(0, len(multipleTFsList)):
            tfsGroup = multipleTFsList[i]
            listOfTFsInGroup = multTFsToIndividualTFsDict[tfsGroup]
            print(":) i = ", i, " and group is: ", tfsGroup)
            for tf in listOfTFsInGroup:
                snpsForTF = snpsDF[snpsDF["TF"] == tfsGroup]
                snpsForTF["TF"] = snpsForTF["TF"].replace(tfsGroup, tf) #multTFsToIndividualTFsDict[multipleTFsList[i]][0])
                snpsForTF = snpsForTF.reset_index(drop=True)
                snpsForTF
                anotherDF = anotherDF.append(snpsForTF)
        anotherDF

        snpsDF = snpsDF.append(anotherDF)
        snpsDF
    snpsDF = snpsDF.merge(tfMappingInfoDF, left_on='TF', right_on='TFName') #on = "targetGeneName")
    snpsDF
    snpsDF["tfAndChrom"] = snpsDF["chromosome"] + "_" + snpsDF["TFEntrezId"].values.astype(str)
    snpsDF = snpsDF.drop_duplicates()
    snpsDF
    return snpsDF


def do_widen_intervalFullNet(row, enhancer_buffer_kbp, promoter_buffer_kbp):
    enhancer_buffer_basePairs  = enhancer_buffer_kbp * 1000
    promoter_buffer_basePairs = promoter_buffer_kbp * 1000
    new_start = None
    new_end = None
    
    if row["regulatoryRegion"] == "enhancer":
        new_start = int(row["start"]) - enhancer_buffer_basePairs
        new_end = int(row["end"]) + enhancer_buffer_basePairs
    
    elif row["regulatoryRegion"] == "promoter":
        new_start = int(row["start"]) - promoter_buffer_basePairs
        new_end = int(row["end"]) + promoter_buffer_basePairs
      
    elif row["regulatoryRegion"] == "Enhancer":
        new_start = int(row["start"]) - enhancer_buffer_basePairs
        new_end = int(row["end"]) + enhancer_buffer_basePairs
    
    elif row["regulatoryRegion"] == "Promoter":
        new_start = int(row["start"]) - promoter_buffer_basePairs
        new_end = int(row["end"]) + promoter_buffer_basePairs
    return [new_start, new_end]

def pleaseGetMotifbreakRResultsFromSNPsForAChromosome(chromName, motifbreakRFinalResultsFilePath, includeTheIndividualTFsInGroup_Python):
    snpsDF = pleaseGetMotifbreakRResultsFromSNPs(motifbreakRFinalResultsFilePath, includeTheIndividualTFsInGroup_Python)
    snpDFForChrom = snpsDF[snpsDF["chromosome"] == chromName]  
    return snpDFForChrom




def find_overlap_subNetwork(row, startColName, endColName):#, enhancer_buffer_kbp, promoter_buffer_kbp):
    overlap = None
    difference = None
    overlapInfo = None
    if row["regulatoryRegion"].upper() == "ENHANCER":
        if row['snpPosition'] >= row[startColName]:
            if row['snpPosition'] <= row[endColName]:
                overlap = "enhancer"
                difference = 0
                overlapInfo = "enhancer"
            elif row['snpPosition'] > row[endColName]:
                difference = row['snpPosition'] - row[endColName]
                overlap = str(difference) + " bases from enhancer end"
                overlapInfo = "enhancer_with_buffer"

        elif row['snpPosition'] < row[startColName]:
            difference = row[startColName] - row['snpPosition']
            overlap = str(difference) + " bases from enhancer start"
            overlapInfo = "enhancer_with_buffer"

    elif row["regulatoryRegion"].upper() == "PROMOTER":
        if row['snpPosition'] >= row[startColName]:
            if row['snpPosition'] <= row[endColName]:
                overlap = "promoter"
                difference = 0
                overlapInfo = "promoter"

            elif row['snpPosition'] > row[endColName]:
                difference = row['snpPosition'] - row[endColName]
                overlap = str(difference) + " bases from promoter end"
                overlapInfo = "promoter_with_buffer"

        elif row['snpPosition'] < row[startColName]:
            difference = row[startColName] - row['snpPosition']
            overlap = str(difference) + " bases from promoter start"      
            overlapInfo = "promoter_with_buffer"
            
            
    
    return [overlap, difference, overlapInfo] #[new_start, new_end]



#def pleaseGetSubNetworkForChromosome(chromName, motifbreakRFinalResultsFilePath, includeTheIndividualTFsInGroup_Python, scgrnomOutputFileName, grnResultsCombinedFileNameCSV, parentName, disease, bodyRegion, dataScalingOutputMini, geneToEntrezIDMappingPath, minNumSourcesGeneGRN, enhancer_buffer_kbp, 
#promoter_buffer_kbp, startColName = "start_x", endColName = "end_x", progressBarPythonIterations = 1000):
def pleaseGetSubNetworkForChromosome(chromName, motifbreakRFinalResultsFilePath, includeTheIndividualTFsInGroup_Python, grnResultsCombinedFileNameCSV, parentName, disease, bodyRegion, dataScalingOutputMini, geneToEntrezIDMappingPath, minNumSourcesGeneGRN, enhancer_buffer_kbp, 
promoter_buffer_kbp, startColName = "start_x", endColName = "end_x", progressBarPythonIterations = 1000):
    outputName = pleaseGetOutputNameForFinalFullNetwork(chromName, parentName, disease, bodyRegion, dataScalingOutputMini, minNumSourcesGeneGRN)


    fullNetDF = pd.read_csv(outputName) #fullNetworkDFForChrom(scgrnomOutputFileName, grnResultsCombinedFileNameCSV, parentName, disease, bodyRegion, dataScalingOutputMini, geneToEntrezIDMappingPath, minNumSourcesGeneGRN, enhancer_buffer_kbp, promoter_buffer_kbp, progressBarPythonIterations = 1000)
    fullNetDF

    if 'Unnamed: 0' in fullNetDF.columns:
        fullNetDF = fullNetDF.drop(['Unnamed: 0'], axis = 1)

    fullNetDF
        # perhaps another function to read in the fullnetwork after construction once
    snpDFForChrom = pleaseGetMotifbreakRResultsFromSNPsForAChromosome(chromName, motifbreakRFinalResultsFilePath, includeTheIndividualTFsInGroup_Python)
    snpDFForChrom
    subNetworkForChrom = fullNetDF.merge(snpDFForChrom, left_on='tfAndChrom', right_on='tfAndChrom') #on = "targetGeneName")
    subNetworkForChrom = subNetworkForChrom.drop_duplicates()
    subNetworkForChrom = subNetworkForChrom.drop_duplicates(subset=None, keep='first', inplace=False)

    subNetworkForChrom = subNetworkForChrom.rename(columns = {"start_y":"snpPosition", "chromosome_x":"chromosome", "TFEntrezId_x:":"TFEntrezId"})

    if 'end_y' in subNetworkForChrom.columns:
        subNetworkForChrom = subNetworkForChrom.drop(['end_y'], axis = 1)
    if 'TFEntrezId_y' in subNetworkForChrom.columns:
        subNetworkForChrom = subNetworkForChrom.drop(['TFEntrezId_y'], axis = 1)
    if 'chromosome_y' in subNetworkForChrom.columns:
        subNetworkForChrom = subNetworkForChrom.drop(['chromosome_y'], axis = 1)    
    subNetworkForChrom
    subNetworkForChrom = subNetworkForChrom[(subNetworkForChrom['snpPosition'] >= subNetworkForChrom['new_start']) & (subNetworkForChrom['snpPosition'] <= subNetworkForChrom['new_end'])]

    #startColName = "start_x"
    #endColName = "end_x"
    subNetworkForChrom[["overlapInfo", "difference", "overlap"]] = subNetworkForChrom.apply(lambda row: pd.Series(find_overlap_subNetwork(row, startColName, endColName)), 
                                               axis=1)
    
    outputName = pleaseGetOutputNameForSubNetworkForChromosome(chromName, parentName, disease, bodyRegion, dataScalingOutputMini, minNumSourcesGeneGRN)
    
    subNetworkForChrom.to_csv(outputName)
    print(":) please note the outputName is: ", outputName)
    return subNetworkForChrom
    
    
    
    
def pleaseGetRegGeneMappingInfoDF(geneToEntrezIDMappingPath):

    infoDF = addingAdditionalGeneInfo(geneToEntrezIDMappingPath)
    infoDF

    regGeneMappingInfoDF = infoDF
    regGeneMappingInfoDF = regGeneMappingInfoDF[['targetGeneEntrezId']]
    regGeneMappingInfoDF
    return regGeneMappingInfoDF

def pleaseGetTFMappingInfoDF(geneToEntrezIDMappingPath):

    infoDF = addingAdditionalGeneInfo(geneToEntrezIDMappingPath)
    infoDF

    tfMappingInfoDF = infoDF
    tfMappingInfoDF['TFName'] = tfMappingInfoDF['Symbols']
    tfMappingInfoDF['TFEntrezId'] = tfMappingInfoDF['entrezID']
    tfMappingInfoDF = tfMappingInfoDF[['TFName',  'TFEntrezId']]
    tfMappingInfoDF = tfMappingInfoDF.set_index(['TFName'])
    
    return tfMappingInfoDF
    
    
    
def pleaseGetOutputNameForSubNetworkForChromosome(chromName, snpParentName, disease, bodyRegion, dataScalingOutputMini, minNumSourcesGeneGRN):
    print("minNumSourcesGeneGRN:", minNumSourcesGeneGRN)
    if int(minNumSourcesGeneGRN) > 1:
        outputName = snpParentName + "fullSNPSubNetwork_"  + str(chromName) + "_" + disease + "_" + bodyRegion + "_" + dataScalingOutputMini + "_" + str(minNumSourcesGeneGRN) + "_minSources.csv"

    else:
        outputName = snpParentName + "fullSNPSubNetwork_"  + str(chromName) + "_" + disease + "_" + bodyRegion + "_" + dataScalingOutputMini + ".csv"
    
    return outputName


def pleaseGetOutputNameForFinalComboGRNGeneExpression(parentName, disease, bodyRegion, dataScalingOutputMini, minNumSourcesGeneGRN):
    if minNumSourcesGeneGRN > 1:
    #finalComboGeneExpRegNet_GenieRTNTrenaTrrust2_TFtoGene_
        outputName = parentName + "finalComboGeneExpRegNet_GenieRTNTrenaTrrust2_TFtoGene_"  + disease + "_" + bodyRegion +  "_" + dataScalingOutputMini + "_" + str(minNumSourcesGeneGRN) + ".0_minSources.csv"

        #outputName = parentName + "finalComboGeneExpressionRegNetwork_GenieRTNTrenaTrrust2_RegStatus_TFtoGenes" + "_" + disease + "_" + bodyRegion +  "_" + dataScalingOutputMini + "_" + str(minNumSourcesGeneGRN) + "_minSources.csv"
    else:
        outputName = parentName + "finalComboGeneExpRegNet_GenieRTNTrenaTrrust2_TFtoGene_" + disease + "_" + bodyRegion +  "_" + dataScalingOutputMini + ".csv"
    return outputName


def organizingOutputFileForTFToGeneRelationship(grnResultsCombinedFileNameCSV, parentName, disease, bodyRegion, dataScalingOutputMini, geneToEntrezIDMappingPath, minNumSourcesGeneGRN, progressBarPythonIterations):#= 1000):
    df = pd.read_csv(grnResultsCombinedFileNameCSV)
    outputName = pleaseGetOutputNameForFinalComboGRNGeneExpression(parentName, disease, bodyRegion, dataScalingOutputMini, minNumSourcesGeneGRN)
    print(":) please note that the output name for the updated source file is: ", outputName)
    if "TF" in df.columns.tolist():
        tfList = df["TF"].tolist()
    else: 
        tfList = df["currentTF"].tolist()
    if "RegulatedGene" in df.columns.tolist():
        regGeneList = df["RegulatedGene"].tolist()
    elif "targetGene" in df.columns.tolist(): 
        regGeneList = df["targetGene"].tolist()
    
    #regGeneList = df["RegulatedGene"].tolist()
    #sourceList = df["sources"].tolist()#df["Source"].tolist()
    #infoList = df["info"].tolist() # df["Info"].tolist()
    
    
    
    
    
    sourceList = df["Source"].tolist()
    infoList = df["Info"].tolist()
    ##################################
#     if ("currentTF" in df.columns.tolist() &  "targetGene" in df.columns.tolist()):
#         df["combinedName"] = df["currentTF"] + " || " +  df["targetGene"]
#     elif ("currentTF" in df.columns.tolist() &  "RegulatedGene" in df.columns.tolist()):
#         df["combinedName"] = df["currentTF"] + " || " +  df["regulatedGene"]
#     elif ("TF" in df.columns.tolist() &  "targetGene" in df.columns.tolist()):
#         df["combinedName"] = df["TF"] + " || " +  df["targetGene"]
#     elif ("TF" in df.columns.tolist() &  "RegulatedGene" in df.columns.tolist()):
#         df["combinedName"] = df["TF"] + " || " +  df["regulatedGene"]       

#############################################
    comboNameList = df["CombinedName"].tolist()

    tfTGInfo_Dict = {} # please note that the key is combo and the value is list of Info
    tfTGSource_Dict = {} # please note that the key is combo and the value is list of Sources

    for i in range(0, len(tfList)):
        currentTF = tfList[i]
        currentRegGene = regGeneList[i]
        source = sourceList[i]
        info = infoList[i]
        combo = comboNameList[i]
        if i % progressBarPythonIterations == 0:
            print(":) i = ", i)
        if combo in tfTGInfo_Dict.keys():
            currentListInfo = tfTGInfo_Dict[combo]
            currentListInfo.append(info)
            tfTGInfo_Dict[combo] = currentListInfo

            currentListSource = tfTGSource_Dict[combo]
            currentListSource.append(source)
            tfTGSource_Dict[combo] = currentListSource   
        else:
            tfTGInfo_Dict[combo] = [info]
            tfTGSource_Dict[combo] = [source]

    listOfListsToAdd = [] # [tf, reg gene, info, source, # of sources]

    for comboName in tfTGInfo_Dict.keys():
        splitList = comboName.split(" || ")
        currentTF = splitList[0]
        currentRegGene = splitList[1]

        info = tfTGInfo_Dict[comboName]
        infoToAdd = ' || '.join(info)

        source = list(set(tfTGSource_Dict[comboName])) #tfTGSource_Dict[comboName]
        sourceToAdd = ' || '.join(source)

        numSources = len(set(source)) #len(info)

        listOfListsToAdd.append([currentTF, currentRegGene, infoToAdd, sourceToAdd, numSources])

    print(":) Please note a heading:")
    print(listOfListsToAdd[1:5])
    masterDF = pd.DataFrame(listOfListsToAdd, columns = ["currentTF", "regulatedGene", "info", "sources", "numSources"])
    masterDF["region"] = repeatVal(bodyRegion, masterDF.shape[0])
    masterDF = masterDF.drop_duplicates()
    print(":) Please note that the # of rows in the masterDF for TF-Gene relationships: ", masterDF.shape[0])
    print(masterDF.head())
    
    ############
    infoDF = addingAdditionalGeneInfo(geneToEntrezIDMappingPath)
    infoDF

    regGeneMappingInfoDF = infoDF
    regGeneMappingInfoDF = regGeneMappingInfoDF[['targetGeneEntrezId']]
    regGeneMappingInfoDF

    tfMappingInfoDF = infoDF
    tfMappingInfoDF['TFName'] = tfMappingInfoDF['Symbols']
    tfMappingInfoDF['TFEntrezId'] = tfMappingInfoDF['entrezID']
    tfMappingInfoDF = tfMappingInfoDF[['TFName',  'TFEntrezId']]
    tfMappingInfoDF = tfMappingInfoDF.set_index(['TFName'])

    tfMappingInfoDF


    masterDF["targetGeneName"] = masterDF["regulatedGene"].tolist()
    masterDF["TFName"] = masterDF["currentTF"].tolist()

    newmasterDF = masterDF.merge(tfMappingInfoDF, left_on='currentTF', right_on='TFName') #on = "targetGeneName")

    newmasterDF = newmasterDF.merge(regGeneMappingInfoDF, left_on='regulatedGene', right_on='targetGeneName') #on = "targetGeneName")

    newmasterDF["CombinedEntrezName"] = newmasterDF['TFEntrezId'] + ' || ' + newmasterDF['targetGeneEntrezId'] #.agg(' || '.join, axis=1)


    newmasterDF[newmasterDF["currentTF"] == newmasterDF["TFName"]]

    newmasterDF = newmasterDF[newmasterDF["currentTF"] == newmasterDF["TFName"]]
    newmasterDF = newmasterDF[newmasterDF["regulatedGene"] == newmasterDF["targetGeneName"]]
    newmasterDF = newmasterDF.drop(columns = ['targetGeneName', 'TFName'], axis = 1)
    newmasterDF = newmasterDF.drop_duplicates()
    newmasterDF # 1686046  to 1685344 

    #newmasterDF = retrieveEntrezIDs(geneToEntrezIDMappingPath, masterDF, 'currentTF', 'regulatedGene')

    print(":) Please note that the # of rows in the newmasterDF for TF-Gene relationships (After merging for the Entrez IDs is): ", masterDF.shape[0])

    if minNumSourcesGeneGRN > 1:
        newmasterDF = newmasterDF[newmasterDF["numSources"] >= minNumSourcesGeneGRN]
        print(":) after we filter for only edges that are found in at least ", minNumSourcesGeneGRN, " sources, we have: ", newmasterDF.shape[0], " total rows left")

    newmasterDF.to_csv(outputName)
    print(":) please note the output has been written here: ", outputName)

    currentDate = pleaseGetCurrentDateInfo()
    outputName2 = parentName + "finalComboGeneXpressRegNet_GenieRTNTrenaTrrust2_RegStatus_" + disease + "_" + bodyRegion + dataScalingOutputMini + "_" + currentDate + ".csv"
    newmasterDF.to_csv(outputName2)
    print(":) and more specifically here as well: ", outputName2)
    
    return(newmasterDF)


def gettingSNPSubNetwork(chromatinRegNetFileNameList, grnResultsCombinedFileNameCSV, motifbreakRFinalResultsFilePath, parentName, disease, bodyRegion, dataScalingOutputMini, geneToEntrezIDMappingPath, minNumSourcesGeneGRN, enhancer_buffer_kbp, promoter_buffer_kbp, includeTheIndividualTFsInGroup_Python, SNPOutputNameStem, geneAssignmentsModAssignmentsPhenosWithKmeansOutputNameCSV, powerEstimate, progressBarPythonIterations = 1000, filePathForFullNetwork = "", hasFinalFullNetworkForChromosome = "FALSE"):
    combinedSubsForAllChromsDF = pd.DataFrame()
    snpDF = pleaseGetMotifbreakRResultsFromSNPs(motifbreakRFinalResultsFilePath, includeTheIndividualTFsInGroup_Python, geneToEntrezIDMappingPath)

    chromeNameList = gettingchromeNameList(chromatinRegNetFileNameList)
    for i in range(0, len(chromeNameList)):
        chromName = chromeNameList[i]
        print(":) i = ", i, " and chromosome = ", chromName)

        scgrnomOutputFileName = chromatinRegNetFileNameList[i]
        snpForChrom = snpDF[snpDF["chromosome"] == chromName]

        if hasFinalFullNetworkForChromosome == "TRUE":
            fullNetDF = pd.read_csv(filePathForFullNetwork)
        else:
            fullNetDF = fullNetworkDFForChromFunction(scgrnomOutputFileName, grnResultsCombinedFileNameCSV, parentName, disease, bodyRegion, dataScalingOutputMini, geneToEntrezIDMappingPath, minNumSourcesGeneGRN, enhancer_buffer_kbp, promoter_buffer_kbp, progressBarPythonIterations = 1000)

        if fullNetDF.shape[0] == 0:
            print(":( Please note we skipped this ", chromName, " since 0 rows")
            continue

        df_forChrom = fullNetDF[fullNetDF["chromosome"] == chromName]

        if 'Unnamed: 0' in df_forChrom.columns:
            df_forChrom = df_forChrom.drop(['Unnamed: 0'], axis = 1)
        if 'TFEntrezId_y' in df_forChrom.columns:
            df_forChrom = df_forChrom.drop(['TFEntrezId_y'], axis = 1)
        if 'TFEntrezId.1' in df_forChrom.columns:
            df_forChrom = df_forChrom.drop(['TFEntrezId.1'], axis = 1)
        combinedDFForSubNetworkChrom = pd.DataFrame()
        start = 0
        intervalJump = 20000
        end = start + (intervalJump - 1)   
        listOfVals = [[start, end]]
        while end <= df_forChrom.shape[0]:
            print(":) start = ", start, " and end = ", end) 
            start = end + 1
            end = start + (intervalJump - 1)
            if end >= df_forChrom.shape[0]:
                end = df_forChrom.shape[0]
            if start >= end:
                break
            listOfVals.append([start, end])
        listOfVals
        print(":) listOfVals = ", listOfVals) 
        for listVal in listOfVals:
            start = listVal[0]
            end = listVal[1]
            print(":) start = ", start, " and end = ", end)
            finalFullNetworkDF_part = df_forChrom.loc[start:end]
            subNetworkForChromosome = finalFullNetworkDF_part.merge(snpForChrom, left_on='tfAndChrom', right_on='tfAndChrom')
            subNetworkForChromosome = subNetworkForChromosome.drop_duplicates()
            subNetworkForChromosome = subNetworkForChromosome.rename(columns = {"start_y":"snpPosition", "chromosome_x":"chromosome", "TFEntrezId_x:":"TFEntrezId"})

            if 'end_y' in subNetworkForChromosome.columns:
                subNetworkForChromosome = subNetworkForChromosome.drop(['end_y'], axis = 1)
            if 'TFEntrezId_y' in subNetworkForChromosome.columns:
                subNetworkForChromosome = subNetworkForChromosome.drop(['TFEntrezId_y'], axis = 1)
            if 'chromosome_y' in subNetworkForChromosome.columns:
                subNetworkForChromosome = subNetworkForChromosome.drop(['chromosome_y'], axis = 1)    
            subNetworkForChromosome
            subNetworkForChromosome = subNetworkForChromosome.rename(columns={"start_y": "snpPosition"})
            #subNetworkForChromosome = subNetworkForChromosome.drop(['end_y'], axis = 1)
            subNetworkForChromosome = subNetworkForChromosome[(subNetworkForChromosome['snpPosition'] >= subNetworkForChromosome['new_start']) & (subNetworkForChromosome['snpPosition'] <= subNetworkForChromosome['new_end'])]
            subNetworkForChromosome
            startColName = "start_x"
            endColName = "end_x"
            sourcesToRemoveList = sourcesToRemoveListFunction(minNumSourcesGeneGRN)
            if subNetworkForChromosome.shape[0] > 0:

                subNetworkForChromosome[["overlapInfo", "difference", "overlap"]] = subNetworkForChromosome.apply(lambda row: pd.Series(find_overlap_subNetwork(row, startColName, endColName)), 
                                                           axis=1)
                combinedDFForSubNetworkChrom = combinedDFForSubNetworkChrom.append(subNetworkForChromosome)    
                if len(sourcesToRemoveList) > 0:
                    subNetworkForChromosome = subNetworkForChromosome[~subNetworkForChromosome["sources"].isin(sourcesToRemoveList)]
        outputName = SNPOutputNameStem + bodyRegion + "_" + str(chromName) + "_noTrrust2_atLeast" + str(minNumSourcesGeneGRN) + "Sources.csv"
        outputName
        print(chromName, " # of rows: ", combinedDFForSubNetworkChrom.shape[0])
        combinedDFForSubNetworkChrom.to_csv(outputName)
        combinedSubsForAllChromsDF = combinedSubsForAllChromsDF.append(combinedDFForSubNetworkChrom)
    outputNameALL = SNPOutputNameStem + bodyRegion + "_ALLChromosomes_noTrrust2_atLeast" + str(minNumSourcesGeneGRN) + "Sources.csv"
    
    
    #geneModuleAssignmentsDF.to_csv("F://organizedAlzheimers//GeneExpressionPreparation//outputs//LateralTemporalLobe//LateralTemporalLobe_wgcna_with_kmeans_Power10finalCombined_GeneAndModuleInfo_June9.csv")


    combinedSubsForAllChromsDF = combinedSubsForAllChromsDF.rename({"TFEntrezId_x": "TFEntrezId", "TF":"currentTF"}, axis = 1)
    combinedSubsForAllChromsDF = combinedSubsForAllChromsDF.loc[:,~combinedSubsForAllChromsDF.columns.duplicated()]
    if "Refpvalue" in combinedSubsForAllChromsDF.columns.tolist():
        combinedSubsForAllChromsDF = combinedSubsForAllChromsDF.drop(["Refpvalue"], axis = 1)
    if "Altpvalue" in combinedSubsForAllChromsDF.columns.tolist():
        combinedSubsForAllChromsDF = combinedSubsForAllChromsDF.drop(["Altpvalue"], axis = 1)
    combinedSubsForAllChromsDF
    combinedSubsForAllChromsDF.rename(columns = {'regulatedGene':'targetGene'}, inplace = True)
    combinedSubsForAllChromsDF

    geneModuleAssignmentsDF = pd.read_csv(geneAssignmentsModAssignmentsPhenosWithKmeansOutputNameCSV)
    if 'Unnamed: 0' in geneModuleAssignmentsDF.columns.tolist():
        geneModuleAssignmentsDF = geneModuleAssignmentsDF.drop(['Unnamed: 0'], axis = 1)
    geneModuleAssignmentsDF

    combinedSubsForAllChromsDF = combinedSubsForAllChromsDF.merge(geneModuleAssignmentsDF, on = 'targetGene', how = 'left')

    combinedSubsForAllChromsDF["powerUsed"] = repeatVal(powerEstimate, combinedSubsForAllChromsDF.shape[0])

    combinedSubsForAllChromsDF


    combinedSubsForAllChromsDF = combinedSubsForAllChromsDF.rename(columns = {'targetGeneModule_x':'targetGeneModule',
                                                  'targetGeneSimplePositivePhenotypes_x':'targetGeneSimplePositivePhenotypes',
                                                  'targetGeneMoreInfoPosPhenotypes_x':'targetGeneMoreInfoPosPhenotypes',
                                                  'numOfPositivePhenotypesVec_x':'numOfPositivePhenotypesVec',
                                                  'powerUsed_x':'powerUsed',
                                                  'targetGeneModuleSimplePosPhenotypes_x':'targetGeneModuleSimplePosPhenotypes',
                                                  'moreInfoOnModulePosPhenotypes_x':'moreInfoOnModulePosPhenotypes',
                                                  'targetGeneModule_y':'targetGeneModule',
                                                  'targetGeneSimplePositivePhenotypes_y':'targetGeneSimplePositivePhenotypes',
                                                  'targetGeneMoreInfoPosPhenotypes_y':'targetGeneMoreInfoPosPhenotypes',
                                                  'numOfPositivePhenotypesVec_y':'numOfPositivePhenotypesVec',
                                                  'powerUsed_y':'powerUsed_y',
                                                  'targetGeneModuleSimplePosPhenotypes_y':'targetGeneModuleSimplePosPhenotypes',
                                                  'moreInfoOnModulePosPhenotypes_y':'moreInfoOnModulePosPhenotypes',
                                                  'powerUsed_x':'powerUsed'
                                                 })
    combinedSubsForAllChromsDF = combinedSubsForAllChromsDF.loc[:,~combinedSubsForAllChromsDF.columns.duplicated()]
    combinedSubsForAllChromsDF
    
    
    combinedSubsForAllChromsDF.to_csv(outputNameALL) 
    #combinedSubsForAllChromsDF = pd.read_csv(outputNameALL)
    print(":) YAY! WE are all DONE with the SNP Subnetwork!! :)")
    print(combinedSubsForAllChromsDF.shape)
    print(combinedSubsForAllChromsDF.head())
    return combinedSubsForAllChromsDF

def gettingchromeNameList(chromatinRegNetFileNameList):
    #please note that this function helps us get the list of chromosomes is proper order from the scgrnom output:
    chromeNameList = []
    for i in range(0, len(chromatinRegNetFileNameList)):
        chromName = "chr" + chromatinRegNetFileNameList[i].split("chr")[1].split("_")[0]
        chromeNameList.append(chromName)
    print(":) please note the chromosomes are: ", chromeNameList)
    return chromeNameList


def pleaseGetGenesModulesPhenotypesAssignmentsOrganized(powerEstimate, geneAssignmentsWithKmeansOutputNameCSV, updatedModuleTraitPhenotypeDFFilePath, updatedGeneTraitPhenotypeDFFilePath, geneAssignmentsModAssignmentsPhenosWithKmeansOutputNameCSV):
    print(":) Please note that here, we organize the power estimate, gene names, gene modules, gene significant phenotypes, and significant module phenotypes in a neat table")
    geneModuleAssignmentsDF = pd.read_csv(geneAssignmentsWithKmeansOutputNameCSV)
    geneModuleAssignmentsDF.rename(columns = {'Unnamed: 0':'targetGene', 'module':'targetGeneModule'}, inplace = True) 
    geneModuleAssignmentsDF

    #signifModTraitCorrelationOutputNameCSV = "F://organizedAlzheimers//ADSNPhenoOutputs_NEW//Alzheimers_Brain_Mini_LTL_Demo//WGCNA_and_kMeansAfter_Outputs//WGCNAThenKMeans//PhenoCorrs//ModulePhenoCorrs//WGCNAWithKMeanssigned_pow20_StatSignifModulePhenoCorr_Thresh0.05_log2InputOnly.csv"

    modulePhenotypesDF = pd.read_csv(updatedModuleTraitPhenotypeDFFilePath)
    # checking that multiple modules start with ME
    #if modulePhenotypesDF["module"].tolist()[0][0:2] == "ME":
        #if modulePhenotypesDF["module"].tolist()[1][0:2] == "ME":
            #modulePhenotypesDF["module"] = modulePhenotypesDF["module"].apply(lambda row: row.replace("ME", ""))

    if 'Unnamed: 0' in modulePhenotypesDF.columns.tolist():
        modulePhenotypesDF = modulePhenotypesDF.drop(['Unnamed: 0'], axis = 1)
    if 'moduleEigenegene' in modulePhenotypesDF.columns.tolist():
        modulePhenotypesDF = modulePhenotypesDF.drop(['moduleEigenegene'], axis = 1)    
        modulePhenotypesDF = modulePhenotypesDF[['moduleName', 'powerUsed', 'simplePositivePhenotypesVec', 'moreInfoPositivePhenotypesVec']]
    modulePhenotypesDF.rename(columns = {'moduleName':'targetGeneModule','module':'targetGeneModule', 'simplePositivePhenotypesVec':'targetGeneModuleSimplePosPhenotypes', 'moreInfoPositivePhenotypesVec': 'moreInfoOnModulePosPhenotypes'}, inplace = True) 

    modulePhenotypesDF

    genePhenoFilePath = updatedGeneTraitPhenotypeDFFilePath #signifPositiveGeneTraitCorrelationOutputNameCSV

    genePhenoInfoDF = pd.read_csv(genePhenoFilePath)
    genePhenoInfoDF.rename(columns = {'Unnamed: 0':'targetGene'}, inplace = True) 
    genePhenoInfoDF = genePhenoInfoDF[['targetGene', 'simplePositivePhenotypesVec', 'moreInfoPositivePhenotypesVec', 'numOfPositivePhenotypesVec']] #.rename(columns = {'Unnamed: 0':'targetGene', 'module'}, inplace = True) 
    genePhenoInfoDF.rename(columns = {'simplePositivePhenotypesVec':'targetGeneSimplePositivePhenotypes', 'moreInfoPositivePhenotypesVec':'targetGeneMoreInfoPosPhenotypes'}, inplace = True) 

    genePhenoInfoDF

    geneModuleAssignmentsDF = geneModuleAssignmentsDF.merge(genePhenoInfoDF, on = 'targetGene', how = 'left')
    geneModuleAssignmentsDF


    geneModuleAssignmentsDF = geneModuleAssignmentsDF.merge(modulePhenotypesDF, on = 'targetGeneModule', how = 'left')
    geneModuleAssignmentsDF
    geneModuleAssignmentsDF = geneModuleAssignmentsDF.drop_duplicates()
    geneModuleAssignmentsDF.to_csv(geneAssignmentsModAssignmentsPhenosWithKmeansOutputNameCSV)

    print(":) Please note we organized output for genes, modules, phenotypes and assignments here: ", geneAssignmentsModAssignmentsPhenosWithKmeansOutputNameCSV)
    print(geneModuleAssignmentsDF.shape)
    print(geneModuleAssignmentsDF.head())
    return geneModuleAssignmentsDF



def pleaseGetOrganizedFinalGeneRegulatoryNetwork(geneAssignmentsModAssignmentsPhenosWithKmeansOutputNameCSV, finalFullNetworkChromatinAndGeneExpressCSV, parentName, disease, bodyRegion, dataScalingOutputMini, minNumSourcesGeneGRN, chromatinRegNetFileNameList, grnResultsCombinedFileNameCSV, geneToEntrezIDMappingPath,
    enhancer_buffer_kbp, promoter_buffer_kbp, powerEstimate, progressBarPythonIterations, hasFinalFullNetworkForChromosome = "FALSE", filePathForFullNetwork = ""):

    print(":) Please note that this function pleaseGetOrganizedFinalGeneRegulatoryNetwork helps us organize our Full Gene Regulatory Network (combining Chromatin Epigenomics data with gene expression data)")
    geneModuleAssignmentsDF = pd.read_csv(geneAssignmentsModAssignmentsPhenosWithKmeansOutputNameCSV)
    if 'Unnamed: 0' in geneModuleAssignmentsDF.columns.tolist():
        geneModuleAssignmentsDF = geneModuleAssignmentsDF.drop(['Unnamed: 0'], axis = 1)
    geneModuleAssignmentsDF

    if hasFinalFullNetworkForChromosome == "FALSE": # 2821761 
        finalGRN = pleaseGetFinalFullNetworkDFForAllChromosomes(parentName, disease, bodyRegion, dataScalingOutputMini, minNumSourcesGeneGRN, chromatinRegNetFileNameList, grnResultsCombinedFileNameCSV, geneToEntrezIDMappingPath,
        enhancer_buffer_kbp, promoter_buffer_kbp, progressBarPythonIterations)
        finalGRN.rename(columns = {'regulatedGene':'targetGene'}, inplace = True)
        finalGRN


        finalGRN = finalGRN.merge(geneModuleAssignmentsDF, on = 'targetGene', how = 'left')

        finalGRN["powerUsed"] = repeatVal(powerEstimate, finalGRN.shape[0])

        finalGRN


    else:
        finalGRN = pd.read_csv(filePathForFullNetwork)
        finalGRN["region"] = repeatVal(bodyRegion, finalGRN_Hipp.shape[0])

        finalGRN.rename(columns = {'regulatedGene':'targetGene'}, inplace = True)
        finalGRN


        finalGRN = finalGRN.merge(geneModuleAssignmentsDF, on = 'targetGene', how = 'left')

        finalGRN["powerUsed"] = repeatVal(powerEstimate, finalGRN.shape[0])

        finalGRN
        finalGRN.to_csv(filePathForFullNetwork)
    finalGRN.to_csv(finalFullNetworkChromatinAndGeneExpressCSV)
    print(":) Please note we output the final full network (gene regulatory network based on Gene Expression and Chromatin Data) file here: finalFullNetworkChromatinAndGeneExpressCSV = ", finalFullNetworkChromatinAndGeneExpressCSV)
    print(finalGRN.shape)
    print(finalGRN.head())
    return finalGRN
    
    
    
####  FULL SOLUTION:    
print(":) Please note we have finished reading in this function:  pleaseGetOrganizedFinalGeneRegulatoryNetwork")
print(":) Please note we have finished reading in this function:  gettingSNPSubNetwork(chromatinRegNetFileNameList, grnResultsCombinedFileNameCSV, motifbreakRFinalResultsFilePath, parentName, disease, bodyRegion, dataScalingOutputMini, geneToEntrezIDMappingPath, minNumSourcesGeneGRN, enhancer_buffer_kbp, promoter_buffer_kbp, includeTheIndividualTFsInGroup_Python, SNPOutputNameStem, progressBarPythonIterations = 1000, filePathForFullNetwork = '', hasFinalFullNetworkForChromosome = FALSE")
print(":) Please note we have finished reading in this function: organizingOutputFileForTFToGeneRelationship(grnResultsCombinedFileNameCSV, parentName, disease, bodyRegion, dataScalingOutputMini, geneToEntrezIDMappingPath, minNumSourcesGeneGRN, progressBarPythonIterations)!")



 
 
 
 
 