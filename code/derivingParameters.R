

fullDiseaseName = paste0(disease, "Disease")

######### Please note the function
if (log2transformInputData == "TRUE"){
  if (scaleInputData == "TRUE"){
    outputAddOn = "_log2TransformedAndThenScaledInput"
    print(":) Please note we will apply a Log2(x+1) transformation on the input gene expression data and then will scale that.")
    print(":) This final data set will be our input for WGCNA")
    dataScaling = "log2TransformedAndScaledGeneExpression" #  <-- will come from WGCNA
    
  } else {
    outputAddOn = "_log2TransformedInputOnly"
    print(":) Please note we will apply a Log2(x+1) transformation on the input gene expression data.")
    print(":) This final data set will be our input for WGCNA")
    dataScaling = "log2TransformedGeneExpressionOnly" #  <-- will come from WGCNA
    
  }
} else if (scaleInputData == "TRUE"){
  outputAddOn = "_ScaledInputOnly"
  print(":) Please note we will apply a Scale transformation on the input gene expression data.")
  print(":) This final data set will be our input for WGCNA")
  dataScaling = "ScaledGeneExpressionOnly" #  <-- will come from WGCNA
  
} else {
  outputAddOn = "_OriginalInputData"
  print(":) Please note we will use our original input gene expression dataset for WGCNA.")
  dataScaling = "originalGeneExpression" #  <-- will come from WGCNA
  
}
