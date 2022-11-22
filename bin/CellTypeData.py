import os
import pandas
from math import sqrt

class CellTypeData:
    def __init__(self, name, inputData, geneStdThreshold, umiStdThreshold, residualThreshold):
        self.name = name

        #Calculating Mean and std for genes and umi. Also assigning X to Umis and Y to genes
        n = len(inputData["nGenes"])

        X = inputData["nGenes"]
        meanX = sum(X)/n
        stdDeviationX = sqrt(sum(list([(x - meanX) ** 2 for x in X]))/n)
       
        Y = inputData["nUMIs"]
        meanY = sum(Y)/n
        stdDeviationY = sqrt(sum(list([(y - meanY) ** 2 for y in Y]))/n)

        #Calculate linear regression model for cell type data using Least Squares Method.
        #A "blast from the past" of statistics class

        #Calculate slope for linear model
        regressionSlope = ((n * (sum(list([x * y for x,y in zip(X,Y)])))) - (sum(X) * sum(Y))) / \
                          ((n * sum(list([x ** 2 for x in X])) - (sum(X) ** 2)))
        
        #calculate coefficient for the linear model
        regressionCoeff = (sum(Y) - (regressionSlope * sum(X))) / \
                                                   n
        #regression function should be positive as data is has positive correlation 
        if regressionSlope <= 0:
            print(f"Warning! {self.name} distribution is not as expected by a positive linear model")

        #Define Linear regression function
        def predict(x):
            return (regressionSlope * x) + regressionCoeff

        #Calculate standardized residual to be able to use std to classify distance from prediction
        self.standardizedResidualsY = list([(y - predict(x))/stdDeviationY for x,y in zip(X,Y)])
        
        #Define a dictionary "params" that holds the values for cutoff thresholds
        self.params = {}
        self.params["residual_threshold"] = residualThreshold
        self.params["min_genes"] = meanY - (stdDeviationY * geneStdThreshold)
        self.params["max_genes"] = meanY + (stdDeviationY * geneStdThreshold)
        self.params["min_umis"] =  meanX - (stdDeviationX * umiStdThreshold)
        self.params["max_umis"] =  meanX + (stdDeviationX * umiStdThreshold)
        self.params["meanX"] = meanX
        self.params["meanY"] = meanY
        self.params["stdDeviationY"] = stdDeviationY
        self.params["stdDeviationX"] = stdDeviationX
        self.params["regressionSlope"] = regressionSlope
        self.params["regressionCoeff"] = regressionCoeff
        self.params["maxX"] = max(X)
        self.params["maxY"] = max(Y)
        #Iterate over entries and filter them according to the thresholds defined in params. 
        filterColumn = []
        for x,y,yr in zip(X,Y,self.standardizedResidualsY):
            if abs(yr) > self.params["residual_threshold"]:
                filterColumn.append("bad gene/umi ratio")
            elif y >= self.params["max_genes"]:
                filterColumn.append("to many genes")
            elif y <= self.params["min_genes"]:
                filterColumn.append("not enough genes")
            elif x <= self.params["min_umis"]:
                filterColumn.append("not enough UMIs")
            elif x >= self.params["max_umis"]:
                filterColumn.append("to many UMIs")
            else:
                filterColumn.append("ok")

        #Data holding object for this Cell type
        self.data = inputData
        self.data["quality_control"] = filterColumn
        self.data["standardized_residuals"] = self.standardizedResidualsY

    def formatGraph(self)

        def formatLine(varName, p1, p2, rgbArray = (150, 150, 150), lineWidth = 2):
            return f"\nvar {varName} = \{\nx: [{p1[0]} ,{p2[0]}],\ny: [{p1[0]} ,{p2[0]}],\nmode: 'lines',\ntype: 'scatter',\nline: \{\ncolor: 'rgb{rgbArray}',\nwidth: {lineWidth}\n\}\n\};\n"

        def formatData(varName, X, Y, rgbArray, markerSize):
            return f"\nvar {varName} = \{\nx: {X},\ny: {Y},\nmode: 'markers',\ntype: 'scatter',\nmarker : \{\ncolor: 'rgb{rgbArray}',\nsize: {markerSize}\n\}\n\};\n"
        
        out = ""
        out += formatLine("lower_nGenes_threshold", (0, self.params["min_ngenes"]), (self.params["maxX"], self.params["min_ngenes"])))
        out += formatLine("upper_nGenes_threshold", (0, self.params["max_ngenes"]), (self.params["maxX"], self.params["max_ngenes"])))
        out += formatLine("lower_numis_threshold", (self.params["min_numis"], 0), (self.params["min_numis"], self.params["maxY"]))
        out += formatLine("upper_numis_threshold", (self.params["max_numis"], 0), (self.params["max_numis"], self.params["maxY"]))
        out += formatData("low_error", list([umi for umi, key in zip(self.data["nUmis"], self.data["quality_keys"]) if not key.startswith("bad")])
        out += formatData("high_error", list([umi for umi, key in zip(self.data["nUmis"], self.data["quality_keys"]) if key.startswith("bad")])
