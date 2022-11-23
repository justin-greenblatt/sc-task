import os
import pandas
from math import sqrt

class CellTypeData:
    def __init__(self, name, inputData, geneStdThreshold, umiStdThreshold, errorStdThreshold):
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
       
        #Perform QC and append results to data
        self.data["ratio_error/std"] = standardizedResidualsY

        self.data["bad_ratio_qc"] = list([1 if abs(yr) > self.params["residual_threshold"] else 0 for yr in standardizedResidualsY])
        self.data["to_many_genes_qc"] = list([1 if y >= self.params["max_genes"] else 0 for y in Y])
                
        self.data["not_enough_genes_qc"] = list([1 if y <= self.params["min_genes"] else 0 for y in Y])

        self.data["to_many_umis_qc"] = list([1 if x >= self.params["max_umis"] else 0 for x in X])
            
        self.data["not_enough_umis_qc"] = list([1 if x <= self.params["min_umis"] else 0 for x in X])

        #Data holding object for this Cell type
        self.data = inputData
        self.data["quality_control"] = filterColumn
        self.data["standardized_residuals"] = self.standardizedResidualsY

    def formatGraph(self):


def formatData(varName, X, Y, rgbArray, markerSize):
    return f"""
var cellData = {{
  x: [1, 2, 3, 4],
  y: [10, 15, 13, 17],
  mode: 'markers',
  type: 'scatter'
}};
"""

def formatLine(varName, p1, p2, rgbArray = (150, 150, 150), lineWidth = 2):
    return f"""
var lower_ngenes = {{
  x: [0, 20],
  y: [3, 3],
  mode: 'lines',
  type: 'scatter',
  name: 'nGenes',
  line: {{
    color: 'rgb(170, 170, 170)',
    width: 2
  }}
}};
"""

       
        out = ""
        out += formatLine("lower_nGenes_threshold", (0, self.params["min_ngenes"]), (self.params["maxX"], self.params["min_ngenes"])))
        out += formatLine("upper_nGenes_threshold", (0, self.params["max_ngenes"]), (self.params["maxX"], self.params["max_ngenes"])))
        out += formatLine("lower_numis_threshold", (self.params["min_numis"], 0), (self.params["min_numis"], self.params["maxY"]))
        out += formatLine("upper_numis_threshold", (self.params["max_numis"], 0), (self.params["max_numis"], self.params["maxY"]))
        out += formatData("low_error", list([umi for umi, key in zip(self.data["nUmis"], self.data["quality_keys"]) if not key.startswith("bad")])
        out += formatData("high_error", list([umi for umi, key in zip(self.data["nUmis"], self.data["quality_keys"]) if key.startswith("bad")])
