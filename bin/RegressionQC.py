import os
import pandas
from math import sqrt
from Graph import Graph

class RegressionQC:
    def __init__(self, name, inputData, geneStdThreshold, umiStdThreshold, errorStdThreshold):
        self.name = name
        self.data = inputData
        #Calculating Mean and std for genes and umi. Also assigning X to Umis and Y to genes
        n = len(inputData["nGenes"])

        X = inputData["nUMIs"]
        meanX = sum(X)/n
        stdDeviationX = sqrt(sum(list([(x - meanX) ** 2 for x in X]))/n)
       
        Y = inputData["nGenes"]
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
        standardizedResidualsY = list([(y - predict(x))/stdDeviationY for x,y in zip(X,Y)])
        
        #Define a dictionary "params" that holds the values for cutoff thresholds
        self.params = {}
        self.params["residual_threshold"] = errorStdThreshold
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

    #Define Linear regression function
    def predict(self,x):
        return (self.params["regressionSlope"] * x) + self.params["regressionCoeff"]


    def getMetadata(self):
        return self.params

    def getGraph(self):
        g = Graph(self.name)
        g.addLine("lower_nGenes_threshold", (0, self.params["min_genes"]), (self.params["maxX"], self.params["min_genes"]))
        g.addLine("upper_nGenes_threshold", (0, self.params["max_genes"]), (self.params["maxX"], self.params["max_genes"]))
        g.addLine("lower_numis_threshold", (self.params["min_umis"], 0), (self.params["min_umis"], self.params["maxY"]))
        g.addLine("upper_numis_threshold", (self.params["max_umis"], 0), (self.params["max_umis"], self.params["maxY"]))
        g.addLine("linea_regression", (0, self.params["regressionCoeff"]), (self.params["maxX"], self.predict(self.params["maxX"])), (100,100,250),4)
        g.addData("low_error", list([umi for umi, k in zip(self.data["nUMIs"], self.data["bad_ratio_qc"]) if k == 0]),\
                               list([gene for gene, k in zip(self.data["nGenes"], self.data["bad_ratio_qc"]) if k == 0]))
        g.addData("high_error", list([umi for umi, k in zip(self.data["nUMIs"], self.data["bad_ratio_qc"]) if k == 1]),\
                                list([gene for gene, k in zip(self.data["nGenes"], self.data["bad_ratio_qc"]) if k == 1]),(250, 80,80))
        g.addMetadata(self.getMetadata())
        return g
