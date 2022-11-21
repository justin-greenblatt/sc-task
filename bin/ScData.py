import os
import pandas
from math import sqrt
from . import config

class CellTypeData:
    def __init__(self, name, data):
        self.name = name
        self.data = data

        n = len(data[0])
        self.n = n

        X = data[0]
        self.meanX = sum(X)/n
        self.stdDeviationX = sqrt(sum(list([(x - self.meanX) ** 2 for x in X]))/n)
       
        Y = data[1]
        self.meanY = sum(Y)/n
        self.stdDeviationY = sqrt(sum(list([(y - self.meanY) ** 2 for y in Y]))/n)

        #Calculate linear regression model for cell type data using Least Squares Method.
        #A "blast from the past" of statistics class

        #Calculate slope for linear model
        self.regressionSlope = ((n * (sum(list([x * y for x,y in zip(X,Y)])))) - (sum(X) * sum(Y))) / \
                          ((n * sum(list([x ** 2 for x in X])) - (sum(X) ** 2)))
        #calculate coefficient for the linear model
        self.regressionCoeff = (sum(Y) - (self.regressionSlope * sum(X))) / \
                                                   n
        if self.regressionSlope <= 0:
            print(f"Warning! {self.name} distribution is not as expected by a positive linear model")

        def predict(x):
            return (self.regressionSlope * x) + self.regressionCoeff

        self.standardizedResidualsY = list([(y - predict(x))/self.stdDeviationY for x,y in zip(X,Y)])
        self.filterd = list([1 for x,y,yr in zip(X,Y,self.standardizedResidualsY) if \
                                             abs(yr) < float(config["filter"]["residual_threshold"]) \
                                             and y >= float(config["filter"]["min_genes"]) \ 
                                             and y <= float(config["filter"]["max_genes"]) \ 
                                             and x >= float(config["filter"]["min_umis"])  \
                                             and x <= float(config["filter"]["max_umis"])  \
                              else 0])

class ScData:

    def __init__(self, dataPath : str):
        self.inputDataPath = dataPath
        self.cellTypes = []
        if os.paths.isfile(self.inputData):
            self.data = pandas.readtable(self.inputDataPath)
            for cellType in self.data["cell_types"].unique():
                self.cellTypes.append(CellTypeData(self.data[self.data.cell_type == cellType]))





