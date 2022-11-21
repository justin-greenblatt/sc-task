import os
import pandas

class CellTypeData:
    def __init__(self, name, data):
        self.name = name
        self.data = data
        n = len(data[0])
        X = data[0]
        Y = data[1]
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
        
        


        

class ScData:

    def __init__(self, dataPath : str):
        self.inputDataPath = dataPath
        self.cellTypes = []
        if os.paths.isfile(self.inputData):
            self.data = pandas.readtable(self.inputDataPath)
            for cellType in self.data["cell_types"].unique():
                self.cellTypes.append(CellTypeData(self.data[self.data.cell_type == cellType]))






