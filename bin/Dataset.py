class ScData:

    def __init__(self, dataPath : str):
        self.inputDataPath = dataPath
        self.cellTypes = []
        if os.paths.isfile(self.inputData):
            self.data = pandas.readtable(self.inputDataPath)
            for cellType in self.data["cell_types"].unique():
                self.cellTypes.append(CellTypeData(self.data[self.data.cell_type == cellType]))





