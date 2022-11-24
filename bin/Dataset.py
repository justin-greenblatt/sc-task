from RegressionQC import RegressionQC
from Report import Report
class Dataset:
    def __init__(self, tsvPath, separator = '\t'):
        self.name = tsvPath.split('/')[-1].split('.')[0]
        handler = open(tsvPath, 'r')
        #loading data to python
        rawData = list([l.strip('\n').split(separator) for l in handler])  
        handler.close()
        #looking at headers
        headers = {h:n for n,h in enumerate(rawData[0])}
        #finding different cell types in the dataset
        cells = set({c[headers["cell_type"]] for c in rawData[1:]})
        #separating Dataset according to cells
        numerical = ['nGenes', 'nUMIs']
        self.data = dict({c : dict({h : list([int(line[n]) if h in numerical else line[n] for line in rawData[1:] if line[headers["cell_type"]] == c]) for h,n in headers.items()}) for c in cells})
        self.analisys = list([RegressionQC(k, self.data[k]) for k in self.data])

    def writeReport(self, reportPath):
        self.report = Report(self.name)
        for cell in self.analisys:
            self.report.addGraph(cell.getGraph())
        self.report.write(reportPath)
