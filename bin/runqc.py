from Dataset import Dataset
from sys import argv
import os
import configparser

if __name__ == "__main__":

    d = Dataset(argv[1])

    if len(argv) > 3:
        try:
            standardDeviationCuttof = float(argv[3])
            d.runQC(standardDeviationCuttof)
        except:
            print("Incorrect Standard deviation parameter")

    else:
        d.runQC()

    d.writeReport(f"{argv[2]}.html")
