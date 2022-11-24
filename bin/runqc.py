from Dataset import Dataset
from sys import argv, exit
import os
import configparser

if __name__ == "__main__":

    d = Dataset(argv[1])

    if len(argv) > 3:
        try:
            standardDeviationCuttof = float(argv[3])
            d.runQC(standardDeviationCuttof)
        except:
            #try:
            d.runQC(configPath = argv[3])
            #except:
                
            #    print(f"Incorrect config Path or config file: {argv[3]}\nplease reffer to sc-task/exampleConfig.ini for config format.\nAborting analisys :(")
            #    exit(1)

    else:
        d.runQC()

    d.writeReport(f"{argv[2]}.html")
