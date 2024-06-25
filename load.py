from joblib import load
from statistics import statisticsClass
import argparse
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
import os

parser = argparse.ArgumentParser()
parser.add_argument("--o", type=str, help="The File Name")
parser.add_argument("--m", type=str, help="Model Name")

args = parser.parse_args()

inputFileStatistics = statisticsClass()

fileName = "oneSampIn"
if (args.o):
    fileName = str(args.o)
else:
    print("WARNING:main: No filename provided.  Using oneSampIn")

modelName = "oneSampIn"
if (args.m):
    modelName = str(args.m)
else:
    print("WARNING:main: No filename provided.  Using oneSampIn")

# t = time.time()
inputFileStatistics.readData(fileName)
inputFileStatistics.filterIndividuals(0.2)
inputFileStatistics.filterLoci(0.2)
#if (args.n):
#    inputFileStatistics.filterMonomorphicLoci()
inputFileStatistics.test_stat1()
inputFileStatistics.test_stat2()
inputFileStatistics.test_stat3()
inputFileStatistics.test_stat5()
inputFileStatistics.test_stat4()

inputStatsList = [str(inputFileStatistics.stat1), str(inputFileStatistics.stat2), str(inputFileStatistics.stat3),
             str(inputFileStatistics.stat4), str(inputFileStatistics.stat5)]

inputStatsList = pd.DataFrame([inputStatsList], columns=['Emean_exhyt', 'Fix_index', 'Mlocus_homozegosity_mean', 'Mlocus_homozegosity_variance', 'Gametic_disequilibrium'])
Z = np.array(inputStatsList[['Emean_exhyt','Fix_index','Mlocus_homozegosity_mean','Mlocus_homozegosity_variance','Gametic_disequilibrium']])
name, ext = os.path.splitext(modelName)

model = load(modelName)

res = model.predict(Z)
print(res)
print("chosen Alpha: ", model.alpha_)
