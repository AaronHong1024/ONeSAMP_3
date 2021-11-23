### TO DO:
### -> Regression?
### -> Filter for individuals with > 20% missing data & lcoi with > 20% missing data
### -> Handle missing data

import sys
import argparse
import numpy
import os
#import numpy as np
#from sklearn.linear_model import LinearRegression
import time
from statistics import statisticsClass

start_time = time.time()

DEBUG = 0       ## BOUCHER: Change this to 1 for debuggin mode
OUTPUTFILENAME = "priors.txt"
POPULATION_GENERATOR = "./refactor_main"
FINAL_R_ANALYSIS = "./r_analysis.R"

#Creating argument parser

parser = argparse.ArgumentParser()
parser.add_argument("--m", type = float, help="Minimum Allele Frequency")
parser.add_argument("--r", type = float, help="Mutation Rate")
parser.add_argument("--lNe", type = int, help="Lower of Ne Range")
parser.add_argument("--uNe", type = int, help="Upper of Ne Range")
parser.add_argument("--lT", type = float, help="Lower of Theta Range")
parser.add_argument("--uT", type = float, help="Upper of Theta Range")
parser.add_argument("--s", type = int, help="Number of OneSamp Trials")
parser.add_argument("--lD", type = float, help="Lower of Duration Range")
parser.add_argument("--uD", type = float, help="Upper of Duration Range")
parser.add_argument("--o", type = str, help="The File Name")

args = parser.parse_args()

minAlleleFreq = 0.005
if (args.m):
    minAlleleFreq = float(args.m)

mutationRate = 0.000000012
if (args.r):
    mutationRate = float(args.r)

lowerNe = 10
if (args.lNe):
    lowerNe = int(args.lNe)

upperNe = 500
if (args.uNe):
    upperNe = int(args.uNe)

if (int(lowerNe) > int(upperNe)) :
    print("ERROR:main:lowerNe > upperNe. Fatal Error")
    exit()

if(int(lowerNe) < 1) :
    print("ERROR:main:lowerNe must be a positive value. Fatal Error")
    exit()

if(int(upperNe) < 1) :
    print("ERROR:main:upperNe must be a positive value. Fatal Error")
    exit()

rangeNe = "%d,%d" % (lowerNe, upperNe)

lowerTheta = 0
if (args.lT):
    lowerTheta = float(args.lT)

upperTheta = 10
if (args.uT):
    upperTheta = float(args.uT)

rangeTheta = "%d,%d" % (lowerTheta, upperTheta)

numOneSampTrials = 50000
if (args.s):
    numOneSampTrials = int(args.s)

lowerDuration = 2
if (args.lD):
    lowerDuration = float(args.lD)

upperDuration = 8
if (args.uD):
    upperDuration = float(args.uD)

rangeDuration = "%d,%d" % (lowerDuration, upperDuration)

fileName = "oneSampIn"
if (args.o):
    fileName = str(args.o)
else:
    print("WARNING:main: No filename provided.  Using oneSampIn")

if(DEBUG) :
    print("Start calculation of statistics for input population")

rangeTheta = "%d,%d" % (lowerTheta, upperTheta)

inputFileStatistics = statisticsClass()
inputFileStatistics.readData(fileName)
inputFileStatistics.stat1()
inputFileStatistics.stat2()
inputFileStatistics.stat3()
inputFileStatistics.stat4()
inputFileStatistics.stat5()
numLoci = inputFileStatistics.numLoci
sampleSize = inputFileStatistics.sampleSize


if(DEBUG) :
    print("Finish calculation of statistics for input population")


############ Boucher: calculate statistics for numOneSampTrials
############        : call refactor to generate a new population and
############        : calcualte the statistics for new population, store in an array

if(DEBUG) :
    print("Start calculation of statistics for ALL populations")

statistics1 = []
statistics2 = []
statistics3 = []
statistics4 = []
statistics5 = []

statistics1 = [0 for x in range(numOneSampTrials)]
statistics2 = [0 for x in range(numOneSampTrials)]
statistics3 = [0 for x in range(numOneSampTrials)]
statistics4 = [0 for x in range(numOneSampTrials)]
statistics5 = [0 for x in range(numOneSampTrials)]


for x in range(numOneSampTrials) :

    loci = inputFileStatistics.numLoci
    sampleSize = inputFileStatistics.sampleSize
    intermediateFilename = "currRefactorFile"
    cmd = "%s -u%d -v%s -rC -l%d -i%d -d%s -s -t1 -b%s -f%d -o1 -p > %s" % (POPULATION_GENERATOR, mutationRate, rangeTheta, loci, sampleSize, rangeDuration, rangeNe, minAlleleFreq, intermediateFilename)

    if(DEBUG) :
        print(cmd)

    returned_value = os.system(cmd)  # returns the exit code in unix

    if returned_value :
        print("ERROR:main:Refactor did not run")
        exit()

    refactorFileStatistics = statisticsClass()
    refactorFileStatistics.readData(intermediateFilename)
    refactorFileStatistics.stat1()
    refactorFileStatistics.stat2()
    refactorFileStatistics.stat3()
    refactorFileStatistics.stat4()
    refactorFileStatistics.stat5()
    statistics1[x] = refactorFileStatistics.stat1
    statistics2[x] = refactorFileStatistics.stat2
    statistics3[x] = refactorFileStatistics.stat3
    statistics4[x] = refactorFileStatistics.stat4
    statistics5[x] = refactorFileStatistics.stat5


outputFile = open(OUTPUTFILENAME, "w")
for x in range(numOneSampTrials) :
    outputline = "%d %d %d %d %d \n" % (statistics1[x], statistics2[x], statistics3[x], statistics4[x], statistics5[x])
    outputFile.write(outputline)
outputFile.close()

if (DEBUG):
    print("Start calculation of statistics for ALL populations")

if (DEBUG):
    print("Start linear regression")


#model = LinearRegression()
#model = LinearRegression().fit(statistics1, statistics2)

## Complete R analysis on output.
## Assumes the output is in "priors.txt"
returned_value = os.system("module load R")
if (returned_value)
    print("ERROR:main: Could not Load R.  FATAL ERROR.")
    exit()

returned_value = os.system(FINAL_R_ANALYSIS)
if (returned_value)
    print("ERROR:main: Could not run R code. Fatal ERROR.")
if (DEBUG):
    print("Finish linear regression")

print("--- %s seconds ---" % (time.time() - start_time))

