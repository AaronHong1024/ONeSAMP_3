### TO DO:
### 0. Statistic 1 fixed?
### 1. Run refactor with correct parameters
### 2. Check we are using all input parameters, remove those that are unused
### 3. Regression?
### 4. Filter for individuals with > 20% missing data & lcoi with > 20% missing data
### 5. Handle missing data
### 6.


import sys
import argparse
import numpy
#import oneSampStatistics   ## WHAT IS THIS?
import os
import numpy as np
from sklearn.linear_model import LinearRegression
import time
start_time = time.time()

DEBUG = 0       ## BOUCHER: Change this to 1 for debuggin mode

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

#est variables
#set defaults for variables and reads in what user provides

minAlleleFreq = 0.005
if (args.m):
    minAlleleFreq = float(args.m)

#print(minAlleleFreq)
#print("__________________")
#print(mutationRate)
#print(lowerNe)
#print("__________________")
#print(upperNe)
#print("__________________")
#print(lowerTheta)
#print("__________________")
#print(upperTheta)
#print("__________________")
#print(numOneSampTrials)
#print("__________________")
#print(durationLower)
#print("__________________")
#print(durationUpper)
#print("__________________")
#print(fileName)
#print("__________________")



mutationRate = 0.000000012
if (args.r):
    mutationRate = float(args.r)

lowerNe = 10
if (args.lNe):
    lowerNe = int(args.lNe)

upperNe = 500
if (args.uNe):
    upperNe = int(args.uNe)

lowerTheta = 0
if (args.lT):
    lowerTheta = float(args.lT)

upperTheta = 10
if (args.uT):
    upperTheta = float(args.uT)

numOneSampTrials = 50000
if (args.s):
    numOneSampTrials = int(args.s)

durationLower = 2
if (args.lD):
    durationLower = float(args.lD)

durationUpper = 8
if (args.uD):
    durationUpper = float(args.uD)

fileName = "oneSampIn"
if (args.o):
    fileName = str(args.o)
else:
    print("WARNING:main: No filename provided.  Using oneSampIn")


if(DEBUG) :
    print("Start calculation of statistics for input population")

currStatistics = oneSampStatistics.statisticsClass()
currStatistics.readData(fileName)
#Something wrong w func stat1
#currStatistics.stat1()
currStatistics.stat2()
currStatistics.stat3()
currStatistics.stat4()
currStatistics.stat5()
numLoci = currStatistics.numLoci
sampleSize = currStatistics.sampleSize


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
        NeVals = 256 # This needs to be fixed
        intermediateFilename = "currRefactorFile"
        cmd = "./refactor.main -t1 -rC -b%s -d1 -u%s -v%s  -s -l%s -i%s -o1 -f%s -p > %s" % (NeVals, mutationRate, lowerTheta, numLoci, sampleSize, minAlleleFreq, intermediateFilename)
        returned_value = os.system(cmd)  # returns the exit code in unix

        if returned_value :
            print("ERROR:main:Refactor did not run")
            exit()

        else :
            refactorFileStatistics = oneSampStatistics.statisticsClass()
            refactorFileStatistics.readData( intermediateFilename )
            # Something wrong w func stat1
            # currStatistics.stat1()
            refactorFileStatistics.stat2()
            refactorFileStatistics.stat3()
            refactorFileStatistics.stat4()
            refactorFileStatistics.stat5()
            #statistics1[x] =  refactorFileStatistics.stat1
            statistics2[x] = refactorFileStatistics.stat2
            statistics3[x] = refactorFileStatistics.stat3
            statistics4[x] = refactorFileStatistics.stat4
            statistics5[x] = refactorFileStatistics.stat5

if (DEBUG):
    print("Start calculation of statistics for ALL populations")

        ############ Boucher: setting things up for linear regression

if (DEBUG):
    print("Start linear regression")


model = LinearRegression()
model = LinearRegression().fit(x, y)


if (DEBUG):
    print("Finish linear regression")

print("--- %s seconds ---" % (time.time() - start_time))

