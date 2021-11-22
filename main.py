### TO DO:
### -> Run refactor with correct parameters
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

rangeNe = "%d,%d" % (lowerNe, upperNe)
#print(rangeNe)

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

#    cmd = "./refactor_main -t1 -rC -b%s -d1 -u%s -v%s  -s -l%s -i%s -o1 -f%s -p > %s" % (NeVals, mutationRate, lowerTheta, numLoci, sampleSize, minAlleleFreq, intermediateFilename)

rangeTheta = "%d,%d" % (lowerTheta, upperTheta)

loci = 10
sampleSize = 200
intermediateFilename = "currRefactorFile"
cmd = "./refactor_main -u%d -v%s -rC -l%d -i%d -d%s -s -t1 -b%s -f%d -o1 -p > %s" % (mutationRate, rangeTheta, loci, sampleSize, rangeDuration, rangeNe, minAlleleFreq, intermediateFilename)   #  - b$reducedSize   - f$ONESAMP2COAL_MINALLELEFREQUENCY - o1 - g < $OUTPUT > $OUTPUT$suffix
print(cmd)

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

    #    cmd = "./refactor_main -t1 -rC -b%s -d1 -u%s -v%s  -s -l%s -i%s -o1 -f%s -p > %s" % (NeVals, mutationRate, lowerTheta, numLoci, sampleSize, minAlleleFreq, intermediateFilename)

     #   cmd = "-u%d" % (rangeDuration,)$mutationRate - v$theta - rC - l$loci - i$outputSampleSize - b$reducedSize - d$duration -$microsatsOrSNPs - t1 - f$ONESAMP2COAL_MINALLELEFREQUENCY - o1 - g < $OUTPUT > $OUTPUT$suffix

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

if (DEBUG):
    print("Start calculation of statistics for ALL populations")

        ############ Boucher: setting things up for linear regression

if (DEBUG):
    print("Start linear regression")


#model = LinearRegression()
#model = LinearRegression().fit(statistics1, statistics2)


if (DEBUG):
    print("Finish linear regression")

print("--- %s seconds ---" % (time.time() - start_time))

