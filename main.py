#!/usr/bin/python
### TO DO:
### -> Regression?
### -> Filter for individuals with > 20% missing data & lcoi with > 20% missing data
### -> Handle missing data
import subprocess
import sys
import argparse

import statistics
import os
import numpy as np
#from sklearn.linear_model import LinearRegression
import time
from statistics import statisticsClass

NUMBER_OF_STATISTICS = 5
DEBUG = 1       ## BOUCHER: Change this to 1 for debuggin mode
OUTPUTFILENAME = "priors.txt"

BASE_PATH = os.path.dirname(os.path.realpath(__file__))

POPULATION_GENERATOR = "/blue/boucher/ishayooseph/RefactorRemote/ONeSAMP/build/OneSamp"
FINAL_R_ANALYSIS = "./r_analysis.R"


#############################################################
## Helper functions
#############################################################
def mean(data):
    n = len(data)
    mean = sum(data) / n
    return mean

def variance(data):
    n = len(data)
    mean = sum(data) / n
    deviations = [(x - mean) ** 2 for x in data]
    variance = sum(deviations) / n
    return variance

def stdev(data):
    import math
    var = variance(data)
    std_dev = math.sqrt(var)
    return std_dev

def normalization(data, mean, variance):
    from operator import truediv
    normalized_data = [0 for a in range(len(str(data)))]
    cnt = 0
    for x in range(len(str(data))):
        normalized_data[cnt] = truediv((x - mean), float(np.sqrt(variance)))
        cnt += 1
    return normalized_data

#############################################################
start_time = time.time()

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
parser.add_argument("--i", type = float, help="Missing data for individuals")
parser.add_argument("--l", type = float, help="Missing data for loci")
parser.add_argument("--o", type = str, help="The File Name")
#'i' for indiv (float) and 'l' for loci (float) [for missing data]
#default: 0.2 for both

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

lowerTheta = 1
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


indivMissing = .2
if (args.i):
    indivMissing = float(args.i)

lociMissing = .2
if (args.l):
    lociMissing = float(args.l)

rangeDuration = "%d,%d" % (lowerDuration, upperDuration)

fileName = "oneSampIn"
if (args.o):
    fileName = str(args.o)
else:
    print("WARNING:main: No filename provided.  Using oneSampIn")

if(DEBUG) :
    print("Start calculation of statistics for input population")

rangeTheta = "%d,%d" % (lowerTheta, upperTheta)

#data = np.array([7, 5, 4, 9, 12, 45])
#normalizedData = normalization(data, mean(data), stdev(data))
#for x in range(len(data)):
#    formatted_float = "{:.2f}".format(normalizedData[x])
#    print(formatted_float)
   # print("the value of data is % d " % "{:.2f}".format(normalizedData[x]))


inputFileStatistics = statisticsClass()
inputFileStatistics.readData(fileName)
#MISSIND DATA FUNCTION IS RUN (filter for indiv and loci)
#Output data to make sure its being removed
inputFileStatistics.filterIndividuals(indivMissing)
inputFileStatistics.filterLoci(lociMissing)
inputFileStatistics.stat1()
inputFileStatistics.newStat4()
inputFileStatistics.stat2()
inputFileStatistics.stat3()
inputFileStatistics.stat4()
inputFileStatistics.stat5()
numLoci = inputFileStatistics.numLoci
sampleSize = inputFileStatistics.sampleSize

textList = [str(inputFileStatistics.stat1), str(inputFileStatistics.stat2), str(inputFileStatistics.stat3), str(inputFileStatistics.stat4), str(inputFileStatistics.stat5)]
with open('/blue/boucher/ishayooseph/inputPopStats','w') as file:
    file.write('\t'.join(textList[0:]) + '\t')
file.close()

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

file = open('/blue/boucher/ishayooseph/RemoteProjects/OneSamp3.0/allPopStats', 'w+')
for x in range(numOneSampTrials) :

    loci = inputFileStatistics.numLoci
    sampleSize = inputFileStatistics.sampleSize
    intermediateFilename = "/blue/boucher/ishayooseph/genePopTiny"

    cmd = "%s -u%d -v%s -rC -l%d -i%d -d%s -s -t1 -b%s -f%d -o1 -p > %s" % (POPULATION_GENERATOR, mutationRate, rangeTheta, loci, sampleSize, rangeDuration, rangeNe, minAlleleFreq, intermediateFilename)
#Change command line, then output the generated NE value to our text file (read in then save and output)

    if(DEBUG) :
        print(cmd)

    returned_value = os.system(cmd)  # returns the exit code in unix

    if returned_value:
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


    textList = []
    textList = [str(refactorFileStatistics.NE_VALUE), str(refactorFileStatistics.stat1), str(refactorFileStatistics.stat2),
                    str(refactorFileStatistics.stat3),
                    str(refactorFileStatistics.stat4), str(refactorFileStatistics.stat5)]
    file.writelines('\t'.join(textList) + '\n')



#if (DEBUG):
    #print("Start calculation of statistics for ALL populations")

###################### Normalization and linear regression


#Comment lines 256-337



res = subprocess.call(["module load R && Rscript /blue/boucher/ishayooseph/RefactorRemote/ONeSAMP/refactor/release/rScript.r < /blue/boucher/ishayooseph/allPopStats"], shell = True)
res



if (returned_value):
    print("ERROR:main: Could not Load R.  FATAL ERROR.")
    exit()

returned_value = os.system(FINAL_R_ANALYSIS)
if (returned_value):
    print("ERROR:main: Could not run R code. Fatal ERROR.")

if (DEBUG):
    print("Finish linear regression")

print("--- %s seconds ---" % (time.time() - start_time))

