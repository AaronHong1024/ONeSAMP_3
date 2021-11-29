### TO DO:
### -> Regression?
### -> Filter for individuals with > 20% missing data & lcoi with > 20% missing data
### -> Handle missing data
import sys
import argparse

import statistics
import numpy
import os
import numpy as np
#from sklearn.linear_model import LinearRegression
import time
from statistics import statisticsClass

NUMBER_OF_STATISTICS = 5
DEBUG = 0       ## BOUCHER: Change this to 1 for debuggin mode
OUTPUTFILENAME = "priors.txt"
POPULATION_GENERATOR = "./refactor_main"
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
    normalized_data = [0 for a in range(len(data))];
    cnt = 0;
    for x in data:
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

#data = np.array([7, 5, 4, 9, 12, 45])
#normalizedData = normalization(data, mean(data), stdev(data))
#for x in range(len(data)):
#    formatted_float = "{:.2f}".format(normalizedData[x])
#    print(formatted_float)
   # print("the value of data is % d " % "{:.2f}".format(normalizedData[x]))


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


if (DEBUG):
    print("Start calculation of statistics for ALL populations")

###################### Normalization and linear regression

normalizedStatistics1 = normalization(statistics1, mean(statistics1), stdev(statistics1))
normalizedStatistics2 = normalization(statistics2, mean(statistics2), stdev(statistics2))
normalizedStatistics3 = normalization(statistics3, mean(statistics3), stdev(statistics3))
normalizedStatistics4 = normalization(statistics4, mean(statistics4), stdev(statistics4))
normalizedStatistics5 = normalization(statistics5, mean(statistics5), stdev(statistics5))

normalizedInputFileStatistic1 = normalization(inputFileStatistics.stat1, mean(statistics1), stdev(statistics1))
normalizedInputFileStatistic2 = normalization(inputFileStatistics.stat2, mean(statistics2), stdev(statistics2))
normalizedInputFileStatistic3 = normalization(inputFileStatistics.stat3, mean(statistics3), stdev(statistics3))
normalizedInputFileStatistic4 = normalization(inputFileStatistics.stat4, mean(statistics4), stdev(statistics4))
normalizedInputFileStatistic5 = normalization(inputFileStatistics.stat5, mean(statistics5), stdev(statistics5))

#dist < - sqrt((scaled.sumstat[, 1]-target.s[1]) ^ 2 +
#              (scaled.sumstat[, 2]-target.s[2]) ^ 2 +
#              (scaled.sumstat[, 3]-target.s[3]) ^ 2 +
#              (scaled.sumstat[, 4]-target.s[4]) ^ 2 +
#              (scaled.sumstat[, 5]-target.s[5]) ^ 2 +
#              (scaled.sumstat[, 6]-target.s[6]) ^ 2 +
#              (scaled.sumstat[, 7]-target.s[7]) ^ 2 +
#              (scaled.sumstat[, 8]-target.s[8]) ^ 2)

length = len(data)
dist = [0 for a in length];
distTransform  = [0 for a in length];
for i in range(length):
    dist[i] = np.sqrt(normalizedStatistics1[i] - normalizedInputFileStatistic1) \
                + np.sqrt(normalizedStatistics2[i] - normalizedInputFileStatistic2) \
                + np.sqrt(normalizedStatistics3[i] - normalizedInputFileStatistic3) \
                + np.sqrt(normalizedStatistics4[i] - normalizedInputFileStatistic4) \
                + np.sqrt(normalizedStatistics5[i] - normalizedInputFileStatistic5)
    #dist[!gwt] <- floor(max(dist[gwt])+10)
    #distTransform[i] = np.floor()

# abstol <- quantile(dist,tol)
abstol =  10 #numpy.quantile(dist, 0.01)

#wt1 <- dist < abstol
wt1 = [0 for a in length]
for i in range(length):
    if(dist[i] < abstol) :
        wt1 = 1


#regwt < - 1 - dist[wt1] ^ 2 / abstol ^ 2
#x1 < - scaled.sumstat[, 1][wt1]
#x2 < - scaled.sumstat[, 2][wt1]
#x3 < - scaled.sumstat[, 3][wt1]
#x4 < - scaled.sumstat[, 4][wt1]
regwt = [0 for a in length]
x1 = [0 for a in length]
x2 = [0 for a in length]
x3 = [0 for a in length]
x4 = [0 for a in length]
x5 = [0 for a in length]
for i in range(length):
    if(wt1[i]) :
        regwt[i] = 1 - np.sqrt(dist[i]) / np.sqrt(abstol)
        x1.insert(normalizedStatistics1[i])
        x2.insert(normalizedStatistics2[i])
        x3.insert(normalizedStatistics3[i])
        x4.insert(normalizedStatistics4[i])
        x5.insert(normalizedStatistics5[i])

## ADD Ne Value
outputFile = open(OUTPUTFILENAME, "w")
for x in range(length):
    if(wt1[x]) :
        outputline = "%d %d %d %d %d %d \n" % ( regwt[x], statistics1[x], statistics2[x], statistics3[x], statistics4[x], statistics5[x])
        outputFile.write(outputline)
outputFile.close()


#predmean <- predict.lm(fit1,data.frame(
#x1=target.s[1],x2=target.s[2],x3=target.s[3],x4=target.s[4]))


#model = LinearRegression()
#model = LinearRegression().fit(statistics1, statistics2)

## Complete R analysis on output.
## Assumes the output is in "priors.txt"
returned_value = os.system("module load R")
if (returned_value):
    print("ERROR:main: Could not Load R.  FATAL ERROR.")
    exit()

#returned_value = os.system(FINAL_R_ANALYSIS)
if (returned_value):
    print("ERROR:main: Could not run R code. Fatal ERROR.")

if (DEBUG):
    print("Finish linear regression")

print("--- %s seconds ---" % (time.time() - start_time))

