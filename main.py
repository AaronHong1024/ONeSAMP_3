#!/usr/bin/python
import subprocess
import sys
import argparse

import os
import numpy as np
import time
from statistics import statisticsClass

NUMBER_OF_STATISTICS = 5
DEBUG = 0  ## BOUCHER: Change this to 1 for debuggin mode
OUTPUTFILENAME = "priors.txt"

BASE_PATH = os.path.dirname(os.path.realpath(__file__))

POPULATION_GENERATOR = "./build/OneSamp"
FINAL_R_ANALYSIS = "./scripts/rScript.r"

def getName(filename):
    (_, filename) = os.path.split(filename)
    return filename

#############################################################
start_time = time.time()

parser = argparse.ArgumentParser()
parser.add_argument("--m", type=float, help="Minimum Allele Frequency")
parser.add_argument("--r", type=float, help="Mutation Rate")
parser.add_argument("--lNe", type=int, help="Lower of Ne Range")
parser.add_argument("--uNe", type=int, help="Upper of Ne Range")
parser.add_argument("--lT", type=float, help="Lower of Theta Range")
parser.add_argument("--uT", type=float, help="Upper of Theta Range")
parser.add_argument("--s", type=int, help="Number of OneSamp Trials")
parser.add_argument("--lD", type=float, help="Lower of Duration Range")
parser.add_argument("--uD", type=float, help="Upper of Duration Range")
parser.add_argument("--i", type=float, help="Missing data for individuals")
parser.add_argument("--l", type=float, help="Missing data for loci")
parser.add_argument("--o", type=str, help="The File Name")

args = parser.parse_args()

#########################################
# INITIALIZING PARAMETERS
#########################################

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

if (int(lowerNe) > int(upperNe)):
    print("ERROR:main:lowerNe > upperNe. Fatal Error")
    exit()

if (int(lowerNe) < 1):
    print("ERROR:main:lowerNe must be a positive value. Fatal Error")
    exit()

if (int(upperNe) < 1):
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

if (DEBUG):
    print("Start calculation of statistics for input population")

rangeTheta = "%d,%d" % (lowerTheta, upperTheta)

#########################################
# STARTING INITIAL POPULATION
#########################################

inputFileStatistics = statisticsClass()

# t = time.time()
inputFileStatistics.readData(fileName)
inputFileStatistics.filterIndividuals(indivMissing)
inputFileStatistics.filterLoci(lociMissing)
inputFileStatistics.test_stat1()
inputFileStatistics.test_stat2()
inputFileStatistics.test_stat3()
inputFileStatistics.test_stat4()
inputFileStatistics.test_stat5()
# print(f'coast:{time.time() - t:.4f}s')
#
# t = time.time()
# inputFileStatistics.testRead(fileName)
# # inputFileStatistics.filterIndividuals(indivMissing)
# # inputFileStatistics.filterLoci(lociMissing)
# inputFileStatistics.new_stat1()
# inputFileStatistics.stat2()
# inputFileStatistics.stat3()
# inputFileStatistics.newStat4()
#
# inputFileStatistics.stat5()
# print(f'coast:{time.time() - t:.4f}s')
numLoci = inputFileStatistics.numLoci
sampleSize = inputFileStatistics.sampleSize

##Creting input file with intial statistics
textList = [str(inputFileStatistics.stat1), str(inputFileStatistics.stat2), str(inputFileStatistics.stat3),
            str(inputFileStatistics.stat4), str(inputFileStatistics.stat5)]
inputPopStats = "inputPopStats_" + getName(fileName)
with open(inputPopStats, 'w') as fileINPUT:
    fileINPUT.write('\t'.join(textList[0:]) + '\t')
fileINPUT.close()

if (DEBUG):
    print("Finish calculation of statistics for input population")

#############################################
# FINISH STATS FOR INITIAL INPUT  POPULATION
############################################

#########################################
# STARTING ALL POPULATIONS
#########################################

if (DEBUG):
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

allPopStats = "allPopStats_" + getName(fileName)
fileALLPOP = open(allPopStats, 'w+')
for x in range(numOneSampTrials):

    loci = inputFileStatistics.numLoci
    sampleSize = inputFileStatistics.sampleSize
    # change the intermediate file name
    intermediateFilename = "intermediate_" + getName(fileName)

    cmd = "%s -u%d -v%s -rC -l%d -i%d -d%s -s -t1 -b%s -f%d -o1 -p > %s" % (
        POPULATION_GENERATOR, mutationRate, rangeTheta, loci, sampleSize, rangeDuration, rangeNe, minAlleleFreq,
        intermediateFilename)

    if (DEBUG):
        print(cmd)

    returned_value = os.system(cmd)

    if returned_value:
        print("ERROR:main:Refactor did not run")
        exit()

    refactorFileStatistics = statisticsClass()

    refactorFileStatistics.readData(intermediateFilename)
    refactorFileStatistics.test_stat1()
    refactorFileStatistics.test_stat2()
    refactorFileStatistics.test_stat3()
    refactorFileStatistics.test_stat4()
    refactorFileStatistics.test_stat5()

    statistics1[x] = refactorFileStatistics.stat1
    statistics2[x] = refactorFileStatistics.stat2
    statistics3[x] = refactorFileStatistics.stat3
    statistics4[x] = refactorFileStatistics.stat4
    statistics5[x] = refactorFileStatistics.stat5

    # Making file with stats from all populations
    textList = []
    textList = [str(refactorFileStatistics.NE_VALUE), str(refactorFileStatistics.stat1),
                str(refactorFileStatistics.stat2),
                str(refactorFileStatistics.stat3),
                str(refactorFileStatistics.stat4), str(refactorFileStatistics.stat5)]
    # print (textList)
    fileALLPOP.write('\t'.join(textList) + '\n')
fileALLPOP.close()

#########################################
# FINISHING ALL POPULATIONS
########################################
# STARTING RSCRIPT
#########################################
# TODO double check there
ALL_POP_STATS_FILE = allPopStats

rScriptCMD = "Rscript %s < %s" % (FINAL_R_ANALYSIS, ALL_POP_STATS_FILE)
print(rScriptCMD)
res = os.system(rScriptCMD)

if (res):
    print("ERROR:main: Could not run Rscript.  FATAL ERROR.")
    exit()

if (DEBUG):
    print("Finish linear regression")

print("--- %s seconds ---" % (time.time() - start_time))

# Deleting temporary files
delete1 = "rm " + inputPopStats
delete_INPUTPOP = os.system(delete1)

delete2 = "rm " + allPopStats
delete_ALLPOP = os.system(delete2)

##########################
# END
#########################
