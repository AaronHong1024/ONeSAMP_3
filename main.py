import sys
import argparse
import numpy
import oneSampStatistics

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

if (args.m):
    minAlleleFreq = float(args.m)
else:
    minAlleleFreq = 0.005

print(minAlleleFreq)
print("__________________")

if (args.r):
    mutationRate = float(args.r)
else:
    mutationRate = 0.000000012

print(mutationRate)
print("__________________")

if (args.lNe):
    lowerNe = int(args.lNe)
else:
    lowerNe = 10

print(lowerNe)
print("__________________")

if (args.uNe):
    upperNe = int(args.uNe)
else:
    upperNe = 500

print(upperNe)
print("__________________")

if (args.lT):
    lowerTheta = float(args.lT)
else:
    lowerTheta = 0

print(lowerTheta)
print("__________________")

if (args.uT):
    upperTheta = float(args.uT)
else:
    upperTheta = 10

print(upperTheta)
print("__________________")

if (args.s):
    numOneSampTrials = int(args.s)
else:
    numOneSampTrials = 50000

print(numOneSampTrials)
print("__________________")

if (args.lD):
    durationLower = float(args.lD)
else:
    durationLower = 2

print(durationLower)
print("__________________")

if (args.uD):
    durationUpper = float(args.uD)
else:
    durationUpper = 8

print(durationUpper)
print("__________________")

fileName = ""
if (args.o):
    fileName = str(args.o)
else:
    fileName = "None Provided"
    #Close program?
    #No file no data

print(fileName)
print("__________________")



currStatistics = oneSampStatistics.statisticsClass()
currStatistics.readData(fileName)
#Something wrong w func stat1
#currStatistics.stat1()
currStatistics.stat2()
currStatistics.stat3()
currStatistics.stat4()
currStatistics.stat5()
