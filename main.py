# IDLE CODE
import sys
import argparse
import numpy


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
    durationLower = 0

print(durationLower)
print("__________________")

if (args.uD):
    durationUpper = float(args.uD)
else:
    durationUpper = 0

print(durationUpper)
print("__________________")

if (args.o):
    fileName = str(args.o)
else:
    fileName = "None Provided"

print(fileName)
print("__________________")



#Opening file
myfile = open(fileName, "r")
line = myfile.readline()

#Reading number of columns
numCol = -1;
while line:
    #print(line)
    if line == "Pop\n":
        break
    line = myfile.readline()
    numCol = numCol + 1

print('The number of columns is:', numCol)

#Reading number of rows
numRow = -1
while line:
    if line == "":
        break
    line = myfile.readline()
    numRow = numRow + 1

print('The number of rows is: ', numRow)
myfile.close()

#Reopening file to create matrix
matrixFile = open(fileName, "r")
line = matrixFile.readline()

#Read until "Pop" in file
while line:
    if line == "Pop\n":
        break
    line = matrixFile.readline()

line = matrixFile.readline() #Reads "Pop"

#Starts creating data matrix
data = []
i = 0;
for i in range(numRow):
    temp = []
    temp = line.split()
    temp.pop(1) #Getting rid of comma in array
    data.append(temp)
    line = matrixFile.readline()

#Stats 2
homoCount = 0
totalCount = 0
for i in range(numRow):
    for j in (data[i]):
        totalCount = totalCount + 1
        if (j == '0101' or j == '0202' or j == '0303' or j == '0404'): #What was actg?
                homoCount = homoCount + 1

print('Homozygosity Count: ',homoCount)
print("Total number of alleles", totalCount)
stat2 = homoCount/totalCount
print ("Stats2 is ", stat2)
#print(arr);

#Stats3
homoDiff = 0
totalHomoDiff = 0
#Total count is same as above stats
for i in range(numRow):
    for j in (data[i]):
        totalCount = totalCount + 1
        #AA CC TT GG
        if (j == '0101' or j == '0202' or j == '0303' or j == '0404'): #What was actg?
            homoDiff = int(j) - stat2
        totalHomoDiff = totalHomoDiff + homoDiff

#'print(k)' prints 0303 meaning it correctly converts from str to int
stat3 = totalHomoDiff/(totalCount-1)
print ("Stats3 is ", stat3)

#Resetting value of data to iterate through columns
#data = data[0]
#print (data)

#Stat5 and stat4

#print(data[0][1][:2])
#iterating by columns
for i in range(numCol):
    for j in range(numRow):
            print(data[j][i])
