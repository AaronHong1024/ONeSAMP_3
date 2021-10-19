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
    durationLower = 2

print(durationLower)
print("__________________")

if (args.uD):
    durationUpper = float(args.uD)
else:
    durationUpper = 8

print(durationUpper)
print("__________________")

if (args.o):
    fileName = str(args.o)
else:
    fileName = "None Provided"
    #Close program?
    #No file no data

print(fileName)
print("__________________")



#Opening file
myfile = open(fileName, "r")
line = myfile.readline()

#Reading number of columns
numCol = 0;
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
totalAlleles = (numCol-1) * numRow
for i in range(numRow):
    for j in (data[i]):
        if (j == '0101' or j == '0202' or j == '0303' or j == '0404'): #What was actg?
                homoCount = homoCount + 1

print('Homozygosity Count: ',homoCount)
print("Total number of alleles", totalAlleles)
stat2 = homoCount/totalAlleles
print ("Stats2 is ", stat2)
#print(arr);

#Stats3
homoDiff = 0
totalHomoDiff = 0
#Total count is same as above stats
for i in range(numRow):
    for j in (data[i]):
        #AA CC TT GG
        if (j == '0101' or j == '0202' or j == '0303' or j == '0404'): #What was actg?
            homoDiff = int(j) - stat2
        totalHomoDiff = totalHomoDiff + homoDiff

#'print(k)' prints 0303 meaning it correctly converts from str to int
stat3 = totalHomoDiff/(totalAlleles-1)
print ("Stats3 is ", stat3)


#Stat5 and stat4

a = 0
c = 0
t = 0
g = 0
num = []
stat5 = 0
for i in range(numCol):
    # Setting to zero
    num *= 0
    addStat5 = 0
    a = 0
    c = 0
    t = 0
    g = 0
    for j in range(numRow):
            #Checking freq of first two numbers
            if (data[j][i][2:] == '01'):
                a = a + 1
            elif (data[j][i][2:] == '02'):
                c = c + 1
            elif (data[j][i][2:] == '03'):
                t = t + 1
            elif (data[j][i][2:] == '04'):
                g = g + 1

            #Checking last two numbers
            if (data[j][i][:2] == '01'):
                a = a + 1
            elif (data[j][i][:2] == '02'):
                c = c + 1
            elif (data[j][i][:2] == '03'):
                t = t + 1
            elif (data[j][i][:2] == '04'):
                g = g + 1

    divisor = numCol * 2
    a = a / divisor
    c = c / divisor
    t = t / divisor
    g = g / divisor

    if (a > 0):
        num.append(float(a))
    if (c > 0):
        num.append(float(c))
    if (t > 0):
        num.append(float(t))
    if (g > 0):
        num.append(float(g))

    if (num):
        for i in num:
            addStat5 = float(addStat5) + float(i * i)
        y = 1 - addStat5
        stat5 = stat5 + y

    #stat 5 divided by number of loci? ie 20 for genepoptinytiny??

print('Stat5 is ', stat5)

#Stat4
stat4 = 0.0
tempstat4 = 0.0
for i in range(numCol):
    # Setting to zero
    num *= 0
    addStat5 = 0.0
    stat4Math = 0.0
    a = 0
    c = 0
    t = 0
    g = 0
    for j in range(numRow):
            #Checking freq of first two numbers
            if (data[j][i][2:] == '01'):
                a = a + 1
            elif (data[j][i][2:] == '02'):
                c = c + 1
            elif (data[j][i][2:] == '03'):
                t = t + 1
            elif (data[j][i][2:] == '04'):
                g = g + 1

            #Checking last two numbers
            if (data[j][i][:2] == '01'):
                a = a + 1
            elif (data[j][i][:2] == '02'):
                c = c + 1
            elif (data[j][i][:2] == '03'):
                t = t + 1
            elif (data[j][i][:2] == '04'):
                g = g + 1

            if ((len(data[j][i])  == 4) and (data[j][i][2:] != data[j][i][:2])):
                stat4Math = stat4Math + 1
            
    divisor = numCol * 2
    a = a / divisor
    c = c / divisor
    t = t / divisor
    g = g / divisor

    if (a > 0):
        num.append(float(a))
    if (c > 0):
        num.append(float(c))
    if (t > 0):
        num.append(float(t))
    if (g > 0):
        num.append(float(g))

    if (num):
        for i in num:
            addStat5 = float(addStat5) + float(i * i)

    if (addStat5 > 0.0  and stat4Math > 0.0):
        addStat5 = addStat5/(numCol*2)
        stat4Math = stat4Math/numCol
        tempstat4 = stat4Math / addStat5
    stat4 = stat4 + tempstat4

#REMEBER to divide by 1-(1/L) ask abt  value of L <<<<<<<<<<<<
print ('Stat4 is ',stat4)
