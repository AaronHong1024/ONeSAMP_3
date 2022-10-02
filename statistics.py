#!/usr/bin/python

import math
import os
import re

import numpy as np
import pandas as pd
from decimal import *


class statisticsClass:
    ####### Members of Class

    ARRAY_MISSINGIndiv = []
    ARRAY_MISSINGLoci = []
    data = []  ## Matrix is [individuals by loci]
    stat1 = 0
    stat2 = 0
    stat3 = 0
    stat4 = 0
    stat5 = 0
    numLoci = 0
    sampleSize = 0  ##Indivduals
    NE_VALUE = 0
    DEBUG = 0
    result = []
    allcnt = []

    ######################################################################
    # writeStatistics                                                   ##
    ######################################################################
    # def writeStatistics(self, myFileName):
    # outputFile = open(myFileName, "w")

    ######################################################################
    # readData                                                          ##
    ######################################################################
    def readData(self, myFileName):
        print("filename: ", myFileName)
        matrixFile = open(myFileName, "r")
        line = matrixFile.readline()
        sampleSize = 0

        # Read until "Pop" in file
        popReached = 0;
        while line:
            if ((line == "Pop\n") or (line == "POP\n") or (line == "pop\n")):
                popReached = 1
                sampleSize = sampleSize - 1
                break
            sampleSize = sampleSize + 1
            line = matrixFile.readline()

        if (popReached == 0):
            print("ERROR:statistics.py:line 32:: POP not contained in file. Fatal error")
            exit()

        # Starts creating data matrix
        data = self.data
        # df = pd.DataFrame(data)
        # print(df)
        # print(data)
        temp = []
        temp = line.split()
        line = matrixFile.readline()
        temp = line.split()
        temp.pop(1)
        currLociCnt = len(temp)  # Getting length of first line of data
        print("currLociCnt: ", currLociCnt)
        for i in range(sampleSize):
            line = re.sub(",", "", line)
            # print("line: ", line)
            temp = line.split()
            # print("temp: ",temp)
            # temp.pop(1)  # Getting rid of comma in array
            data.append(temp)
            line = matrixFile.readline()
            if len(temp) != currLociCnt:
                print(len(temp))
                # print(currLociCnt)
                print("ERROR:statistics.py:line 56:: There is an incorrect number of loci. Fatal Error.")
                exit()

        # Reading in NE Value
        NE_VALUEtemp = 0
        if (line):
            temp = line.split()
            NE_VALUEtemp = temp[0]

        self.data = data
        self.numLoci = currLociCnt - 1  # Subracting 1 bc number takes into acct the indiviudal's number
        self.sampleSize = len(self.data)
        # print("data: ", self.data)
        # print("sample size: ", self.sampleSize)
        self.NE_VALUE = NE_VALUEtemp
        print("NE value: ", NE_VALUEtemp)
        print("-----------------------------------------------------------------")

    def testRead(self, myFileName):
        print("test read filename: ", myFileName)
        matrixFile = open(myFileName, "r")
        line = matrixFile.readline()
        sampleSize = 0

        # Read until "Pop" in file
        popReached = 0;
        while line:
            if ((line == "Pop\n") or (line == "POP\n") or (line == "pop\n")):
                popReached = 1
                sampleSize = sampleSize - 1
                break
            sampleSize = sampleSize + 1
            line = matrixFile.readline()
        print(sampleSize)
        if (popReached == 0):
            print("ERROR:statistics.py:line 32:: POP not contained in file. Fatal error")
            exit()
        # Import the data matrix
        test = []
        last_line = self.get_file_last_line(os.path.abspath(myFileName))
        skipfooter = 0
        if "," not in last_line:
            skipfooter = 1

        test = last_line.strip().split(" ")

        # data = pd.read_csv(myFileName, skiprows=sampleSize + 2, header=None, sep=' , | ', dtype='str', engine='python',
        #                    keep_default_na=False, error_bad_lines=False, warn_bad_lines=False)
        data = pd.read_csv(myFileName, skiprows=sampleSize + 2, header=None, sep=' , | ', dtype='str', engine='python',
                           keep_default_na=False, skipfooter=skipfooter, warn_bad_lines=False)
        print(data)

        print(test)
        NE_VALUEtemp = 0
        if len(test) != data.shape[1] + 1:
            NE_VALUEtemp = test[0]

        self.numLoci = data.shape[1] - 1

        self.NE_VALUE = NE_VALUEtemp
        self.sampleSize = data.shape[0]
        self.data = data
        #print(data)
        # print("-----------------------------------------------------------------")

    def get_file_last_line(self, inputfile):
        with open(inputfile, 'rb') as f:
            first_line = f.readline()
            offset = -50
            while True:
                f.seek(offset, 2)
                lines = f.readlines()
                if len(lines) >= 2:
                    last_line = lines[-1]
                    break
                offset *= 2

            return last_line.decode()

    ######################################################################
    # filterIndividuals                                                 ##
    # Filters the data for all individuals that have > 20% missing data ##
    ######################################################################
    def filterIndividuals(self, PERCENT_MISSINGIndiv):
        for i in range(self.sampleSize):
            individual = self.data.loc[i]
            numMissing = 0
            for j in (individual):
                if (j == '0100' or j == '0001' or j == '0000'):
                    numMissing += 1
            if (numMissing > PERCENT_MISSINGIndiv * self.numLoci):
                print("Deleted:", self.data[j])
                self.ARRAY_MISSINGIndiv.append(self.data[j])
                del self.data[j]
                self.numLoci = self.numLoci - 1
                self.sampleSize = len(self.data) - 1

                # only that place., but double check with David

    def testfilerIndividuals(self, PERCENT_MISSINGIndiv):
        data = self.result
        for row in data.itertuples():
            num = row.count("0100") + row.count("0001") + row.count("0000")
            if num > PERCENT_MISSINGIndiv * self.numLoci:
                self.ARRAY_MISSINGIndiv.append(row)
                print(data.drop(row.Index))
                self.sampleSize = data.shape[0]
            # print(row.Index)

    ######################################################################
    # filterLoci                                                        ##
    ######################################################################
    def filterLoci(self, PERCENT_MISSINGLoci):
        for i in range(self.numLoci):
            individual = self.data[i]
            numMissing = 0
            for j in (individual):
                if (j == '0100' or j == '0001' or j == '0000'):
                    numMissing += 1
            if (numMissing > PERCENT_MISSINGLoci * self.numLoci):
                print("Deleted:", self.data[j])
                self.ARRAY_MISSINGLoci.append(self.data[j])
                del self.data[j]
                self.numLoci = self.numLoci - 1
                self.sampleSize = len(self.data) - 1

 ######################################################################
    # stat1_v2 BW Estimator REVISION                                          ##
    ######################################################################
    def stat1_v2(self):
        if (self.DEBUG):
            print("printing for stat1 begin:")

        data = self.data
        numLoci = len(alleleA[0])
        sampleSize = self.sampleSize

        alleleA = []  # same dimensions as data but holds the first allele in a locus
        alleleB = []  # same dimensions as data but holds the second allele in a locus


       for i in range(len(data)):  # fills up alleleA and alleleB
            temA = []
            temB = []
            for j in range(1, len(data[i])):
                temA.append(data[i][j][:2])
                temB.append(data[i][j][2:])

            alleleA.append(temA)
            alleleB.append(temB)


        allcnt = []  # a 1D array, each element is a dictionary of alleles with frequency counts per loci
    #    homoloci = []  # maintains frequency counts of homologous alleles per loci
        denominator = float(numLoci)

        di = []  # a 1D array that holds the departures of each loci from Hardy-Weinberg equilibrium
        totspots = float(len(data) * 2)  # total number of alleles per locus. used in computing allele freq per locus

        for i in sampleSize:  # fills up allcnt
            allele_counts = [0,0,0,0]
            homozygous_count = 0
            for j in numLoci:

                # check if it is homozygous
                if (alleleA[i][j] == alleleB[i][j]):
                    homozygous_count += 1

                allele_counts[int(alleleA[i][j]) - 1] +=1
                allele_counts[int(alleleB[i][j]) - 1] +=1

                #if (alleleA[i][j] in newdic):
                 #   newdic[alleleA[i][j]] += 1
                #else:
                 #   newdic[alleleA[i][j]] = 1

                #if (alleleB[i][j] in newdic):
                 #   newdic[alleleB[i][j]] += 1
                #else:
                 #   newdic[alleleB[i][j]] = 1
                # END INNER FOR LOOP

            currHomoLoci = (homozygous_count/ denominator)
            currDeparture = ( numLoci -homozygous_count ) / totspots
           # homoloci.append( currHomoLoci )
            allcnt.append(allele_counts)

            # vals = []
            # for key, value in allcnt[i].items():
            #    vals.append(float(value) / float(totspots))
            di.append( (currHomoLoci  *  currHomoLoci) - (currDeparture  *  currDeparture ))
        # END OUTER FOR LOOP
        ###############################

        running_sum = 0
        ai = 0
        alA = 0
        # Loop through Allele A
        for i in sampleSize:

            for l in allcnt[i]:
                if (allcnt[i][l] > 0 )
                    alA = allcnt[i][l]
                    ai = allcnt[i][l] / totspots
                    break

            # Loop through Allele B
            for j in sampleSize:

                keysj = []
                for l in allcnt[j]:
                    if (allcnt[j][l] > 0)
                        keysj.append(allcnt[j][l])

                if (len(keysj) >= 2):

                    alB = keysj[1]
                    bj = allcnt[j][alB] / totspots
                    hits = 0

                    for k in numLoci :
                        if ((alleleA[i][k] == alA or alleleB[i][k] == alA) and (alleleA[j][k] == alB or alleleB[j][k] == alB )):
                            hits += 1

                    x = (hits / (sampleSize - ai *  bj)) /  ((ai * (1 - ai) + di[i]) * (bj * (1 - bj) + di[j]))
                    running_sum += math.sqrt(abs(x))

        self.stat1 = 2*running_sum/(numloci*(numloci-1))

        if (self.DEBUG):
            print("printing for stat1 end   ---->", self.stat1)

    ######################################################################
    # stat1 BW Estimator                                                ##
    ######################################################################
    # def test_stat1(self):
    #     if (self.DEBUG):
    #         print("printing for stat1 begin:")
    #     data = self.result
    #     alleleA = []  # same dimensions as data but holds the first allele in a locus
    #     alleleB = []  # same dimensions as data but holds the second allele in a locus
    #
    #     # fills up alleleA and alleleB
    #     matrix = np.array(data)
    #     a = matrix[:, 1:]
    #
    #     print(a)

    def new_stat1(self):
        if self.DEBUG:
            print("printing for stat1 begin: ")

        data = self.data
        alleleA = []
        alleleB = []

        for row in data.itertuples(index=False):
            temA = []
            temB = []
            for col in range(1, len(row)):
                temA.append(row[col][0:2])
                temB.append(row[col][2:])

            alleleA.append(temA)
            alleleB.append(temB)
        print("alleleA: ", alleleA)
        print("alleleB: ", alleleB)
        allcnt = []
        homoloci = []

        for i in range(len(alleleA[0])):
            newdic = {}
            temp = 0
            for j in range(len(alleleA)):
                if alleleA[j][i] == alleleB[j][i] and alleleA[j][i] == alleleA[0][i]:
                    temp += 1
                if alleleA[j][i] in newdic:
                    newdic[alleleA[j][i]] += 1
                else:
                    newdic[alleleA[j][i]] = 1
                if alleleB[j][i] in newdic:
                    newdic[alleleB[j][i]] += 1
                else:
                    newdic[alleleB[j][i]] = 1

            homoloci.append(float(temp) / float(self.sampleSize))
            allcnt.append(newdic)

        print("homoloci: ", homoloci)
        print("allcnt: ", allcnt)
        self.allcnt = allcnt
        di = []  # a 1D array that holds the departures of each loci from Hardy-Weinberg equilibrium
        totspots = float(
            self.sampleSize * 2)  # total number of alleles per locus. used in computing allele freq per locus
        # print("totspots: ", totspots)
        for i in range(len(allcnt)):  # fills up di
            vals = []
            for key, value in allcnt[i].items():
                vals.append(float(value) / float(totspots))
            print("vals: ", vals)
            di.append(homoloci[i] - math.pow(vals[0], 2))

        print("di: ", di)
        hits = 0

        running_sum = 0

        for i in range(self.numLoci):
            #skip the same column
            if di[i] == 0:
                continue
            for j in range(i + 1, self.numLoci):
                # skip the same column
                if di[j] == 0:
                    continue
                alA = next(iter(allcnt[i]))
                alB = next(iter(allcnt[j]))

                hits = 0
                # the hit is wrong need to change
                for k in range(self.sampleSize):
                    # if alleleA[k][i] == alA and alleleA[k][j] == alB:
                    #     hits += 1
                    # if alleleB[k][i] == alA and alleleB[k][j] == alB:
                    #     hits += 1
                    if alleleA[k][j] == alB:
                        hits += 1
                    if alleleB[k][j] == alB:
                        hits += 1

                ai = float(allcnt[i][alA]) / float(totspots)
                bj = float(allcnt[j][alB]) / float(totspots)
                # NEVER GOES IN HITS???? -- did I fix it correctly? Email King
                # x = (float((float(hits) / float(len(data) - ai * bj)))) / (
                #     float(((ai * float(1 - ai) + float(di[i])) * (float(bj) * float(1 - bj) + di[j]))))
                # print(hits, alA, alB)
                if ai * (1 - ai) + di[i] == 0 or bj * (1 - bj) + di[j] == 0:
                    print("denominator is 0")
                    continue
                # x = (float(hits) / float(allcnt[i][alA]) - ai * bj) / (
                #     (ai*(1-ai) + di[i]) * (bj*(1-bj) + di[j]))
                # print(ai,bj,i,j)
                x = (float(hits) / float(2 * self.sampleSize) - ai * bj) / (
                        (ai * (1 - ai) + di[i]) * (bj * (1 - bj) + di[j]))

                running_sum += math.sqrt(abs(x))

        numloci = len(alleleA[0])
        self.stat1 = 2 * running_sum / (numloci * (numloci - 1))

        if (self.DEBUG):
            print("printing for stat1 end   ---->", self.stat1)
        # print("-----------------------------------------------------------------")


    ######################################################################
    # stat2 First Moment of Multilocus Homozygosity                     ##
    ######################################################################
    def stat2(self):
        # taking average homozygosity for each indiv then adding that all up and dividing by number of indivudls, so basically avg homozygosity over all indiv

        homozygosityCnt = 0
        tempVarStat2 = 0

        for i in range(self.sampleSize):
            for j in range(self.numLoci + 1):
                j = self.data.iloc[i, j]
                if (j == '0101' or j == '0202' or j == '0303' or j == '0404'):
                    homozygosityCnt = homozygosityCnt + 1
            tempStat2 = float(homozygosityCnt) / float(self.numLoci)
            tempVarStat2 = tempVarStat2 + tempStat2
            tempStat2 = 0
            homozygosityCnt = 0

        self.stat2 = float(tempVarStat2) / float(self.sampleSize)

        if (self.DEBUG):
            print("(First moment of homozygosity) Stats2 is ", self.stat2)

    ######################################################################
    # stat3 Second Moment of Multilocus Homozygosity                    ##
    ######################################################################
    def stat3(self):
        homozygosityCnt = 0
        totalHomozygosityDiff = 0
        # Total count is same as above stats
        for i in range(self.sampleSize):
            for j in range(self.numLoci + 1):
                j = self.data.iloc[i, j]
                # AA CC TT GG
                if (j == '0101' or j == '0202' or j == '0303' or j == '0404'):
                    homozygosityCnt = homozygosityCnt + 1
            if (homozygosityCnt > 0):
                homozygosityCnt = float(homozygosityCnt) / float(self.numLoci)
                difference = homozygosityCnt - self.stat2
                homozygosityCnt = 0
                totalHomozygosityDiff = totalHomozygosityDiff + (difference * difference)

        self.stat3 = float(totalHomozygosityDiff) / float(self.sampleSize - 1)

        if (self.DEBUG):
            print("(Second moment of multilocus homozygosity) Stats3 is ", self.stat3)

    ######################################################################
    # stat4 Wrights                                                     ##
    ######################################################################
    def stat4(self):
        tempstat4 = 0.0
        expected = 0
        num = []
        for i in range(self.numLoci):
            num *= 0
            homozygosityCnt = 0
            a = 0
            c = 0
            t = 0
            g = 0
            for j in range(self.sampleSize):
                # Checking freq of first two numbers by column
                if (self.data[j][i][2:] == '01'):
                    a = a + 1
                elif (self.data[j][i][2:] == '02'):
                    c = c + 1
                elif (self.data[j][i][2:] == '03'):
                    t = t + 1
                elif (self.data[j][i][2:] == '04'):
                    g = g + 1

                # Checking last two numbers
                if (self.data[j][i][:2] == '01'):
                    a = a + 1
                elif (self.data[j][i][:2] == '02'):
                    c = c + 1
                elif (self.data[j][i][:2] == '03'):
                    t = t + 1
                elif (self.data[j][i][:2] == '04'):
                    g = g + 1

                if (self.data[j][i] == '0101' or self.data[j][i] == '0202' or self.data[j][i] == '0303' or self.data[j][
                    i] == '0404'):
                    homozygosityCnt = homozygosityCnt + 1

            homozygosityCnt = float(homozygosityCnt) / float(self.sampleSize)
            homozygosityCnt = homozygosityCnt * homozygosityCnt
            temp = 1 - homozygosityCnt
            expected = expected + temp

            divisor = self.numLoci * 2
            a = a / float(divisor)
            c = c / float(divisor)
            t = t / float(divisor)
            g = g / float(divisor)

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
                    addStat4 = float(i) / float(expected)

                    tempstat4 = tempstat4 + addStat4

        tempstat4 = float(tempstat4) / float(self.numLoci)
        self.stat4 = 1 - tempstat4

        if (self.DEBUG):
            print('(Wrights) Stat4 is ', self.stat4)

    ######################################################################
    # stat5 Expected Heterozygosity                                     ##
    ######################################################################
    def stat5(self):
        tempstat5 = 0
        for i in range(self.numLoci + 1):
            homozygosityCnt = 0
            for j in range(self.sampleSize):
                j = self.data.iloc[j, i]
                if (j == '0101' or j == '0202' or j == '0303' or j == '0404'):
                    homozygosityCnt = homozygosityCnt + 1
            if (homozygosityCnt > 0):
                homozygosityCnt = float(homozygosityCnt) / float(self.sampleSize)
                totalhomozygosityCnt = homozygosityCnt * homozygosityCnt
                temp = 1 - totalhomozygosityCnt
                tempstat5 = tempstat5 + temp  # New heterozygosity value

        tempstat5 = float(tempstat5) / float(self.numLoci)
        self.stat5 = tempstat5

        if (self.DEBUG):
            print("(Expected heterozygosity) stat5 is ", self.stat5)

    ######################################################################
    # stat4 Updated after meeting w Dav                                 ##
    ######################################################################
    def newStat4(self):
        # https://academic.oup.com/jhered/article/106/3/306/2961865
        expected = 0
        homozygosityCnt = 0
        newstat4 = 0
        allcnt = self.allcnt
        totalNum = self.sampleSize * 2
        #Not sure why need to add 1? email isha to ask about this
        for i in range(self.numLoci + 1):
            for j in range(self.sampleSize):
                j = self.data.iloc[j, i]
                # print(j)
                if (j == '0101' or j == '0202' or j == '0303' or j == '0404'):
                    homozygosityCnt = homozygosityCnt + 1

            if (homozygosityCnt > 0):
                homozygosityCnt = float(homozygosityCnt) / float(self.sampleSize)
                tempObs = 1 - float(homozygosityCnt)
                # print(allcnt[i-1].values())
                valA = list(allcnt[i-1].values())[0]/totalNum
                if valA == 1:
                    continue
                valB = list(allcnt[i-1].values())[1]/totalNum
                # print("valA, valB", valA, valB)

                expected = 1 - (valA * valB)
                newstat4 = newstat4 + float(tempObs / expected)
                homozygosityCnt = 0

        newstat4 = 1 - (float(newstat4 / self.numLoci))
        self.stat4 = newstat4
        if (self.DEBUG):
            print("New stat4:   ", newstat4)
