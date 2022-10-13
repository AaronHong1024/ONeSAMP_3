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
    homoLoci = []

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
        # print("currLociCnt: ", currLociCnt)
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
        # print("NE value: ", NE_VALUEtemp)
        # print("-----------------------------------------------------------------")

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
        # print(data)
        #
        # print(test)
        NE_VALUEtemp = 0
        if len(test) != data.shape[1] + 1:
            NE_VALUEtemp = test[0]

        self.numLoci = data.shape[1] - 1

        self.NE_VALUE = NE_VALUEtemp
        self.sampleSize = data.shape[0]
        self.data = data
        # print(data)
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
    # stat1 BW Estimator                                                ##
    ######################################################################

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
        # print("alleleA: ", alleleA)
        # print("alleleB: ", alleleB)
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

        # print("homoloci: ", homoloci)
        # print("allcnt: ", allcnt)
        self.allcnt = allcnt
        di = []  # a 1D array that holds the departures of each loci from Hardy-Weinberg equilibrium
        totspots = float(
            self.sampleSize * 2)  # total number of alleles per locus. used in computing allele freq per locus
        # print("totspots: ", totspots)
        for i in range(len(allcnt)):  # fills up di
            vals = []
            for key, value in allcnt[i].items():
                vals.append(float(value) / float(totspots))
            # print("vals: ", vals)
            di.append(homoloci[i] - math.pow(vals[0], 2))

        # print("di: ", di)
        hits = 0

        running_sum = 0

        for i in range(self.numLoci):
            # skip the same column
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

        tempVarStat2 = 0
        data = self.data
        # store the homo number into a list (use more space but save more time)
        homoLoci = [0 for _ in range(self.numLoci)]
        i = 0

        for row in data.itertuples():
            homozygosityCnt = 0
            homozygosityCnt += row.count('0101')
            homozygosityCnt += row.count('0202')
            homozygosityCnt += row.count('0303')
            homozygosityCnt += row.count('0404')
            homoLoci[i] = homozygosityCnt
            i += 1
            tempVarStat2 += float(homozygosityCnt) / float(self.numLoci)

        self.homoLoci = homoLoci
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
        data = self.data
        homoLoci = self.homoLoci
        for i in range(self.sampleSize):
            homozygosityCnt = homoLoci[i]
            if homozygosityCnt > 0:
                homozygosityCnt = float(homozygosityCnt) / float(self.numLoci)
                difference = homozygosityCnt - self.stat2
                totalHomozygosityDiff = totalHomozygosityDiff + np.power(difference, 2)
        self.stat3 = float(totalHomozygosityDiff) / float(self.sampleSize - 1)

        if (self.DEBUG):
            print("(Second moment of multilocus homozygosity) Stats3 is ", self.stat3)

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
        # Not sure why need to add 1? email isha to ask about this
        for i in range(self.numLoci + 1):
            for j in range(self.sampleSize):
                j = self.data.iloc[j, i]
                if (j == '0101' or j == '0202' or j == '0303' or j == '0404'):
                    homozygosityCnt = homozygosityCnt + 1
            if (homozygosityCnt > 0):
                homozygosityCnt = float(homozygosityCnt) / float(self.sampleSize)
                observed = homozygosityCnt
                # print(allcnt[i-1].values())
                valA = list(allcnt[i - 1].values())[0] / totalNum
                if valA == 1:
                    continue
                valB = list(allcnt[i - 1].values())[1] / totalNum
                # print("valA, valB", valA, valB)

                expected = 1 - math.pow(valA, 2) - math.pow(valB, 2)
                newstat4 = newstat4 + float(observed / expected)
                homozygosityCnt = 0

        newstat4 = 1 - (float(newstat4 / self.numLoci))
        self.stat4 = newstat4
        if (self.DEBUG):
            print("New stat4:   ", newstat4)

        ######################################################################
        # stat5 Expected Heterozygosity                                     ##
        ######################################################################
        # slightly different form paper. Need to be reviewed.


    def stat5(self):
        tempstat5 = 0
        data =self.data
        allCnt = self.allcnt
        Individual_Number = self.sampleSize * 2
        FreqCount = 0
        for i in range(len(allCnt)):

            #Get the first allele's number, then divide the whole individual number
            FreqH = list(allCnt[i].values())[0] / Individual_Number
            FreqCount += FreqH
            # print(FreqH)
        # tempstat5 = FreqCount / self.numLoci
        self.stat5 = FreqCount / self.numLoci

        if (self.DEBUG):
            print("(Expected heterozygosity) stat5 is ", self.stat5)

