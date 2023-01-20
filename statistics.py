#!/usr/bin/python

import math
import os
import re

import numpy as np
import pandas as pd
import collections
from decimal import *

m = {'0101': [0, 0],
     '0102': [0, 1],
     '0103': [0, 2],
     '0104': [0, 3],
     '0201': [1, 0],
     '0202': [1, 1],
     '0203': [1, 2],
     '0204': [1, 3],
     '0301': [2, 0],
     '0302': [2, 1],
     '0303': [2, 2],
     '0304': [2, 3],
     '0401': [3, 0],
     '0402': [3, 1],
     '0403': [3, 2],
     '0404': [3, 3],
     '0000': [-1,-1],
     '0100': [0,-1],
     '0001': [-1,0],
     '0200': [1,-1],
     '0300': [2,-1],
     '0400': [3,-1],
     '0002': [-1,1],
     '0003': [-1,2],
     '0004': [-1,3]}


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

    # m is the transfer array

    ######################################################################
    # writeStatistics                                                   ##
    ######################################################################
    # def writeStatistics(self, myFileName):
    # outputFile = open(myFileName, "w")

    ######################################################################
    # readData                                                          ##
    ######################################################################
    def readData(self, myFileName):
        print("test read filename: ", myFileName)
        with open(myFileName, 'r') as f:
            lines = f.readlines()
        result = []
        for line in lines[1:]:
            if (len(line.split()) > 10):
                result.append([m[i] for i in line.split()[2:]])
        self.data = np.asarray(result)
        last_line = self.get_file_last_line(os.path.abspath(myFileName))
        test = last_line.strip().split(" ")
        NE_VALUEtemp = 0
        if len(test) != self.data.shape[1] + 2:
            NE_VALUEtemp = test[0]
        self.NE_VALUE = NE_VALUEtemp
        self.numLoci = self.data.shape[1]
        self.sampleSize = self.data.shape[0]


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
        data = self.data
        deleteRow = []
        # to each row, if there exist 0*00 delete this row.
        for i in range(self.sampleSize):
            numMissing = 0
            temp = data[i,:,:]
            numMissing += np.sum(np.logical_or(temp[:,0] == -1, temp[:,1] == -1)) - np.sum((temp==[-1,-1]).all())
            if numMissing > 0:
                deleteRow.append(i)
        if len(deleteRow) != 0:
            newData = np.delete(data, deleteRow, axis = 0)
            self.data = newData
            self.sampleSize = self.sampleSize - 1
            print("filter individuals")



        # for i in range(self.sampleSize):
        #     individual = self.data.loc[i]
        #     numMissing = 0
        #     for j in (individual):
        #         if (j == '0100' or j == '0001' or j == '0000'):
        #             numMissing += 1
        #     if (numMissing > PERCENT_MISSINGIndiv * self.numLoci):
        #         print("Deleted:", self.data[j])
        #         self.ARRAY_MISSINGIndiv.append(self.data[j])
        #         del self.data[j]
        #         self.numLoci = self.numLoci - 1
        #         self.sampleSize = len(self.data) - 1

                # only that place., but double check with David

    ######################################################################
    # filterLoci                                                        ##
    ######################################################################
    def filterLoci(self, PERCENT_MISSINGLoci):
        data = self.data
        deleteCol = []

        for i in range(self.numLoci):
            temp = data[:,i,:]
            numMissing = 0
            numMissing += np.sum(np.logical_or(temp[:,0] == -1, temp[:,1] == -1)) - np.sum((temp==[-1,-1]).all())

            if numMissing > 0:
                deleteCol.append(i)

        if len(deleteCol) != 0:
            newData = np.delete(data, deleteCol, axis = 1)
            self.data = newData
            self.numLoci = self.numLoci - 1
            print("filter loci")
        #
        #
        # for i in range(self.numLoci):
        #     individual = self.data[i]
        #     numMissing = 0
        #     for j in (individual):
        #         if (j == '0100' or j == '0001' or j == '0000'):
        #             numMissing += 1
        #     if (numMissing > PERCENT_MISSINGLoci * self.numLoci):
        #         print("Deleted:", self.data[j])
        #         self.ARRAY_MISSINGLoci.append(self.data[j])
        #         del self.data[j]
        #         self.numLoci = self.numLoci - 1
        #         self.sampleSize = len(self.data) - 1

    ######################################################################
    # stat1 BW Estimator                                                ##
    ######################################################################

    def test_stat1(self):
        if self.DEBUG:
            print("printing for stat1 begin: ")

        data = self.data
        # compute the frequency of each allele
        allcnt = []
        # compute the homoloci number for data matrix
        homolociArray = []
        di = []
        totalspots = self.sampleSize * 2
        running_sum = 0
        numloci = self.numLoci
        for i in range(self.numLoci):
            temp = data[:, i, :]
            # Can be optimized
            homoloci = np.sum(np.logical_and(temp[:, 1] == temp[:, 0], temp[:, 1] == temp[0][0],
                                             temp[:, 0] == temp[0][0]) == True) / self.sampleSize
            homolociArray.append(homoloci)
            currCnt = np.sum(temp == temp[0][0])
            allcnt.append(currCnt)
            currDi = homoloci - ((currCnt / totalspots) ** 2)
            di.append(currDi)


        for i in range(self.numLoci):
            if di[i] == 0:
                continue
            LociA = data[:,i,:]
            for j in range(i + 1, self.numLoci):
                if di[j] == 0:
                    continue
                LociB = data[:, j, :]
                index_A = np.sum((LociA == LociA[0][0]).astype(int),axis=1)
                index_B = np.sum((LociB == LociB[0][0]).astype(int),axis=1)
                hits = np.sum(index_A*index_B)
                ai = allcnt[i] / totalspots
                bj = allcnt[j] / totalspots
                if ai * (1 - ai) + di[i] == 0 or bj * (1 - bj) + di[j] == 0:
                    # print("denominator is 0")
                    continue
                    # x = (float(hits) / float(allcnt[i][alA]) - ai * bj) / (
                    #     (ai*(1-ai) + di[i]) * (bj*(1-bj) + di[j]))
                    # print(ai,bj,i,j)
                x = (float(hits) / float(2*self.sampleSize) - ai * bj) / (
                        (ai * (1 - ai) + di[i]) * (bj * (1 - bj) + di[j]))
                running_sum += math.sqrt(abs(x))

        self.allcnt = allcnt
        self.stat1 = 2 * running_sum / (numloci * (numloci - 1))

        if (self.DEBUG):
            print("printing for teststat1 end   ---->", self.stat1)


    # def new_stat1(self):
    #     if self.DEBUG:
    #         print("printing for stat1 begin: ")
    #
    #     data = self.data
    #     alleleA = []
    #     alleleB = []
    #
    #     for row in data.itertuples(index=False):
    #         temA = []
    #         temB = []
    #         for col in range(1, len(row)):
    #             temA.append(row[col][0:2])
    #             temB.append(row[col][2:])
    #
    #         alleleA.append(temA)
    #         alleleB.append(temB)
    #
    #     allcnt = []
    #     homoloci = []
    #
    #     for i in range(len(alleleA[0])):
    #         newdic = {}
    #         temp = 0
    #         for j in range(len(alleleA)):
    #             if alleleA[j][i] == alleleB[j][i] and alleleA[j][i] == alleleA[0][i]:
    #                 temp += 1
    #             if alleleA[j][i] in newdic:
    #                 newdic[alleleA[j][i]] += 1
    #             else:
    #                 newdic[alleleA[j][i]] = 1
    #             if alleleB[j][i] in newdic:
    #                 newdic[alleleB[j][i]] += 1
    #             else:
    #                 newdic[alleleB[j][i]] = 1
    #
    #         homoloci.append(float(temp) / float(self.sampleSize))
    #         allcnt.append(newdic)
    #
    #     self.allcnt = allcnt
    #     di = []  # a 1D array that holds the departures of each loci from Hardy-Weinberg equilibrium
    #     totspots = float(
    #         self.sampleSize * 2)  # total number of alleles per locus. used in computing allele freq per locus
    #     # print("totspots: ", totspots)
    #     for i in range(len(allcnt)):  # fills up di
    #         vals = []
    #         for key, value in allcnt[i].items():
    #             vals.append(float(value) / float(totspots))
    #         di.append(homoloci[i] - math.pow(vals[0], 2))
    #
    #     hits = 0
    #
    #     running_sum = 0
    #
    #     for i in range(self.numLoci):
    #         # skip the same column
    #         if di[i] == 0:
    #             continue
    #         for j in range(i + 1, self.numLoci):
    #             # skip the same column
    #             if di[j] == 0:
    #                 continue
    #             alA = next(iter(allcnt[i]))
    #             alB = next(iter(allcnt[j]))
    #             hits = 0
    #             # the hit is wrong need to change
    #             for k in range(self.sampleSize):
    #                 # if alleleA[k][i] == alA and alleleA[k][j] == alB:
    #                 #     hits += 1
    #                 # if alleleB[k][i] == alA and alleleB[k][j] == alB:
    #                 #     hits += 1
    #                 if alleleA[k][j] == alB:
    #                     hits += 1
    #                 if alleleB[k][j] == alB:
    #                     hits += 1
    #             ai = float(allcnt[i][alA]) / float(totspots)
    #             bj = float(allcnt[j][alB]) / float(totspots)
    #             # NEVER GOES IN HITS???? -- did I fix it correctly? Email King
    #             # x = (float((float(hits) / float(len(data) - ai * bj)))) / (
    #             #     float(((ai * float(1 - ai) + float(di[i])) * (float(bj) * float(1 - bj) + di[j]))))
    #             # print(hits, alA, alB)
    #             if ai * (1 - ai) + di[i] == 0 or bj * (1 - bj) + di[j] == 0:
    #                 print("denominator is 0")
    #                 continue
    #             # x = (float(hits) / float(allcnt[i][alA]) - ai * bj) / (
    #             #     (ai*(1-ai) + di[i]) * (bj*(1-bj) + di[j]))
    #             # print(ai,bj,i,j)
    #             x = (float(hits) / float(2 * self.sampleSize) - ai * bj) / (
    #                     (ai * (1 - ai) + di[i]) * (bj * (1 - bj) + di[j]))
    #
    #             running_sum += math.sqrt(abs(x))
    #
    #     numloci = len(alleleA[0])
    #     self.stat1 = 2 * running_sum / (numloci * (numloci - 1))
    #
    #     if (self.DEBUG):
    #         print("printing for stat1 end   ---->", self.stat1)
    #     # print("-----------------------------------------------------------------")

    def test_stat2(self):
        data = self.data
        homolociRow = []
        for i in range(self.sampleSize):
            temp = data[i,:,:]
            homoloci = np.sum(temp[:, 1] == temp[:, 0])
            homolociRow.append(homoloci)

        self.homoLoci = homolociRow
        self.stat2 = np.mean(homolociRow)

        if (self.DEBUG):
            print("(First moment of homozygosity) test Stats2 is ", self.stat2)


    ######################################################################
    # stat2 First Moment of Multilocus Homozygosity                     ##
    ######################################################################
    # def stat2(self):
    #     # taking average homozygosity for each indiv then adding that all up and dividing by number of indivudls, so basically avg homozygosity over all indiv
    #
    #     tempVarStat2 = 0
    #     data = self.data
    #     # store the homo number into a list (use more space but save more time)
    #     homoLoci = [0 for _ in range(self.sampleSize)]
    #     i = 0
    #
    #     for row in data.itertuples():
    #         homozygosityCnt = 0
    #         homozygosityCnt += row.count('0101')
    #         homozygosityCnt += row.count('0202')
    #         homozygosityCnt += row.count('0303')
    #         homozygosityCnt += row.count('0404')
    #         homoLoci[i] = homozygosityCnt
    #         i += 1
    #         # tempVarStat2 += homozygosityCnt
    #     self.homoLoci = homoLoci
    #     self.stat2 = np.mean(homoLoci)
    #     # self.stat2 = float(tempVarStat2) / float(self.sampleSize)
    #
    #     if (self.DEBUG):
    #         print("(First moment of homozygosity) Stats2 is ", self.stat2)

    ######################################################################
    # stat3 Second Moment of Multilocus Homozygosity                    ##
    ######################################################################

    def test_stat3(self):
        homolociRow = self.homoLoci
        self.stat3 = np.var(homolociRow, ddof=1)

        if self.DEBUG:
            print("(Second moment of multilocus homozygosity) Stats3 is ", self.stat3)

    # def stat3(self):
    #     homoLoci = self.homoLoci
    #     ## Compute the sample variance
    #     self.stat3 = np.var(homoLoci, ddof=1)
    #
    #     if self.DEBUG:
    #         print("(Second moment of multilocus homozygosity) Stats3 is ", self.stat3)

    ######################################################################
    # stat4 Updated after meeting w Dav                                 ##
    ######################################################################

    def test_stat4(self):
        data = self.data
        allcnt = self.allcnt
        totalNum = self.sampleSize * 2
        newstat4 = 0
        for i in range(self.numLoci):
            temp = data[:,i, :]
            # Can be optimized
            homoloci = np.sum(temp[:, 1] == temp[:, 0]) / self.sampleSize

            if homoloci > 0:
                observed = homoloci
                valA = allcnt[i] / totalNum
                if valA == 1:
                    continue
                valB = 1 - valA
                expected = 1 - math.pow(valA, 2) - math.pow(valB, 2)
                newstat4 = newstat4 + float(observed / expected)


        newstat4 = 1 - (float(newstat4 / self.numLoci))
        self.stat4 = newstat4
        if (self.DEBUG):
            print("New stat4:   ", newstat4)


    # def newStat4(self):
    #     # https://academic.oup.com/jhered/article/106/3/306/2961865
    #     expected = 0
    #     homozygosityCnt = 0
    #     newstat4 = 0
    #     allcnt = self.allcnt
    #     totalNum = self.sampleSize * 2
    #     # Not sure why need to add 1? email isha to ask about this
    #     for i in range(self.numLoci + 1):
    #         for j in range(self.sampleSize):
    #             j = self.data.iloc[j, i]
    #             if (j == '0101' or j == '0202' or j == '0303' or j == '0404'):
    #                 homozygosityCnt = homozygosityCnt + 1
    #         if (homozygosityCnt > 0):
    #             homozygosityCnt = float(homozygosityCnt) / float(self.sampleSize)
    #             observed = homozygosityCnt
    #             # print(allcnt[i-1].values())
    #             valA = list(allcnt[i - 1].values())[0] / totalNum
    #             if valA == 1:
    #                 continue
    #             valB = list(allcnt[i - 1].values())[1] / totalNum
    #             expected = 1 - math.pow(valA, 2) - math.pow(valB, 2)
    #             newstat4 = newstat4 + float(observed / expected)
    #             homozygosityCnt = 0
    #
    #     newstat4 = 1 - (float(newstat4 / self.numLoci))
    #     self.stat4 = newstat4
    #     if (self.DEBUG):
    #         print("New stat4:   ", newstat4)

        ######################################################################
        # stat5 Expected Heterozygosity                                     ##
        ######################################################################
        # slightly different form paper. Need to be reviewed.

    def test_stat5(self):
        tempstat5 = 0
        data = self.data
        allcnt = np.asarray(self.allcnt)
        totalNum = self.sampleSize * 2
        freqA = (allcnt / totalNum)
        freqB = 1 - freqA
        freqRes = 1 - freqA ** 2 - freqB ** 2
        # stat5 = np.sum(freqRes) / self.numLoci
        self.stat5 = np.sum(freqRes) / self.numLoci

        if (self.DEBUG):
            print("(Expected heterozygosity) stat5 is ", self.stat5)



    # def stat5(self):
    #     tempstat5 = 0
    #     data = self.data
    #     allCnt = self.allcnt
    #     totalNum = self.sampleSize * 2
    #
    #     for i in range(self.numLoci):
    #         freqP = 0
    #         for value in allCnt[i].values():
    #             freqP = freqP + (value / totalNum) ** 2
    #         tempstat5 = tempstat5 + (1 - freqP)
    #     self.stat5 = tempstat5 / self.numLoci
    #
    #     if (self.DEBUG):
    #         print("(Expected heterozygosity) stat5 is ", self.stat5)
