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
        with open(myFileName, 'r') as f:
            lines = f.readlines()
        result = []
        for line in lines[1:]:
            if (len(line.split()) > 10):
                result.append([m[i] for i in line.split()[2:]])
        data = np.asarray(result)
        self.data = data
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
            if numMissing > self.numLoci * PERCENT_MISSINGIndiv:
                deleteRow.append(i)
        if len(deleteRow) != 0:
            newData = np.delete(data, deleteRow, axis = 0)
            self.data = newData
            self.sampleSize = self.sampleSize - len(deleteRow)
            print("filter individuals")



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

            if numMissing > self.sampleSize * PERCENT_MISSINGLoci:
                deleteCol.append(i)

        if len(deleteCol) != 0:
            newData = np.delete(data, deleteCol, axis = 1)
            self.data = newData
            self.numLoci = self.numLoci - len(deleteCol)
            print("filter loci")


    ######################################################################
    # stat1 BW Estimator                                                ##
    ######################################################################

    def test_stat1(self):
        if self.DEBUG:
            print("printing for stat1 begin: ")

        data = self.data
        deletecol = []
        for i in range(len(data[0])):
            temp = data[:, i, :]
            unique = np.unique(temp)
            if len(unique) == 1:
                deletecol.append(i)

        if len(deletecol) != 0:
            data = np.delete(data, deletecol, axis=1)
        numloci = data.shape[1]
        sampleSize = data.shape[0]
        # compute the frequency of each allele
        allcnt = []
        # compute the homoloci number for data matrix
        homolociArray = []
        di = []
        totalspots = sampleSize * 2
        running_sum = 0
        sampCorrection = 2 / (numloci * (numloci - 1))
        r = 0
        for i in range(numloci):
            temp = data[:, i, :]
            # Can be optimized
            homoloci = np.sum(np.logical_and(temp[:, 1] == temp[:, 0], temp[:, 1] == temp[0][0],
                                             temp[:, 0] == temp[0][0]) == True) / sampleSize
            homolociArray.append(homoloci)
            currCnt = np.sum(temp == temp[0][0])
            allcnt.append(currCnt)
            currDi = homoloci - ((currCnt / totalspots) ** 2)
            di.append(currDi)


        for i in range(numloci):
            if di[i] == 0:
                continue
            LociA = data[:,i,:]
            for j in range(i + 1, numloci):
                if di[j] == 0:
                    continue
                LociB = data[:, j, :]
                index_A = np.sum((LociA == LociA[0][0]).astype(int),axis=1)
                index_B = np.sum((LociB == LociB[0][0]).astype(int),axis=1)
                hits = np.sum(index_A*index_B)
                ai = float(allcnt[i] / totalspots)
                bj = float(allcnt[j] / totalspots)
                if ai * (1 - ai) + di[i] == 0 or bj * (1 - bj) + di[j] == 0:
                    # print("denominator is 0")
                    continue
                    # x = (float(hits) / float(allcnt[i][alA]) - ai * bj) / (
                    #     (ai*(1-ai) + di[i]) * (bj*(1-bj) + di[j]))
                    # print(ai,bj,i,j)
                jointAB = float(hits / (2 * totalspots))

                denominator = (ai * (1 - ai) + di[i]) * (bj * (1 - bj) + di[j])

                r_intermdediate = 4*float((jointAB - ai*bj) ** 2) / denominator

                r += r_intermdediate


        running_sum = r

        # self.allcnt = allcnt

        self.stat1 = running_sum * sampCorrection


        if (self.DEBUG):
            print("printing for teststat1 end   ---->", self.stat1)



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
        totalNum = self.sampleSize * 2
        newstat4 = 0
        allcnt = []
        for i in range(self.numLoci):
            temp = data[:,i, :]
            # Can be optimized
            homoloci = np.sum(temp[:, 1] == temp[:, 0]) / self.sampleSize
            currCnt = np.sum(temp == temp[0][0])
            allcnt.append(currCnt)
            if homoloci > 0:
                observed = 1 - homoloci
                valA = currCnt / totalNum
                if valA == 1:
                    continue
                valB = float(1 - valA)
                expected = 1 - math.pow(valA, 2) - math.pow(valB, 2)
                newstat4 = newstat4 + float(observed / expected)


        newstat4 = 1 - (float(newstat4 / self.numLoci))
        self.stat4 = newstat4
        self.allcnt = allcnt
        if (self.DEBUG):
            print("New stat4:   ", newstat4)

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


