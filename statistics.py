#!/usr/bin/python
import sys

import math
import os
import re

import numpy as np

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

elements = np.array([[0, 0], [1, 1], [2, 2], [3, 3]])


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
    homoLociCol = []
    hexp = []
    hob = []

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
            NE_VALUEtemp = int(float(test[0]))
        self.NE_VALUE = NE_VALUEtemp
        self.numLoci = self.data.shape[1]
        if self.numLoci > 5000:
            print("error: loci size should be smaller than 5000")
            sys.exit(1)
        self.sampleSize = self.data.shape[0]


    def filterMonomorphicLoci(self):
        data = self.data
        deleteCol = []
        for i in range(self.numLoci):
            temp = data[:,i,:]
            if len(np.unique(temp)) == 1:
                deleteCol.append(i)
        print(deleteCol)
        if len(deleteCol) != 0:
            newData = np.delete(data, deleteCol, axis = 1)
            self.data = newData
            self.numLoci = self.numLoci - len(deleteCol)
            print("filter monomorphic loci")


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


    ######################################################################
    # stat1 BW Estimator                                                ##
    ######################################################################
#TODO: when we have small individuals the r result is slightly different from the LDNE
    def test_stat1(self):
        if self.DEBUG:
            print("printing for stat1 begin: ")

        data = self.data

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
        index = 0
        deletecol=[]
        # can be optimized. Merge it into the next for loop
        for i in range(numloci):
            temp = data[:, i, :]
            # Can be optimized
            # test = (temp == [0,0] + temp == [1,1] + temp == [2,2] + temp == [3,3]).any()

            homoloci = np.sum(np.logical_and(temp[:, 1] == temp[:, 0], temp[:, 1] == temp[0][0],
                                             temp[:, 0] == temp[0][0]) == True) / sampleSize

            allHomoloci = np.sum(np.count_nonzero((temp[:, np.newaxis, :] == elements).all(axis=2), axis=0)) / sampleSize
            homolociArray.append(allHomoloci)
            currCnt = np.sum(temp == temp[0][0])
            if currCnt == totalspots:
                deletecol.append(i)
            allcnt.append(currCnt)
            currDi = homoloci - ((currCnt / totalspots) ** 2)
            di.append(currDi)


# delete the poly column and recalculate the numloci and sample size
        if len(deletecol) != 0:
            data = np.delete(data, deletecol, axis=1)
            numloci = data.shape[1]
            sampleSize = data.shape[0]
            new_di = [di[i] for i in range(len(di)) if i not in deletecol]
            di = new_di
            sampCorrection = 2 / (numloci * (numloci - 1))
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
                currCntA = np.sum(LociA == LociA[0][0])
                currCntB = np.sum(LociB == LociB[0][0])
                ai = float(currCntA / totalspots)
                bj = float(currCntB / totalspots)

                if ai * (1 - ai) + di[i] == 0 or bj * (1 - bj) + di[j] == 0:
                    continue
                jointAB = float(hits / (2 * totalspots))

                denominator = (ai * (1 - ai) + di[i]) * (bj * (1 - bj) + di[j])

                r_intermdediate = float((jointAB - ai*bj) ** 2) / denominator

                r += r_intermdediate


        running_sum = r
        self.allcnt = allcnt
        self.stat1 = running_sum * sampCorrection
        self.homoLociCol = homolociArray


        if (self.DEBUG):
            print("printing for teststat1 end   ---->", self.stat1)


    ######################################################################
    # stat2 First Moment of Multilocus Homozygosity                     ##
    ######################################################################

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
    # stat3 Second Moment of Multilocus Homozygosity                    ##
    ######################################################################

    def test_stat3(self):
        homolociRow = self.homoLoci
        self.stat3 = np.var(homolociRow, ddof=1)

        if self.DEBUG:
            print("(Second moment of multilocus homozygosity) Stats3 is ", self.stat3)


    ######################################################################
    # stat4 Updated after meeting w Dav                                 ##
    ######################################################################

    def test_stat4(self):

        # According to the gene diversities equation,
        # if the observe is 0, then the expected will also be 0

        data = self.data
        totalNum = self.sampleSize * 2
        hexp = np.asarray(self.hexp)
        heob = np.asarray(self.hob)
        # hexp = hexp[hexp != 0]
        # heob = heob[heob != 0]
        # fis = 1 - heob/hexp

        fis = np.zeros_like(hexp)
        mask = hexp != 0
        fis[mask] = 1 - heob[mask] / hexp[mask]
        fis[np.logical_and(hexp == 0, heob == 0)] = 0
        if not mask[0] and heob[0] == 0:
            fis[0] = 0

        self.stat4 = np.sum(fis) / self.numLoci


        if (self.DEBUG):
            print("New stat4:   ", newstat4)

        ######################################################################
        # stat5 Expected Heterozygosity                                     ##
        ######################################################################
        # slightly different form paper. Need to be reviewed.

    #Currently, we have the missing data not filted.
    # thus, slightly different from the David example

    def test_stat5(self):
        sampleCorrection = self.sampleSize / (self.sampleSize - 1)
        totalNum = self.sampleSize * 2
        homoloci = np.asarray(self.homoLociCol)
        data = self.data
        allcnt = np.asarray(self.allcnt)
        freqA = (allcnt / totalNum)
        freqB = 1 - freqA
        h_obser = (1-homoloci)/totalNum
        freqRes = (1 - freqA ** 2 - freqB ** 2 - h_obser)*sampleCorrection
        self.hexp = freqRes
        self.hob = h_obser*totalNum
        self.stat5 = np.sum(freqRes) / self.numLoci

        if (self.DEBUG):
            print("(Expected heterozygosity) stat5 is ", self.stat5)


