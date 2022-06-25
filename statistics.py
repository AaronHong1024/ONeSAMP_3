#!/usr/bin/python

import math

class statisticsClass:

####### Members of Class

    ARRAY_MISSINGIndiv = []
    ARRAY_MISSINGLoci = []
    data = []      ## Matrix is [individuals by loci]
    stat1 = 0
    stat2 = 0
    stat3 = 0
    stat4 = 0
    stat5 = 0
    numLoci = 0
    sampleSize = 0  ##Indivduals
    NE_VALUE = 0
    DEBUG = 0


    ######################################################################
    # writeStatistics                                                   ##
    ######################################################################
    # def writeStatistics(self, myFileName):
    # outputFile = open(myFileName, "w")

    ######################################################################
    # readData                                                          ##
    ######################################################################
    def readData(self, myFileName):
        matrixFile = open(myFileName, "r")
        line = matrixFile.readline()
        sampleSize = 0

        # Read until "Pop" in file
        popReached = 0;
        while line:
            if ((line == "Pop\n") or (line == "POP\n") or (line == "pop\n")) :
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
        temp = []
        temp = line.split()
        line = matrixFile.readline()
        temp = line.split()
        temp.pop(1)
        currLociCnt = len(temp) #Getting length of first line of data


        for i in range(sampleSize):
            temp = line.split()
            temp.pop(1)  # Getting rid of comma in array
            data.append(temp)
            line = matrixFile.readline()
            if (len(temp) != currLociCnt):
                print(len(temp))
                print("ERROR:statistics.py:line 56:: There is an incorrect number of loci. Fatal Error.")
                exit()

        #Reading in NE Value
        NE_VALUEtemp = 0
        if(line):
            temp = line.split()
            NE_VALUEtemp = temp[0]


        self.data = data
        self.numLoci = currLociCnt - 1 #Subracting 1 bc number takes into acct the indiviudal's number
        self.sampleSize = len(self.data)
        self.NE_VALUE = NE_VALUEtemp

    ######################################################################
    # filterIndividuals                                                 ##
    # Filters the data for all individuals that have > 20% missing data ##
    ######################################################################
    def filterIndividuals(self, PERCENT_MISSINGIndiv):
        for i in range(self.sampleSize):
            individual = self.data[i]
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

    ######################################################################
    # filterLoci                                                        ##
    ######################################################################
    def filterLoci(self, PERCENT_MISSINGLoci):
        for i in range(len(self.data)):
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
        homoloci = []  # maintains frequency counts of homologous alleles per loci
        denominator = float(numLoci)

        di = []  # a 1D array that holds the departures of each loci from Hardy-Weinberg equilibrium
        totspots = float(len(data) * 2)  # total number of alleles per locus. used in computing allele freq per locus

        for i in sampleSize:  # fills up allcnt
            newdic = [0,0,0,0]
            temp = 0
            for j in numLoci:

                # check if it is homozygous
                if (alleleA[i][j] == alleleB[i][j]):
                    temp += 1

                newdic[int(alleleA[i][j])] +=1
                newdic[int(alleleB[i][j])] +=1

                #if (alleleA[i][j] in newdic):
                 #   newdic[alleleA[i][j]] += 1
                #else:
                 #   newdic[alleleA[i][j]] = 1

                #if (alleleB[i][j] in newdic):
                 #   newdic[alleleB[i][j]] += 1
                #else:
                 #   newdic[alleleB[i][j]] = 1
                # END INNER FOR LOOP

            currHomoLoci = (temp / denominator)
            currDeparture = ( numLoci - temp ) / totspots
            homoloci.append( currHomoLoci )
            allcnt.append(newdic)

            # vals = []
            # for key, value in allcnt[i].items():
            #    vals.append(float(value) / float(totspots))
            di.append( (currHomoLoci  *  currHomoLoci) - (currDeparture  *  currDeparture ))
        # END OUTER FOR LOOP

        hits = 0

        running_sum = 0

        for i in range(len(alleleA)):
            for j in range(len(alleleB)):
                keysi = []
                keysj = []

                for key, value in allcnt[i].items():
                    keysi.append(key)


                for key, value in allcnt[j].items():
                    keysj.append(key)


                if (len(keysj) < 2):
                    continue

                alA = keysi[0]
                alB = keysj[1]

                hits = 0
                for k in numLoci :
                    if ((alleleA[i][k] == alA or alleleB[i][k] == alA) and
                        (alleleA[j][k] == alB or alleleB[j][k] == alB )):
                        hits += 1
                ai = allcnt[i][alA] / totspots
                bj = allcnt[j][alB] / totspots

                x = (hits / (sampleSize - ai *  bj)) / (((ai * (1 - ai) + di[i]) * (float(bj) * float(1 - bj) + di[j])))
                running_sum += math.sqrt(abs(x))

        self.stat1 = 2*running_sum/(numloci*(numloci-1))

        if (self.DEBUG):
            print("printing for stat1 end   ---->", self.stat1)

    ######################################################################
    # stat1 BW Estimator                                                ##
    ######################################################################
    def stat1(self):
        if (self.DEBUG):
            print("printing for stat1 begin:")

        data = self.data
        alleleA = []  # same dimensions as data but holds the first allele in a locus
        alleleB = []  # same dimensions as data but holds the second allele in a locus

        #print(data)

        for i in range(len(data)):  # fills up alleleA and alleleB
            temA = []
            temB = []
            for j in range(1, len(data[i])):
                temA.append(data[i][j][:2])
                temB.append(data[i][j][2:])

            alleleA.append(temA)
            alleleB.append(temB)



        allcnt = []  # a 1D array, each element is a dictionary of alleles with frequency counts per loci
        homoloci = []  # maintains frequency counts of homologous alleles per loci
        for i in range(len(alleleA)):  # fills up allcnt
            newdic = {}
            temp = 0
            for j in range(len(alleleA[0])):
                if (alleleA[i][j] == alleleB[i][j] and alleleA[i][j] == alleleA[i][0]):
                    temp += 1

                if (alleleA[i][j] in newdic):
                    newdic[alleleA[i][j]] += 1
                else:
                    newdic[alleleA[i][j]] = 1

                if (alleleB[i][j] in newdic):
                    newdic[alleleB[i][j]] += 1
                else:
                    newdic[alleleB[i][j]] = 1

            homoloci.append(float(temp) / float(len(data)))
            allcnt.append(newdic)

        di = []  # a 1D array that holds the departures of each loci from Hardy-Weinberg equilibrium
        totspots = float(len(data) * 2)  # total number of alleles per locus. used in computing allele freq per locus

        for i in range(len(allcnt)):  # fills up di
            vals = []
            for key, value in allcnt[i].items():
                vals.append(float(value) / float(totspots))
            di.append(float(homoloci[i] * homoloci[i] - float(vals[0] * vals[0])))

        hits = 0

        running_sum = 0

        for i in range(len(alleleA)):
            for j in range(len(alleleB)):
                keysi = []
                keysj = []

                for key, value in allcnt[i].items():
                    keysi.append(key)


                for key, value in allcnt[j].items():
                    keysj.append(key)


                if (len(keysj) < 2):
                    continue

                alA = keysi[0]
                alB = keysj[1]

                hits = 0
                for k in range(len(alleleA[0])):
                    if ((alleleA[i][k] == keysi[0] or alleleB[i][k] == keysi[0]) and (
                            alleleA[j][k] == keysj[1] or alleleB[j][k] == keysj[1])):
                        hits += 1
                ai = float(allcnt[i][alA]) / float(totspots)
                bj = float(allcnt[j][alB]) / float(totspots)


                    #NEVER GOES IN HITS???? -- did I fix it correctly? Email King
                x = (float((float(hits) / float(len(data) - ai *  bj))))/(float(((ai * float(1 - ai) + float(di[i])) * (float(bj) * float(1 - bj) + di[j]))))
                running_sum += math.sqrt(abs(x))

        numloci = len(alleleA[0])
        self.stat1 = 2*running_sum/(numloci*(numloci-1))

        if (self.DEBUG):
            print("printing for stat1 end   ---->", self.stat1)

    ######################################################################
    # stat2 First Moment of Multilocus Homozygosity                     ##
    ######################################################################
    def stat2(self):
        #taking average homozygosity for each indiv then adding that all up and dividing by number of indivudls, so basically avg homozygosity over all indiv

        homozygosityCnt = 0
        tempVarStat2 = 0

        for i in range(self.sampleSize):
            for j in range(self.numLoci + 1):
                j = self.data[i][j]
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
                j = self.data[i][j]
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

                if (self.data[j][i] == '0101' or self.data[j][i] == '0202' or self.data[j][i] == '0303' or self.data[j][i] == '0404'):
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
            print('(Wrights) Stat4 is ',  self.stat4)

    ######################################################################
    # stat5 Expected Heterozygosity                                     ##
    ######################################################################
    def stat5(self):
        tempstat5 = 0
        for i in range(self.numLoci + 1):
            homozygosityCnt = 0
            for j in range(self.sampleSize):
                j = self.data[j][i]
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
        #https://academic.oup.com/jhered/article/106/3/306/2961865
        expected = 0
        homozygosityCnt = 0
        newstat4 = 0

        for i in range(self.numLoci + 1):
            for j in range(self.sampleSize):
                j = self.data[j][i]
                if (j == '0101' or j == '0202' or j == '0303' or j == '0404'):
                    homozygosityCnt = homozygosityCnt + 1


            if (homozygosityCnt > 0):
                homozygosityCnt = float(homozygosityCnt) / float(self.sampleSize)
                tempObs = 1 - float(homozygosityCnt)
                expected = 1 - (homozygosityCnt * homozygosityCnt)
                newstat4 = newstat4 + float(tempObs / expected)
                homozygosityCnt = 0



        newstat4 = 1 - (float(newstat4/self.numLoci))


        if (self.DEBUG):
            print ("New stat4:   " ,newstat4)
