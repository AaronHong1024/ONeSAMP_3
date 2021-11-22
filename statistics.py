import math

class statisticsClass:

####### Members of Class

    PERCENT_MISSING = 0.2
    data = []
#    numCol = 0      ##  = numLoci
 #   numRow = -1     ##  = sampleSize        ## Matrix is [individuals by loci]
    stat1 = 0
    stat2 = 0
    stat3 = 0
    stat4 = 0
    stat5 = 0
    numLoci = 0
    totalAlleles = 0
    sampleSize = 0;


    ######################################################################
    # readData                                                          ##
    ######################################################################
    def readData(self, myFileName):
        matrixFile = open(myFileName, "r")
        line = matrixFile.readline()

        # Read until "Pop" in file
        popReached = 0;
        while line:
            if ((line == "Pop\n") or (line == "POP\n") or (line == "pop\n")) :
                popReached = 1
                break
            line = matrixFile.readline()

        if (popReached == 0):
            print("ERROR:statistics.py:line 32:: POP not contained in file. Fatal error")
            exit()


        # Starts creating data matrix
        data = self.data
        temp = []
        temp = line.split()
        temp.pop(1)  # Getting rid of comma in array
        data.append(temp)
        line = matrixFile.readline()
        currLociCnt = tempLociCnt;

        while(line) :
            temp = line.split()
            temp.pop(1)  # Getting rid of comma in array
            data.append(temp)
            line = matrixFile.readline()
            if (len(temp) != currLociCnt):
                print("ERROR:statistics.py:line 56:: There is an incorrect number of loci. Fatal Error.")
                exit()

        self.data = data
        self.numLoci = currLociCnt
        self.sampleSize = len(self.data)
    #    self.numRow = self.sampleSize
     #   self.numCol = self.numLoci

    ######################################################################
    # filterIndividuals                                                 ##
    # Filters the data for all individuals that have > 20% missing data ##
    ######################################################################
    def filterIndividuals(self):
        for i in range(self.sampleSize):
            individual = self.data[i]
            numMissing = 0
            for j in (individual):
                if (j == '0100' or j == '0001' or j == '0000'):
                    numMissing += 1
            if (numMissing > PERCENT_MISSING*self.numLoci):
                del self.data[j]

    ######################################################################
    # filterLoci                                                        ##
    ######################################################################
    def filterLoci(self):
        for i in range(len(self.data)):
            individual = self.data[i]
            numMissing = 0
            for j in (individual):
                if (j == '0100' or j == '0001' or j == '0000'):
                    numMissing += 1
            if (numMissing > PERCENT_MISSING*self.numLoci):
                del self.data[j]

    ######################################################################
    # stat1                                                       ##
    ######################################################################
    def stat1(self):
        print("printing for stat1 begin:")
        data = self.data
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
        for j in range(len(alleleA[0])):  # fills up allcnt
            newdic = {}
            temp = 0
            for i in range(len(alleleA)):
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

            homoloci.append(temp / len(data))
            allcnt.append(newdic)

        print(allcnt[0])
        di = []  # a 1D array that holds the departures of each loci from Hardy-Weinberg equilibrium
        totspots = len(data) * 2  # total number of alleles per locus. used in computing allele freq per locus

        for i in range(len(allcnt)):  # fills up di
            vals = []
            for key, value in allcnt[i].items():
                vals.append(value / totspots)
            di.append(homoloci[i] * homoloci[i] - vals[0] * vals[0])

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
                    if ((alleleA[i][k] == keysi or alleleB[i][k] == keysi) and (
                            alleleA[j][k] == keysj or alleleB[j][k] == keysj)):
                        hits += 1
                ai = allcnt[i][alA] / totspots
                bj = allcnt[j][alB] / totspots
                x = (hits / len(data) - ai * bj) / ((ai * (1 - ai) + di[i]) * (bj * (1 - bj) + di[j]))
                running_sum += math.sqrt(abs(x))

        numloci = len(alleleA[0])
        print(2 * running_sum / (numloci * (numloci - 1)))
        # return 2*running_sum/(numloci*(numloci-1))
        print("printing for stat1 end")

    ######################################################################
    # stat2                                                       ##
    ######################################################################
    def stat2(self):
        #data = self.data
    #    numCol = self.numCol
    #    numRow = self.numRow
    #    numLoci = self.numLoci

        # Stat 2
        homoCount = 0
        totalAlleles = (numCol - 1) * numRow
        for i in range(numRow):
            for j in (self.data[i]):
                if (j == '0101' or j == '0202' or j == '0303' or j == '0404'):  # What was actg?
                    homoCount = homoCount + 1

        stat2 = homoCount / numRow
        print("(First moment of homozygosity) Stats2 is ", stat2)

        self.stat2 = stat2
        self.totalAlleles = totalAlleles

    ######################################################################
    # stat3                                                        ##
    ######################################################################
    def stat3(self):
        data = self.data
        numCol = self.numCol
        numRow = self.numRow
        numLoci = self.numLoci

        homoDiff = 0
        homoCount = 0
        totalHomoDiff = 0
        # Total count is same as above stats
        for i in range(numRow):
            for j in (data[i]):
                # AA CC TT GG
                if (j == '0101' or j == '0202' or j == '0303' or j == '0404'):  # What was actg?
                    homoCount = homoCount + 1
            homoDiff = homoCount - self.stat2
            homoCount = 0
            totalHomoDiff = totalHomoDiff + (homoDiff * homoDiff)

        stat3 = (totalHomoDiff) / (self.numRow - 1)
        print("(Second moment of multilocus homozygosity) Stats3 is ", stat3)
        self.stat3 = stat3

    def stat4(self):
        data = self.data
        numCol = self.numCol
        numRow = self.numRow
        numLoci = self.numLoci

        # Stat4
        stat4 = 0.0
        tempstat4 = 0.0
        num = []
        for i in range(numCol):
            # Setting to zero
            num *= 0
            addStat4 = 0.0
            homoCount = 0
            stat5 = 0
            a = 0
            c = 0
            t = 0
            g = 0
            for j in range(numRow):
                # Checking freq of first two numbers
                if (data[j][i][2:] == '01'):
                    a = a + 1
                elif (data[j][i][2:] == '02'):
                    c = c + 1
                elif (data[j][i][2:] == '03'):
                    t = t + 1
                elif (data[j][i][2:] == '04'):
                    g = g + 1

                # Checking last two numbers
                if (data[j][i][:2] == '01'):
                    a = a + 1
                elif (data[j][i][:2] == '02'):
                    c = c + 1
                elif (data[j][i][:2] == '03'):
                    t = t + 1
                elif (data[j][i][:2] == '04'):
                    g = g + 1

                if (data[j][i] == '0101' or data[j][i] == '0202' or data[j][i] == '0303' or data[j][i] == '0404'):
                    homoCount = homoCount + 1

            homoCount = homoCount / (numRow)
            homoCount = homoCount * homoCount
            temp = 1 - homoCount
            stat5 = stat5 + temp

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
                    addStat4 = float(i) / stat5

            stat4 = stat4 + addStat4


        stat4 = stat4 / numLoci
        stat4 = 1 - stat4
        print('(Wrights) Stat4 is ',  stat4)

        self.stat4 = stat4

    def stat5(self):
        data = self.data
        numCol = self.numCol
        numRow = self.numRow
        numLoci = self.numLoci

        a = 0
        c = 0
        t = 0
        g = 0
        stat5 = 0
        num = []
        for i in range(numCol):
            num *= 0
            addStat5 = 0
            a = 0
            c = 0
            t = 0
            g = 0
            for j in range(numRow):
                # Checking freq of first two numbers
                if (data[j][i][2:] == '01'):
                    a = a + 1
                elif (data[j][i][2:] == '02'):
                    c = c + 1
                elif (data[j][i][2:] == '03'):
                    t = t + 1
                elif (data[j][i][2:] == '04'):
                    g = g + 1

                # Checking last two numbers
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

        print('stat5 before dividing', stat5)
        stat5 = stat5 / numLoci
        print('Stat5 is ', stat5)

        self.stat5 = stat5

    def stat5(self):
        data = self.data
        numCol = self.numCol
        numRow = self.numRow
        numLoci = self.numLoci
        stat5 = 0
        totalHomoCount = 0

        for i in range(numCol):
            homoCount = 0

            for j in (data[i]):
                if (j == '0101' or j == '0202' or j == '0303' or j == '0404'):
                    homoCount = homoCount + 1
            homoCount = homoCount / (numRow)
            totalHomoCount = homoCount * homoCount
            temp = 1 - totalHomoCount
            stat5 = stat5 + temp

        stat5 = stat5 / numLoci
        self.stat5 = stat5
        print("(Expected heterozygosity) stat5 is ", stat5)
