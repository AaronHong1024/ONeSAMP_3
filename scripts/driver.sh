#!/bin/bash

###########################
###########################
## Simulation Parameters ##
###########################
###########################

############
# simcoal2 #
############

# Number of trials to run simcoal2 (and simuPOP)
export NUMTRIALS=50

# List of different Ne values to be used in simulation. Must have 5 digits apiece.
export NeVals="00256"
# 00128 00256 00512 01024
#export NeVals="00100"
export NeDigits=5

# Do not modify; used for cluster progress calculation
export NeNumberOfValues=`echo $NeVals | wc -w`

# Total number of loci to start with inside simcoal2 (max may be about 500)
#export initialLoci=500
  
# Inside simcoal2, we will use this many diploid individuals.
#export simcoalsample=20
  
# Large presumed historical effective population size in simcoal2.
#export simcoaldeme=1200

# UPDATE

# Use these constants for the simulation instead

#export ONESAMP2COAL_THETA=0.0001
#export ONESAMP2COAL_LOCI=2500
#export ONESAMP2COAL_DIPLOIDINDIVIDUALS=50

###########
# simuPOP #
###########

# simuPOP number of generations
#export simuPOPgens=4

###########
# ONeSAMP #
###########

# Number of individuals simuPOP should input in to Ne estimators. (500)
#export outputSampleSize=100
export outputSampleSize=50

# Number of loci in simulated trials. Only needed if simulating own trials.
export loci=500

# Minimum allele frequency necessary in populations
ONESAMP2COAL_MINALLELEFREQUENCY=0.05

# Mutation rate specified for simuPOP and ONeSAMP
export mutationRate="0.000000012"

# We require a flat distribution of values to supply to the executable.
export rangeNe=100,500

# We require a specified range of theta values to supply to the executable.
export theta=0.000048,0.0048

# Number of populations that ONeSAMP generates
export numOneSampTrials=20000

# Duration of bottlenecks
export duration=2,8

# Unix niceness of processes used by ONeSAMP: 19 is low priority, 0 is normal priority
export processPriority=19

# Unix niceness of file I/O: 7 is low priority, 0 is normal
export filePriority=7

# Do not modify; used for checking length of input files.
export trialsplus1=$(($numOneSampTrials + 1))

# SNPs or microsats: s for SNPs, m for microsatellites
export microsatsOrSNPs=s

# The block size parameter below defines how many iterations ONeSAMP performs
# in one trial. To reduce RAM usage, reduce the block size, but writes to disk
# will increase.
export blocksize=50

# The transfer times parameter below defines the ratio of local disk writes to
# NFS writes. If the ratio is lower, NFS traffic is reduced. If the ratio is
# higher, the progress indcator will take longer to update and the simulation
# may be slightly slower. Rare synchronization errors may also occur if the
# size of the cluster is large, then number of simutions is small, and this
# parameter is small. To reduce the chances of simulation errors, raise this
# parameter. The product of the blocksize and transfertimes variables should
# ideally be 500 or more unless 5 or fewer machines in the cluster are used.
# Empirically edit this parameter and use network diagnostic tools to properly
# set it. So long as the network traffic is low enough, it can be set to
# anything.
export transfertimes=80

#########################
# Export some constants #
#########################
export suffix=".reduced"
export pop="pop"
export number="number"
export gen=".gen"
export PARAMETER=".par"
export ARPSUFFIX="_0.arp"
export DISTSUFFIX=".out"
export LDNeSUFFIX=".genLD.txt"
export myLogin=`whoami`
export hostx=`hostname`
#export numLoggedIn=`who | cut -f1 -d\  | uniq | wc -l`

# Do not modify
export popStrLen=${#pop}
export suffixLen=${#suffix}
export DISTSUFFIXLen=${#DISTSUFFIX}
export genLen=${#gen}

############################################################################
# EDITING STEP 1: enter the filename of this file
export DRIVERNAME=driver.sh
############################################################################

# enter a python path for the simulator that contains the
# appropriate simuPOP library
#export PYTHONPATH=`echo ~`/simupop/lib64/python2.7/site-packages/simuPOP-1.1.6-py2.7-linux-x86_64.egg

# Do not modify
export CURRENTPATH="$( cd "$(dirname "$0")" ; pwd -P )"

############################################################################
# EDITING STEP 2: if you are using a cluster setup, please supply local
# paths to be used for temporary file storage /home/sharique/Downloads/ONeSAMP2/refactor/release/tmp
#export in_bufferdir=/s/${hostx}/a/tmp/${myLogin}_tmp_in_filebuffer
#export out_bufferdir=/s/${hostx}/a/tmp/${myLogin}_tmp_out_filebuffer
#export singleton_out_bufferdir=/s/${hostx}/a/tmp/${myLogin}_tmp_filebuffer_safe
touch temp_in
touch temp_out
touch temp_safe
export in_bufferdir=temp_in
export out_bufferdir=temp_out
export singleton_out_bufferdir=temp_safe
############################################################################

# Do not use
#export SIMCOAL2EXECUTABLE=`echo ~`/simcoal2/simcoal2_1_2

# Do not modify 
export hostFile=$3

##################################
# Define a few programming tools #
##################################

############################################################################
# EDITING STEP 3: please supply the locations of the python, perl, and R
# interpreters, as well as the C and C++ compilers to be used
export PYTHONEXECUTABLE=python
export PERLEXECUTABLE='perl -I/usr/lib64/perl5/vendor_perl -I/usr/share/perl5/vendor_perl'
export RINTERPRETER=Rscript
export CCCOMPILE=g++
export CCOMPILE=gcc
############################################################################

############################################################################
# EDITING STEP 4: please supply the locations of the built ONeSAMP2 executable
# and the coalescent simulator, as well as the ONESAMP2 executable name
export ONESAMP2EXEC=refactor_main
export ONESAMP2=../release/${ONESAMP2EXEC}
export ONESAMP2COAL=./refactor_coalescent_simulator
export RSCRIPT=rScript.r
############################################################################

#######################
# Usage error message #
#######################

if [ $# -eq 0 ]; then
  echo "Usage: Please use one of these sets of options"
  echo "./$DRIVERNAME generateSimulatedPopulationsOfKnownNe"
  echo "./$DRIVERNAME analyzeUnknownPopulation"
  echo "./$DRIVERNAME runOneSamp start"
  echo "./$DRIVERNAME runOneSamp progress"
  echo "./$DRIVERNAME estimateNe"
  echo "./$DRIVERNAME runLDNe"
  echo "./$DRIVERNAME onesampDataAnalytics outFile"
fi



############
############
## STEP 1 ##
############
############

# Generates populations of known Ne
if [ $1 ]; then
  if [ $1 == "generateSimulatedPopulationsOfKnownNe" ]; then
    # This bash script employs Perl, Python and C++ to generate
	# example populations
    # of known effective population size on the local machine.
    
    # Substep: Set the number of trials and the Ne values to generate
    
    # For each combination of trials
    
    for numPOP in $( echo $NeVals ); do
      echo "Generating all test cases for Ne = " $numPOP
      for trial in `seq -w 999 | head -n $NUMTRIALS`; do
        echo "TRIAL "$trial

        # Needed to prevent race condition in ONeSAMP random number initialization
        sleep 2

        #echo "Running SIMCOAL2..."
        echo "Running ONESAMP2 coalescent simulator..."
        # Substep: Establish temporary names
        export genepopdata=ONeSAMP1Temp${numPOP}number${trial}
        export simcoalout=ONeSAMP2Temp${numPOP}number${trial}
        export simupopin=ONeSAMP3Temp${numPOP}number${trial}${gen}
        export executablename=ONeSAMP4Temp${numPOP}number${trial}
        export pythonscript=ONeSAMP5Temp${numPOP}number${trial}.py
        export OUTPUT=${pop}${numPOP}${number}${trial}${gen}
    
        # Substep: Establish simulation parameters
  
        # For generating the expected number of individuals.
        export expectedNe=$numPOP
    
        # Substep: Clean up temporary files

        rm -fr $simcoalout
        rm -fr $simcoalout$PARAMETER
        rm -fr $genepopdata
        rm -fr $executablename
        rm -fr $pythonscript
        rm -fr $OUTPUT
    
        # Substep: Make simcoal2 input file (edit: do not use)
    
#        echo "//Input parameters for the coalescence and recombination simulation program : simcoal2.exe" > $simcoalout$PARAMETER
#        echo "1 samples to simulate" >> $simcoalout$PARAMETER
#        echo "//Deme sizes (haploid number of genes)" >> $simcoalout$PARAMETER
#        echo $simcoaldeme >> $simcoalout$PARAMETER
#        echo "//Sample sizes" >> $simcoalout$PARAMETER
#        echo $(($simcoalsample + $simcoalsample)) >> $simcoalout$PARAMETER
#        echo "//Growth rates" >> $simcoalout$PARAMETER
#        echo "0" >> $simcoalout$PARAMETER
#        echo "//Number of migration matrices : If 0 : No migration between demes" >> $simcoalout$PARAMETER
#        echo "0" >> $simcoalout$PARAMETER
#        echo "//Historical event: time, source, sink, proportion of migrants, new deme size, new growth rate, new migration matrix" >> $simcoalout$PARAMETER
#        echo "0 historical events" >> $simcoalout$PARAMETER
#        echo "//Number of independent (unlinked) chromosomes, and \"chromosome structure\" flag:  0 for identical structure across chromosomes, and 1 for different structures on different chromosomes." >> $simcoalout$PARAMETER
#        echo "1 0" >> $simcoalout$PARAMETER
#        echo "//Number of contiguous linkage blocks in chromosome 1:"  >> $simcoalout$PARAMETER
#        echo "1"  >> $simcoalout$PARAMETER
#        echo "//Per Block: Data type, No. of loci, Recombination rate to the right-side locus, plus optional parameters ***see detailed explanation here***"  >> $simcoalout$PARAMETER
#        echo "SNP "$initialLoci" 0.5 0.05" >> $simcoalout$PARAMETER
#        echo "[EOF]" >> $simcoalout$PARAMETER

        # Substep: Run simcoal2 (edit: do not use)

#        (( echo $simcoalout ; echo 1; echo 1;) | $SIMCOAL2EXECUTABLE ) > /dev/null
    
        # Substep: Convert output of simcoal2 to GENEPOP format
    
#        (echo Auto-generated genepop output from simcoal2; seq $initialLoci; echo Pop; cat $simcoalout/$simcoalout$ARPSUFFIX | tail -n +22 | head -n -10 | $PERLEXECUTABLE -pe 's/\t//g' | paste -sd" \n" | $PERLEXECUTABLE -pe 's/\d//g' | $PERLEXECUTABLE -pe 's/_//g' | $PERLEXECUTABLE -pe 's/  /X/g' | $PERLEXECUTABLE -pe 'use List::MoreUtils qw(zip); my @a = split(/X/, $_); my @b = split(/ /, $a[0]); my @c = split(/ /, $a[1]); $_ = join("", zip(@b, @c));' | $PERLEXECUTABLE -pe 's/(A|C|G|T)(A|C|G|T)/$1$2 /g' | $PERLEXECUTABLE -pe 's/A/01/g' | $PERLEXECUTABLE -pe 's/C/02/g' | $PERLEXECUTABLE -pe 's/G/03/g' | $PERLEXECUTABLE -pe 's/T/04/g' | $PERLEXECUTABLE -pe 's/^/indiv , /g' ) > $genepopdata

        # UPDATE: Replace simcoal2 with a snippet of ONeSAMP2 code to ensure that
        # the assumptions in the frequency distribution are matched between test cases and code
        # UPDATE3:
        #$ONESAMP2COAL SNP $ONESAMP2COAL_THETA $ONESAMP2COAL_LOCI $ONESAMP2COAL_DIPLOIDINDIVIDUALS $ONESAMP2COAL_MINALLELEFREQUENCY > $genepopdata

        # Substep: Make new population from backwards simulator to plug into
        # forwards simulator by making new population with specified frequencies.
        # The number of individuals in the output file will be the number of input
        # individuals in Ne.

        # UPDATE 3: remove
    
        #( echo \#include\ \<iostream\>; echo \#include\ \<cstdio\>; echo \#include\ \<cstdlib\>; echo \#include\ \<ctime\>; echo ; echo \#define\ A\ 0; echo \#define\ C\ 1; echo \#define\ G\ 2; echo \#define\ T\ 3; echo ; echo \#define\ PLOIDY\ 2; echo ; echo using\ namespace\ std\;; echo ; echo int\ checkIfValidDigit\(char\ digit\)\{; echo \ \ return\ digit\ \=\=\ \'1\'\ \|\|\ digit\ \=\=\ \'2\'\ \|\|\ digit\ \=\=\ \'3\'\ \|\|\ digit\ \=\=\ \'4\'\;; echo \}; echo ; echo int\ checkIfValidGenotype\(char\ \*vals\)\{; echo \ \ if\(vals\[0\]\ \!\=\ \'0\'\ \|\|\ vals\[2\]\ \!\=\ \'0\'\)\ return\ \-1\;; echo \ \ if\(checkIfValidDigit\(vals\[1\]\)\ \&\&\ checkIfValidDigit\(vals\[3\]\)\)\ return\ atoi\(vals\)\;; echo \ \ return\ \-1\;; echo \}; echo ; echo \/\*; echo \ \*\ Accepts\ a\ SNP\ GENEPOP\ file\,\ tallies\ frequencies\ for\ each\ locus\,\ and\ prints; echo \ \*\ out\ a\ new\ SNP\ GENEPOP\ file\ with\ specified\ number\ of\ individuals\ and\ a\ sampled; echo \ \*\ number\ of\ individuals\ at\ each\ locus\ with\ frequencies\ specified\ by\ the; echo \ \*\ tallies\.; echo \ \*; echo \ \*\ TODO\:\ Currently\ completes\ minimal\ error\ checking\.; echo \ \*\/; echo int\ main\(int\ argc\,\ char\ \*argv\[\]\)\{; echo \ \ srand\(time\(NULL\)\)\;; echo \ \ int\ i\ \=\ 0\;; echo \ \ int\ j\ \=\ 0\;; echo \ \ int\ k\ \=\ 0\;; echo ; echo \ \ if\(argc\ \!\=\ 4\)\{; echo \ \ \ \ cout\ \<\<\ \"Usage\:\ \.\/\"\ \<\<\ argv\[0\]\ \<\<\ \"\(num\ input\ individuals\)\ \(num\ output\ individuals\)\ \(num\ loci\)\ \<\ \(input\ file\)\\n\"\;; echo \ \ \ \ return\ 1\;; echo \ \ \}; echo ; echo \ \ int\ indivs\ \=\ PLOIDY\ \*\ atoi\(argv\[1\]\)\;; echo \ \ int\ outputIndivs\ \=\ atoi\(argv\[2\]\)\;; echo \ \ int\ loci\ \=\ atoi\(argv\[3\]\)\;; echo ; echo \ \ int\ \*\*counts\ \=\ \(int\ \*\*\)\ malloc\(loci\ \*\ sizeof\(int\ \*\)\)\;; echo \ \ int\ \*\*genes\ \=\ \(int\ \*\*\)\ malloc\(loci\ \*\ sizeof\(int\ \*\)\)\;; echo ; echo \ \ for\(i\ \=\ 0\;\ i\ \<\ loci\;\ i\+\+\)\{; echo \ \ \ \ counts\[i\]\ \=\ \(int\ \*\)\ malloc\(4\ \*\ sizeof\(int\)\)\;; echo \ \ \ \ \/\/floats\[i\]\ \=\ \(double\ \*\)\ malloc\(4\ \*\ sizeof\(double\)\)\;; echo \ \ \ \ counts\[i\]\[A\]\ \=\ 0\;; echo \ \ \ \ counts\[i\]\[C\]\ \=\ 0\;; echo \ \ \ \ counts\[i\]\[G\]\ \=\ 0\;; echo \ \ \ \ counts\[i\]\[T\]\ \=\ 0\;; echo \ \ \ \ genes\[i\]\ \=\ \(int\ \*\)\ malloc\(PLOIDY\ \*\ outputIndivs\ \*\ sizeof\(int\)\)\;; echo \ \ \}; echo ; echo \ \ char\ c\;; echo ; echo \ \ while\(cin\.get\(c\)\)\{; echo \ \ \ \ if\(c\ \=\=\ \'\\n\'\)\ break\;; echo \ \ \}; echo ; echo \ \ char\ first\ \=\ \'\ \'\;; echo \ \ char\ second\ \=\ \'\ \'\;; echo \ \ char\ third\ \=\ \'\ \'\;; echo ; echo \ \ while\(cin\.get\(c\)\)\{; echo \ \ \ \ first\ \=\ second\;; echo \ \ \ \ second\ \=\ third\;; echo \ \ \ \ third\ \=\ c\;; echo \ \ \ \ if\(\(first\ \=\=\ \'P\'\ \|\|\ first\ \=\=\ \'p\'\)\ \&\&\ \(second\ \=\=\ \'O\'\ \|\|\ second\ \=\=\ \'o\'\)\ \&\&\ \(third\ \=\=\ \'P\'\ \|\|\ third\ \=\=\ \'p\'\)\)\ break\;; echo \ \ \}; echo ; echo \ \ char\ vals\[4\]\ \=\ \{\'X\'\,\'X\'\,\'X\'\,\'X\'\}\;; echo ; echo \ \ i\ \=\ 0\;\ \/\/\(loci\ counter\); echo \ \ k\ \=\ 0\;\ \/\/\(genotype\); echo ; echo \ \ while\(cin\.get\(c\)\)\{; echo \ \ \ \ vals\[0\]\ \=\ vals\[1\]\;; echo \ \ \ \ vals\[1\]\ \=\ vals\[2\]\;; echo \ \ \ \ vals\[2\]\ \=\ vals\[3\]\;; echo \ \ \ \ vals\[3\]\ \=\ c\;; echo ; echo \ \ \ \ \/\/\ Count\ number\ of\ As\,\ Cs\,\ Gs\,\ and\ Ts\.; echo \ \ \ \ if\(\(k\ \=\ checkIfValidGenotype\(vals\)\)\ \!\=\ \-1\)\{; echo \ \ \ \ \ \ counts\[i\]\[\(k\ \-\ 1\)\ \%\ 100\]\+\+\;; echo \ \ \ \ \ \ counts\[i\]\[\(k\ \-\ 100\)\ \/\ 100\]\+\+\;; echo \ \ \ \ \ \ i\ \=\ \(i\ \+\ 1\)\ \%\ loci\;; echo \ \ \ \ \}; echo \ \ \}; echo ; echo \ \ for\(i\ \=\ 0\;\ i\ \<\ loci\;\ i\+\+\)\{; echo \ \ \ \ for\(j\ \=\ 0\;\ j\ \<\ PLOIDY\ \*\ outputIndivs\;\ j\+\+\)\{; echo \ \ \ \ \ \ int\ val\ \=\ rand\(\)\ \%\ indivs\;; echo \ \ \ \ \ \ if\(val\ \<\ counts\[i\]\[A\]\)\ genes\[i\]\[j\]\ \=\ A\;; echo \ \ \ \ \ \ else\ if\(val\ \<\ counts\[i\]\[A\]\ \+\ counts\[i\]\[C\]\)\ genes\[i\]\[j\]\ \=\ C\;; echo \ \ \ \ \ \ else\ if\(val\ \<\ counts\[i\]\[A\]\ \+\ counts\[i\]\[C\]\ \+\ counts\[i\]\[G\]\)\ genes\[i\]\[j\]\ \=\ G\;; echo \ \ \ \ \ \ else\ genes\[i\]\[j\]\ \=\ T\;; echo \ \ \ \ \}; echo \ \ \}; echo ; echo \ \ cout\ \<\<\ \"Auto\-generated\ output\ post\-processed\ from\ simcoal2\\n\"\;; echo \ \ for\(i\ \=\ 0\;\ i\ \<\ loci\;\ i\+\+\)\{; echo \ \ \ \ cout\ \<\<\ i\ \<\<\ \"\\n\"\;; echo \ \ \}; echo ; echo \ \ cout\ \<\<\ \"Pop\\n\"\;; echo ; echo \ \ for\(j\ \=\ 0\;\ j\ \<\ outputIndivs\;\ j\+\+\)\{; echo \ \ \ \ cout\ \<\<\ \"indiv\"\ \<\<\ j\ \+\ 1\ \<\<\ \"\,\"\;; echo \ \ \ \ for\(i\ \=\ 0\;\ i\ \<\ loci\;\ i\+\+\)\{; echo \ \ \ \ \ \ cout\ \<\<\ \"\ \"\;; echo \ \ \ \ \ \ for\(k\ \=\ 0\;\ k\ \<\ PLOIDY\;\ k\+\+\)\{; echo \ \ \ \ \ \ \ \ printf\(\"\%02d\"\,\ genes\[i\]\[2\ \*\ j\ \+\ k\]\ \+\ 1\)\;; echo \ \ \ \ \ \ \}; echo \ \ \ \ \}; echo \ \ \ \ cout\ \<\<\ \"\\n\"\;; echo \ \ \}; echo ; echo \ \ return\ 0\;; echo \}) | $CCCOMPILE -x c++ - -o $executablename

        # UPDATE3tmp_in
        #chmod u+x $executablename
        #./$executablename $simcoalsample $expectedNe $initialLoci < $genepopdata > $simupopin
    
        # UPDATE use the number of loci fetched from ONESAMP coalescent simulator
        # ./$executablename $ONESAMP2COAL_DIPLOIDINDIVIDUALS $expectedNe $ONESAMP2COAL_LOCI < $genepopdata > $simupopin

        # UPDATE2 no longer perform the sampling step, simply copy the genotypes over
        #cp $genepopdata $simupopin

        # Substep: Run the input file in the forward simulator for a fixed small
        # number of generations in this python script.
        # UPDATE3
        #echo "Running simuPOP..."
        #( echo \#\ The\ purpose\ of\ this\ script\ is\ to\ generate\ test\ cases\ by\ which\ ONeSAMP\ 2\.0\ can\ be\ tested\.; echo \#\ August\ 19\,\ 2015; echo import\ simuPOP\ as\ simulator; echo from\ simuPOP\.utils\ import\ importPopulation\,\ export; echo import\ sys; echo print\(\"To\ use\ this\ script\ to\ generate\ populations\,\ run\ as\:\ \\\\n\ python\ temp\ inputFile\ numOutputIndividuals\ numOutputLoci\ mutationRate\ outputFile\\\\n\"\); echo \#\ Read\ command\ line\ arguments; echo if\(len\(sys\.argv\)\ \!\=\ 6\)\:; echo \ \ print\(\"Incorrect\ number\ of\ arguments\ specified\.\"\); echo inFile\=sys\.argv\[1\]; echo outputIndivs\=int\(sys\.argv\[2\]\); echo maxLoci\=int\(sys\.argv\[3\]\)\ \-\ 1; echo mRate\=float\(sys\.argv\[4\]\); echo outFile\=sys\.argv\[5\]; echo \#\ Initialize\ constraints; echo \#pop\ \=\ simulator\.Population\(size\=\[individuals\]\,\ ploidy\=2\,\ loci\=nLoc\,\ alleleNames\=\[\'A\'\,\ \'C\'\,\ \'G\'\,\ \'T\'\]\); echo pop\ \=\ simulator\.utils\.importPopulation\(format\=\'GENEPOP\'\,\ filename\=inFile\); echo \#\ initOps\ \=\ \[simulator\.InitSex\(\)\,\ simulator\.InitGenotype\(freq\=\[\.25\,\ \.25\,\ \.25\,\ \.25\]\)\]; echo numGens\ \=\ $simuPOPgens\;; echo \#\ Simulate\ generations; echo pop\.evolve\(initOps\=\[simulator\.InitSex\(\)\]\,\ matingScheme\=simulator\.RandomMating\(ops\=simulator\.Recombinator\(rates\=0\.5\)\)\,\ preOps\=\[simulator\.AcgtMutator\(rate\=\[mRate\]\,\ model\=\'JC69\'\)\,\ simulator\.ResizeSubPops\(sizes\=\[outputIndivs\]\,\ at\=numGens\-1\)\]\,\ postOps\=\[simulator\.Stat\(alleleFreq\=simulator\.ALL\_AVAIL\)\]\,\ gen\=numGens\); echo monomorphs\ \=\ list\(\); echo \#\ Detect\ monomorphic\ loci; echo for\ i\ in\ pop\.dvars\(\)\.alleleFreq\:; echo \ \ if\ pop\.dvars\(\)\.alleleFreq\[i\]\.values\(\)\[0\]\ \=\=\ 1\:; echo \ \ \ \ \ monomorphs\.append\(i\); echo \#\ Remove\ monomorphic\ loci; echo pop\.removeLoci\(monomorphs\); echo \#\ Detect\ extra\ loci; echo extraLoci\ \=\ list\(\); echo for\ i\ in\ pop\.dvars\(\)\.alleleFreq\:; echo \ \ if\ i\ \>\ maxLoci\:; echo \ \ \ \ extraLoci\.append\(i\); echo \#\ Remove\ extra\ loci; echo pop\.removeLoci\(extraLoci\); echo \#\ Save\ population; echo export\(pop\,\ format\=\'genepop\'\,\ output\=outFile\,\ adjust\=0\)) > $pythonscript

        # UPDATE3 chmod u+x $pythonscript

        # Substep: Run python script

        # UPDATE3 $PYTHONEXECUTABLE $pythonscript $simupopin $numPOP $loci $mutationRate $OUTPUT

        # Substep: Clean up temporary files

        rm -fr $ccprog
        rm -fr $simupopin
        rm -fr $simcoalout
        rm -fr $simcoalout$PARAMETER
        rm -fr $genepopdata
        rm -fr $executablename
        rm -fr $pythonscript
    
        # Substep: Reduce the number of samples in the population to a fixed value
        # UPDATE 3
        #echo "Running ONeSAMP population reduction script... "
        #nice -n $processPriority ionice -n $filePriority $ONESAMP2 -u$mutationRate -v$theta -rC -l$loci -i$outputSampleSize -b$reducedSize -d$duration -$microsatsOrSNPs -t1 -f$ONESAMP2COAL_MINALLELEFREQUENCY -o1 -g < $OUTPUT > $OUTPUT$suffix

        # UPDATE 2 simply copy over, throw out extra individuals

        # This step restricts to $reducedSize individuals
        # UPDATE 3
        # export nLines=$((reducedSize+3))
        # UPDATE 3
        #cat $OUTPUT | head -n$nLines > $OUTPUT$suffix

        nice -n $processPriority ionice -n $filePriority $ONESAMP2 -t1 -rC -b$numPOP -d1 -u$mutationRate -v${theta} -$microsatsOrSNPs -l$loci -i$outputSampleSize -o1 -f$ONESAMP2COAL_MINALLELEFREQUENCY -p > $OUTPUT$suffix

        echo "Trial complete"
      done
    done
    
    # Substep: Construct the first line of the analysis file that computes the stats for the original population
    
    for j in `ls *$suffix | cat`; do
      echo "Generating first line in analysis file for "$j
      # Fetch -l and -i flags first
      export flags=`nice -n $processPriority ionice -n $filePriority $ONESAMP2 -x < $j | tr -d '\n'`
      nice -n $processPriority ionice -n $filePriority $ONESAMP2 < $j $flags -rC -$microsatsOrSNPs -d$duration -b$rangeNe -v$theta -u$mutationRate -t1 -f$ONESAMP2COAL_MINALLELEFREQUENCY -o1 -a -w > ${j}$DISTSUFFIX
    done
    
  fi
fi
# Generates first line of analysis file
if [ $1 ]; then
  if [ $1 == "analyzeUnknownPopulation" ]; then
    for j in `ls *$suffix | cat`; do
      echo "Generating first line in analysis file for "$j
      # UPDATE 3
      export flags=`nice -n $processPriority ionice -n $filePriority $ONESAMP2 -x < $j | tr -d '\n'`
      nice -n $processPriority ionice -n $filePriority $ONESAMP2 < $j $flags -rC -$microsatsOrSNPs -d$duration -b$rangeNe -v$theta -u$mutationRate -t1 -f$ONESAMP2COAL_MINALLELEFREQUENCY -o1 -a -w > ${j}$DISTSUFFIX
    done
    
  fi
fi

############
############
## STEP 2 ##
############
############

# Simulate several rounds of ONeSAMP on the cluster using the basic commands
if [ $1 ]; then
  if [ $1 == "runOneSamp" ]; then
     if [ $2 ]; then
       if [ $2 == "start" ]; then
        nice -n $processPriority $CURRENTPATH/driver.sh runOneSamp internal actuator $4 &
       fi

      if [ $2 == "progress" ]; then
        echo "Computing progress on current instance, please wait..."
        export totalPercent=0
        for i in `ls *${gen}${suffix} | shuf`;
        do
          export sizeoffile=`((nice -n $processPriority grep -v ^$ | nice -n $processPriority wc -l) < ${i}$DISTSUFFIX)`
          export percentComplete=`echo "100 * $sizeoffile / $numOneSampTrials / $NUMTRIALS / $NeNumberOfValues" | bc -l`
          export totalPercent=`echo "$percentComplete + $totalPercent" | bc -l`;
        done
        printf "Current progress = %.1f%%\n" $totalPercent
      fi

      if [ $2 == "stop" ]; then
          echo "Shutting down task"
          rm $in_bufferdir
          rm $out_bufferdir
          rm $singleton_out_bufferdir
          export tmpKillSchedCommand="\`ps aux | grep ${myLogin} | grep $DRIVERNAME | grep -v grep | grep -v kill | awk '{ print \$2 }'\`"
          kill $tmpKillSchedCommand
          export tmpKillExecutableCommand="\`ps aux | grep ${myLogin} | grep $ONESAMP2EXEC | grep -v grep | grep -v kill | awk '{ print \$2 }'\`"
          kill $tmpKillExecutableCommand
      fi
     

      # Internal actuator script: 

      if [ $2 == "internal" ]; then
        if [ $3 == "actuator" ]; then
	
          cd $CURRENTPATH
          stop=false
          while [ $stop = false ];
          do

            stop=true
            # of gene files in random order
            for i in `ls *${gen}${suffix} | shuf`;
            do

              sizeoffile=`((nice -n $processPriority grep -v ^$ | nice -n $processPriority wc -l) < ${i}$DISTSUFFIX)`

              if [ $sizeoffile -lt $trialsplus1 ]; then

                stop=false

                # Create temp files for IO
                ionice -n $filePriority rm $in_bufferdir
                ionice -n $filePriority touch $in_bufferdir
                ionice -n $filePriority chmod 600 $in_bufferdir
                ionice -n $filePriority cp $i $in_bufferdir
                ionice -n $filePriority rm $out_bufferdir
                ionice -n $filePriority touch $out_bufferdir
                ionice -n $filePriority chmod 600 $out_bufferdir
    
                # Execute one group of iterations (overwrite temp file)

                # Fetch -l and -i flags first
                export flags=`nice -n $processPriority ionice -n $filePriority $ONESAMP2 -x < $in_bufferdir | tr -d '\n'`
                nice -n $processPriority ionice -n $filePriority $ONESAMP2 < $in_bufferdir $flags -rC -$microsatsOrSNPs -d$duration -b$rangeNe -v$theta -u$mutationRate -t$blocksize -f$ONESAMP2COAL_MINALLELEFREQUENCY -o1 -a -e > $out_bufferdir

                # Execute the rest of the iterations (append to temp file)
    
                for q in `seq 2 $transfertimes`
                do
                  # Fetch -l and -i flags first
                  export flags=`nice -n $processPriority ionice -n $filePriority $ONESAMP2 -x < $in_bufferdir | tr -d '\n'`
                  nice -n $processPriority ionice -n $filePriority $ONESAMP2 < $in_bufferdir $flags -rC -$microsatsOrSNPs -d$duration -b$rangeNe -v$theta -u$mutationRate -t$blocksize -f$ONESAMP2COAL_MINALLELEFREQUENCY -o1 -a -e >> $out_bufferdir
                done

                # Copy trials over NFS


                ionice -n $filePriority cat $out_bufferdir >> ${i}$DISTSUFFIX
                ionice -n $filePriority rm $out_bufferdir
              fi

              # If we have enough trials, end the simulations here

              if [ $sizeoffile -gt $trialsplus1 ]; then
                touch $singleton_out_bufferdir
                chmod 600 $singleton_out_bufferdir
                ionice -n $filePriority cat ${i}$DISTSUFFIX > $singleton_out_bufferdir
                ionice -n $filePriority head -n$trialsplus1 ${i}$DISTSUFFIX > $singleton_out_bufferdir
               (nice -n $processPriority ionice -n $filePriority grep -v ^$ | nice -n $processPriority head -n$trialsplus1) > ${i}$DISTSUFFIX < $singleton_out_bufferdir
                ionice -n $filePriority rm $singleton_out_bufferdir
              fi

              # Delete local trials
              ionice -n $filePriority rm $in_bufferdir
              ionice -n $filePriority rm $out_bufferdir
            done
          done
        fi
      fi
    fi
  fi
fi
############
############
## STEP 3 ##
############
############

# Create file with these columns:
# Known Ne
# ONeSAMP Ne mean
# ONeSAMP Ne median
# ONeSAMP lower confidence value
# ONeSAMP upper confidence value
# LDNe estimate
# LDNe lower confidence bound
# LDNe upper confidence bound

if [ $1 ]; then
  if [ $1 == "onesampDataAnalytics" ]; then
    if [ $2 ]; then
      rm -fr $2
      for dataPoint in `ls *${DISTSUFFIX} | cat`; do
      export prefLength=${#dataPoint}
      export prefLength=$((prefLength-DISTSUFFIXLen))
      export prefLength=$((prefLength-suffixLen))
      export prefLength=$((prefLength-genLen))
      export LDNeFile=${dataPoint:0:$prefLength}${LDNeSUFFIX}
      touch $LDNeFile
      # Run R on the ONeSAMP output and extract;
      # Also extract LDNe output
      ( (echo -n "${dataPoint:popStrLen:NeDigits} " ; (echo -n `${RINTERPRETER} ${RSCRIPT} < $dataPoint | tail -n1`) 2> /dev/null ; echo -n ' '; echo -n `head -n21 $LDNeFile | tail -n1 | awk '{print $4}'`; echo -n ' '; echo -n `head -n24 $LDNeFile | tail -n1 | awk '{print $3}'`; echo -n ' '; echo `head -n25 $LDNeFile | tail -n1 | awk '{print $1}'`)) >> $2
      done
    fi
  fi
fi

if [ $1 ]; then
  if [ $1 == "estimateNe" ]; then
      rm result.txt
      touch result.txt
      echo "Mean Median 95%Min 95%Max" >> result.txt
      for dataPoint in `ls *${DISTSUFFIX} | cat`; do
      # Run R on the ONeSAMP output and extract;
      ( (echo `${RINTERPRETER} ${RSCRIPT} < $dataPoint | tail -n1`) 2> /dev/null) >> result.txt
      done
  fi
fi

############
############
## STEP 4 ##
############
############

if [ $1 ]; then
  if [ $1 == "runLDNe" ]; then
    #echo "LDNe was removed in this release..."
    #exit

# Run LDNe, from Waples and Do

#######################################
# BEGIN CODE INSERTION: WAPLES AND DO #
#######################################

#####################################
# END CODE INSERTION: WAPLES AND DO #
#####################################

    for dataPoint in `ls *${suffix} | cat`; do
    	export prefLength=${#dataPoint}
	export prefLength=$((prefLength-DISTSUFFIXLen))
	export prefLength=$((prefLength-suffixLen))
	export prefLength=$((prefLength-genLen))
	export LDNeFile=${dataPoint:0:$prefLength}${LDNeSUFFIX}
	touch $LDNeFile
      ( ( echo $dataPoint ; echo ; echo n ) | ./ldneEXEC )
    done
  fi
fi
