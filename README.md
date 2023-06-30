# OneSamp3.0
ONeSAMP 3.0

By Aaron Hong1 Kingshuk Mukherjee1 Isha Yooseph1 Marco Oliva1
Mark Heim2 Chris Funk3 David Tallmon4,† Christina Boucher



+++++++++++++++++++++++++++

ABSTRACT

ONeSAMP 3.0 computes the effective population size of gene data sets. The estimation of the genetic effective 
population size (Ne) has been an established challenge in the field of population genetics. This program uses local linear regression and five 
summary statistics to conclude the Ne value. The software will take a file in GENEPOP format and provide the user with an effective population size. 

+++++++++++++++++++++++++++

USAGE OVERVIEW

1. Must have R downloaded in order to run the software
        
   You can download and set up the R environment at this link: 
        
        https://www.tutorialspoint.com/r/r_environment_setup.htm
        
3. Python 3.8 or later

+++++++++++++++++++++++++++

ARUGMENTS

When using ONeSAMP 3.0, the user must input the following arguments at a minimum: **the number of trials and a GENEPOP file**. 

The following are additional arguments with their bounds: 

**--n Flag for Monomorphic loci**
default: False

**--m Minimum Allele Frequency**

Lower Bound: 0

Upper Bound: 1

**--r Mutation Rate**

Lower Bound: 0

Upper Bound: 1

**--lNe Lower of Ne Range**

Lower Bound: 10

Upper Bound: -

**--uNe Upper of Ne Range**

Lower Bound: -

Upper Bound: 500

**--lT Lower of Theta Range**

Lower Bound: 1

Upper Bound: -

**--uT Upper of Theta Range**

Lower Bound: -

Upper Bound: 10

**--s Number of OneSamp Trials**

Lower Bound: 2

Upper Bound: 50,000

**--lD Lower of Duration Range**

Lower Bound: 2

Upper Bound: -

**--uD Upper of Duration Range**

Lower Bound: -

Upper Bound: 8

**--i Missing data for individuals**

Lower Bound: 0

Upper Bound: 1

**--l Missing data for loci**

Lower Bound: 0

Upper Bound: 1

**--o The File Name**

Lower Bound: -

Upper Bound: -

+++++++++++++++++++++++++++

HOW TO RUN

1. Clone or download the repository
        
        cd <OneSamp file address>

2. Give permission to the OneSamp file under the build directory (We are working on applying the gene simulator source to this project)
        
        chmod 777 build/OneSamp
        
3. Go to configurations to add arguments (at least 1000 trails) 

   Example of minimum argument needed:
        
        --s 1000 --o exampleData/genePop10Ix30L
        
Full command will be:

python main --s 1000 --o exampleData/genePop10Ix30L

The input dataset should have a locus size smaller than 5 thousand, Or it will be super time-consuming.

Feel free to add any parameters you want to use. We are working with adding this program to bioconda and docker, It will be much easier to use in future.


 
