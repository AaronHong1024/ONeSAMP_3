# OneSamp3.0
ONeSAMP 3.0

By XXXXXXX

+++++++++++++++++++++++++++

ABSTRACT

ONeSAMP 3.0 computes the effective population size of gene data sets. The estimation of the genetic effective 
population size (Ne) has been an established challenge in the field of population genetics. This program uses local linear regression and five 
summary statistics to conclude the Ne value. The software will take a file in GENEPOP format and provide the user with an effective population size. 

+++++++++++++++++++++++++++

USAGE OVERVIEW

1. Must have R downloaded in order to run the software
        
2. Python 3.8 or later

+++++++++++++++++++++++++++

ARUGMENTS

When using ONeSAMP 3.0, the user must input the following arguments at a minimum: **the number of trials and a GENEPOP file**. 

The following are additional arguments with their bounds: 

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

HOW TO EXECUTE

1. Set up R and Python environments

   How to set up R environment:
        https://linuxize.com/post/how-to-install-r-on-ubuntu-20-04/
        
   How to set up Python environment:
        https://www.digitalocean.com/community/tutorials/how-to-install-python-3-and-set-up-a-programming-environment-on-ubuntu-20-04-quickstart
        
2. Clone or download the repository using the following command:


        git clone git@github.com:IshaYoo/OneSamp3.0.git
        cd /OneSamp3.0
        chmod 777 build/OneSamp
        
3. To run ONeSAMP3.0:


        python3 ./main.py --s 2 --o exampleData/genePop5Ix5L
        
   Adjust configurations accordingly

   Example of minimum argument needed:
        
        --s 2 --o genePop5Ix5L
   (Running 2 trials with GENEPOP file names "genePop5Ix5L")
        
   The example GENEPOP file is provided under the folder “exampleData”

 
