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
        
3. Python 3.8 or later is required to run the program

+++++++++++++++++++++++++++

ARUGMENTS

usage: python main [--s number of trails] [--o input]
```
positional arguments:
   input  input file name

optional arguments:
    --n     Flag for Monomorphic loci (default: False)
    --m     Minimum Allele Frequency (size: 0-1)
    --r     Mutation Rate (size: 0-1)
    --lNe   Lower of Ne Range (size: 10-)
    --uNe   Upper of Ne Range (size: -500)
    --lT    Lower of Theta Range (size: 1-)
    --uT    Upper of Theta Range (size: -10)
    --s     Number of OneSamp Trials (size: 1000-50000)
    --lD    Lower of Duration Range (size: 2-)
    --uD    Upper of Duration Range (size: -8)
    --i     Missing data for individuals (size: 0-1)
    --l     Missing data for loci (size: 0-1)
```

HOW TO RUN

1. Clone or download the repository

        git clone git@github.com:AaronHong1024/OneSamp.git
        cd <OneSamp file address>

2. Give permission to the OneSamp file under the build directory (We are working on applying the gene simulator source to this project)
        
        chmod 777 build/OneSamp
        
3. Go to configurations to add arguments (at least 1,000 trails and the loci size should be smaller than 5,000) 

   Example of minimum argument needed:
        
        --s 1000 --o exampleData/genePop10Ix30L

4. Run the program

        python main --s 1000 --o exampleData/genePop10Ix30L > output.txt


Feel free to add any parameters you want to use. We are working with adding this program to bioconda and docker, It will be much easier to use in future.


 
