
# ONeSAMP3.0  ![GitHub Release](https://img.shields.io/github/v/release/AaronHong1024/ONeSAMP_3) ![GitHub License](https://img.shields.io/github/license/AaronHong1024/ONeSAMP_3)


ONeSAMP 3.0 computes the effective population size of genomic data sets.
This program takes a file in GENEPOP format and computes five summary statistics. 
The software then uses linear regression based on these summary statistics to estimate of effective population size.  

It is strongly recommended that users read the accompanying manuscript before applying ONeSAMP to their data. 



# Usage Overview
1. The program can be executed on MAC OS or LINUX system.

2. Must have R downloaded in order to run the software
        
   You can download and set up the R environment at this link: 
        
        https://www.tutorialspoint.com/r/r_environment_setup.htm

3. Python 3.8 or later is required to run the program

# Getting Started
1. Make a new ONeSAMP directory

        mkdir ONeSAMP
        cd ONeSAMP
2. Clone the repository

        git clone git@github.com:AaronHong1024/OneSamp.git
3. Give the Permission to the ONeSAMP file under the build directory

        chmod 777 build/OneSamp

# Running ONeSAMP

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
    --s     Number of ONeSAMP Trials (size: 1000-50000)
    --lD    Lower of Duration Range (size: 2-)
    --uD    Upper of Duration Range (size: -8)
    --i     Missing data for individuals (size: 0-1)
    --l     Missing data for loci (size: 0-1)
```


Run the program

        python main --s 1000 --o exampleData/genePop10Ix30L > output.txt

# Getting Help

If you have any issues or anyquestions, please feel free to contact us at yu.hong@ufl.edu or through the GitHub Issues.

# Citing ONeSAMP
If you use the ONeSAMP in your research project, please cite: 

Hong, A., Cheek, R.G., Mukherjee, K., Yooseph, I., Oliva, M., Heim, M., Funk, W.C., Tallmon, D. and Boucher, C., 2023. ONeSAMP 3.0: Effective Population Size via SNP Data for One Population Sample.



 
