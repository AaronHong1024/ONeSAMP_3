# OneSamp3.0
ONeSAMP 3.0

By XXXXXXX


## ABSTRACT 

ONeSAMP 3.0 computes the effective population size of gene data sets. The estimation of the genetic effective 
population size (Ne) has been an established challenge in the field of population genetics. This program uses local linear regression and five 
summary statistics to conclude the Ne value. The software will take a file in GENEPOP format and provide the user with an effective population size. 


## DEPENDENCY

1. R 4.1 or later
        
2. Python 3.8 or later


## INPUT & OUTPUT



**INPUT**
##
The minimum argument needed is a GENEPOP file. Example GENEPOP files are provided under the folder “exampleData”. 

Example of minimum argument needed:
        
        --o exampleData/genePop5Ix5L
   (Running ONeSAMP3.0 with a GENEPOP file named "genePop5Ix5L" located under the exampleData folder)
   
   
   
**PARAMETERS**
##
When using ONeSAMP 3.0, the user must input **a GENEPOP file** at minimum. 

The following are additional arguments with their bounds and default value: 

**Minimum Allele Frequency**

Default value: 0.005 ** Lower Bound: 0 ** Upper Bound: 1

Command Line Argument: 

        --m 0.005

**Mutation Rate**

Default value: 0.000000012 ** Lower Bound: 0 ** Upper Bound: 1

Command Line Argument: 

        --r 0.000000012

**Lower of Ne Range**

Default value: 10 ** Lower Bound: 10 ** Upper Bound: -

Command Line Argument: 

        --lNe 10

**Upper of Ne Range**

Default value: 500 ** Lower Bound: - ** Upper Bound: 500

Command Line Argument: 

        --uNe 500

**Lower of Theta Range**

Default value: 1 ** Lower Bound: 1 ** Upper Bound: -

Command Line Argument: 

        --lT 1

**Upper of Theta Range**

Default value: 10 ** Lower Bound: - ** Upper Bound: 10

Command Line Argument: 

        --uT 10

**Number of OneSamp Trials**

Default value: 50000 ** Lower Bound: 2 ** Upper Bound: 50,000

Command Line Argument: 

        --s 50

**Lower of Duration Range**

Default value: 2 ** Lower Bound: 2 ** Upper Bound: -

Command Line Argument: 

        --lD 2

**Upper of Duration Range**

Default value: 8 ** Lower Bound: - ** Upper Bound: 8

Command Line Argument: 

        --uD 8

**Missing data for individuals**

Default value: 0.2 ** Lower Bound: 0 ** Upper Bound: 1

Command Line Argument: 

        --i 0.2

**Missing data for loci**

Default value: 0.2 ** Lower Bound: 0 ** Upper Bound: 1

Command Line Argument: 

        --l 0.2
        
**--o The File Name**

Default value: - ** Lower Bound: - ** Upper Bound: -

Command Line Argument: 

        --o exampleData/genePop5Ix5L



**OUTPUT**
##
ONeSAMP3.0 outputs the mean, median, and 95 credible limits for the posterior distribution of the effective population size using standard out. 
        


## RUN


**HOW TO EXECUTE USING COMMAND LINES**



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


        python3 ./onesamp3.py --s 2 --o exampleData/genePop5Ix5L
        
   Adjust configurations accordingly
 
 
## 


**HOW TO EXECUTE USING DOCKER**



1. Download Docker using the following link:

        https://docs.docker.com/get-docker/

2. Run the following command line:

        docker pull aaronhong10245/onesamp:latest
        
3. Only execute if your computer has an AMD chip: 
        
        export DOCKER_DEFAULT_PLATFORM=linux/amd64   
        
4. To run ONeSAMP3.0:
        
        docker run aaronhong10245/onesamp python3 ./root/OneSamp_python/onesamp3.py --s 2 --o /root/OneSamp_python/exampleData/genePop5Ix5L
    
    Adjust configurations accordingly


