#!/bin/bash
#SBATCH --job-name=oneSamp    # Job name
#SBATCH --mail-type=ALL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=yu.hong@ufl.edu     # Where to send mail
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=1gb                     # Job memory request
#SBATCH --time=24:00:00               # Time limit hrs:min:sec
#SBATCH --output=serial_test_%j.log   # Standard output and error log
pwd; hostname; date

module load R/4.1
chmod +rwx /blue/boucher/yu.hong/OneSampTest/build/OneSamp

echo "Running plot script on a single CPU core"

python /blue/boucher/yu.hong/oneSampTest/main.py --s 20000 --o /blue/boucher/yu.hong/oneSampTest/DataSet/genePop100Ix100L > /blue/boucher/yu.hong/oneSampTest/genePop100Ix100L.out

date
