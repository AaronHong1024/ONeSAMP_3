#!/bin/bash
#SBATCH --job-name=oneSamp    # Job name
#SBATCH --mail-type=ALL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ishayooseph@ufl.edu     # Where to send mail
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=100gb                     # Job memory request
#SBATCH --time=168:00:00               # Time limit hrs:min:sec
#SBATCH --output=serial_test_%j.log   # Standard output and error log
pwd; hostname; date

module load python3
module load R
chmod 777 /blue/boucher/ishayooseph/oneSampTest/build/OneSamp

echo "Running plot script on a single CPU core"

#Running w default duration
python3 /blue/boucher/ishayooseph/oneSampTest/main.py --o /blue/boucher/ishayooseph/genePop500Ix500L > defaultNe.txt

#Running w 
python3 /blue/boucher/ishayooseph/oneSampTest/main.py --lD 2 --uD 8 --o /blue/boucher/ishayooseph/genePop500Ix500L > 4-10-Duration.txt


