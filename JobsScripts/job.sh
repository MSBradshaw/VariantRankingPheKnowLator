#!/bin/bash

#SBATCH --time=24:00:00   # walltime
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --mem=8G   # memory per CPU core
#SBATCH --mail-user=michael.bradshawiii@colorado.edu   # email address
##SBATCH --mail-type=FAIL

/Users/mibr6115/python3/bin/python3 wrwr.py $1 $2 $3 $4 $5 $6 $7 $8 $9 $10
