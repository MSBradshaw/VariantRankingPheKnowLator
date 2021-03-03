#!/bin/bash

#SBATCH --time=167:59:59   # walltime
#SBATCH --ntasks=64   # number of processor cores (i.e. tasks)
#SBATCH --mem=8G   # memory per CPU core
#SBATCH --mail-user=michael.bradshawiii@colorado.edu   # email address
##SBATCH --mail-type=FAIL
#SBATCH -p long

/Users/mibr6115/python3/bin/python3 $1 $2 $3 $4 $5 $6 $7 $8 $9 $10
