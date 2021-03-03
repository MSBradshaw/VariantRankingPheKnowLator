#!/bin/bash

#SBATCH --time=24:00:00   # walltime
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --mem=8G   # memory per CPU core
#SBATCH --mail-user=michael.bradshawiii@colorado.edu   # email address
##SBATCH --mail-type=FAIL

/Users/mibr6115/python3/bin/python3 anti_degree_wrwr.py --edge_list Edgelists/string_hpo_CFTR_direct_HPO_connections_removed_triple_list.txt --disease_sets cftr.txt --threshold 0.000011

