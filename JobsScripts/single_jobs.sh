#!/bin/bash

#SBATCH --time=24:00:00   # walltime
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --mem=8G   # memory per CPU core
#SBATCH --mail-user=michael.bradshawiii@colorado.edu   # email address
##SBATCH --mail-type=FAIL

/Users/mibr6115/python3/bin/python3 wrwr.py --edge_list Edgelists/PKT_Master_Edge_List_NoOntologyData.txt --disease_sets 2_disease.txt --node_mapping Edgelists/PKT_Master_Edge_List_NoOntologyData_LABELS.txt --threshold 0.000010

/Users/mibr6115/python3/bin/python3 wrwr.py --edge_list Edgelists/PKT_Master_Edge_List_NoOntologyData.txt --disease_sets 2_disease_shuffled.txt --node_mapping Edgelists/PKT_Master_Edge_List_NoOntologyData_LABELS.txt --threshold 0.0000100
