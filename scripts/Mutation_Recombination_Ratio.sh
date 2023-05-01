#!/usr/bin/env bash

#Beñat Yañez
#Model a neutral population of Arabidopsis thaliana
#Parameters from Southern Sweden population Huber et al
#Need to run "sconda simulations"

# Grid Engine options (lines prefixed with #$)
#$ -N Mut_RecNum                                                                      # Name of job in 'wstat' list
#$ -V                                                                                           # Pass current environment to job
#$ -cwd                                                                                         # Run file from current working directory
#$ -l h=c3
#$ -pe smp 10                                                                                   # Run array job on this sub-server
#$ -t 1-5
#$ -o /data/hartfield/atsweeps/scripts/output/          # Folder for STDOUT print
#$ -e /data/hartfield/atsweeps/scripts/error/           # Folder for STDERR print

chrmnm=${SGE_TASK_ID}
#Define the number of replicates
replicate=150

python Mut_RecNumber.py "$chrmnm" "$replicate" > Mutation_RecombinationCount_Chrom"$chrmnm".txt

rsync -avz /scratch/byanez/Mutation_RecombinationCount_Chrom"$chrmnm".txt /data/hartfield/atsweeps/analyses/byanez
rsync -avz /scratch/byanez/Chr"$chrmnm"_MutationRecombinationRatio.txt /data/hartfield/atsweeps/analyses/byanez
