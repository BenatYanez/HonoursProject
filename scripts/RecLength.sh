#!/usr/bin/env bash

#Beñat Yañez
#Model a neutral population of Arabidopsis thaliana
#Parameters from Southern Sweden population Huber et al
#Need to run "sconda simulations"

# Grid Engine options (lines prefixed with #$)
#$ -N RecLength                                                                      # Name of job in 'wstat' list
#$ -V                                                                                           # Pass current environm$
#$ -cwd                                                                                         # Run file from current$
#$ -l h=c3
#$ -pe smp 10                                                                                   # Run array job on this$
#$ -t 1-5
#$ -o /data/hartfield/atsweeps/scripts/output/          # Folder for STDOUT print
#$ -e /data/hartfield/atsweeps/scripts/error/           # Folder for STDERR print

chrmnm=${SGE_TASK_ID}
#Define the number of replicates
replicate=150

python RecLength.py "$chrmnm" "$replicate" > LegthofRecombinationEvent_Chrom"$chrmnm".txt

rsync -avz /scratch/byanez/LengthofRecombinationEvent_Chrom"$chrmnm".txt /data/hartfield/atsweeps/analyses/byanez
rsync -avz /scratch/byanez/Chr"$chrmnm"_RecombinationLength.txt /data/hartfield/atsweeps/analyses/byanez
