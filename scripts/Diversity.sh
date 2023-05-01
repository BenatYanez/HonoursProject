#!/bin/sh

# Modified 20/02/2023
#
# Modification of Matthew Fellbaum's script
# Obtain the diversity between two sites at a given distance
# Need to first run 'sconda atsweeps'

# Grid Engine options (lines prefixed with #$)
#$ -N vcf_diveristy                                                                     # Name of job in 'wstat' list
#$ -V                                                                                           # Pass current environm$#$ -cwd                                                                                 # Run file from current working$#$ -l h=c6                                                                      # Run array job on this sub-server
#$ -cwd
#$ -l h=c6
#$ -pe smp 10
#$ -t 1-5                   # Run command for each chromosome
#$ -o /data/hartfield/atsweeps/scripts/output/          # Folder for STDOUT print
#$ -e /data/hartfield/atsweeps/scripts/error/           # Folder for STDERR print
#Create a variable to carry out the script for every chromosome
chrmnm=${SGE_TASK_ID}

#Define size of window
size=100000
# Obtain the diveristy at a given window
vcftools --vcf /data/hartfield/atsweeps/analyses/Vcf_subset_AT.vcf --window-pi "$size" --chr "$chrmnm" \
 --out /scratch/byanez/diversity_windowof_"$size"_chrm"$chrmnm"_Data
# Copy file back to head folder
rsync -avz /scratch/byanez/diversity_windowof_"$size"_chrm"$chrmnm"_Data.windowed.pi /data/hartfield/atsweeps/analyses/byanez/diversity
