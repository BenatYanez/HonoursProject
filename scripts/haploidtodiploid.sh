#!/bin/sh

# Modified 28/03/2023
#
# Need to first run 'sconda simulations'

#Converts haploid simulated data into diploid so that it can be fed to haplotype_count script

# Grid Engine options (lines prefixed with #$)
#$ -N vcf_hap_to_diploid                                                                     # Name of job in 'wstat' list
#$ -V                                                                                           # Pass current environment to job
#$ -cwd                                                                                 # Run file from current working directory
#$ -l h=c1                                                                      # Run array job on this sub-server
#$ -t 1-5                   # Run command for each chromosome
#$ -o /data/hartfield/atsweeps/scripts/output/          # Folder for STDOUT print
#$ -e /data/hartfield/atsweeps/scripts/error/           # Folder for STDERR print

#Create a variable to carry out the script for every chromosome
chrmnm=${SGE_TASK_ID}
#Conversion
bash haploidconversion.py "$chrmnm"

# Copy file back to head folder
rsync -avz /scratch/byanez/Vcf_Chr"$chrmnm"_Neutral_Diploid.vcf /data/hartfield/atsweeps/analyses/byanez/VCF


