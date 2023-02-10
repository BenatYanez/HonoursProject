#!/bin/sh

# Modified 10/02/2023
# 
# Modification of Matthew Fellbaum's script
# Obtain number of unique haplotypes as defined by BED file
# Need to first run 'sconda atsweeps'

# Grid Engine options (lines prefixed with #$)
#$ -N vcf_AT_hap									# Name of job in 'wstat' list
#$ -V												# Pass current environment to job
#$ -cwd												# Run file from current working directory
#$ -l h=c1											# Run array job on this sub-server
#$ -o /data/hartfield/atsweeps/scripts/output/		# Folder for STDOUT print
#$ -e /data/hartfield/atsweeps/scripts/error/		# Folder for STDERR print

# Obtain hapcount as defined by BED file
vcftools --vcf /data/hartfield/atsweeps/analyses/Vcf_subset_AT.vcf --hapcount hapcount_BED_file_20kb.bed --chr 1 --min-alleles 2 --max-alleles 2 --remove-indels --max-missing 1 --min-meanDP 17 --max-meanDP 49 --out /scratch/byanez/hapcount

# Copy file back to head folder
rsync -avz /scratch/byanez/hapcount.hapcount /data/hartfield/atsweeps/analyses/byanez

# EOF
