#!/bin/sh

# Modified 21/02/2023
#
# Modification of Matthew Fellbaum's script
# Obtain number of unique haplotypes as defined by BED files for each chromosome
# Need to first run 'sconda atsweeps'

# Grid Engine options (lines prefixed with #$)
#$ -N vcf_hap_count                                                                     # Name of job in 'wstat' list
#$ -V                                                                                           # Pass current environment to job
#$ -cwd                                                                                 # Run file from current working directory
#$ -l h=c1                                                                      # Run array job on this sub-server
#$ -t 1-5                   # Run command for each chromosome
#$ -o /data/hartfield/atsweeps/scripts/output/          # Folder for STDOUT print
#$ -e /data/hartfield/atsweeps/scripts/error/           # Folder for STDERR print

#Create a variable to carry out the script for every chromosome
chrmnm=${SGE_TASK_ID}
# Obtain hapcount as defined by BED file
vcftools --vcf /data/hartfield/atsweeps/analyses/Vcf_subset_AT.vcf --hapcount /data/hartfield/atsweeps/analyses/byanez/bedfiles/hapcount_BED_file_chrom"$chrmnm"_20kb.bed --chr $chrmnm \
--min-alleles 2 --max-alleles 2 --remove-indels --max-missing 1 --min-meanDP 17 --max-meanDP 49 --out /scratch/byanez/hapcount_chromosome_"$chrmnm"

# Copy file back to head folder
rsync -avz /scratch/byanez/hapcount_chromosome_"$chrmnm".hapcount /data/hartfield/atsweeps/analyses/byanez/hapcount

# EOF
