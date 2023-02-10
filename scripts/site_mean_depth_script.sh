#!/bin/sh

# 15th February 2021
# Matthew Fellbaum
# Modification of Gertjan Bisschop's script for running array job on server
# Obtain mean sequencing depth for all sites in chromosome 1 from the vcf subset file
# Additional filters added from previous site-mean-depth script
# Need to first run 'sconda atsweeps'

# Grid Engine options (lines prefixed with #$)
#$ -N vcf_AT                                            # Name of job in 'wstat' list
#$ -V                                                   # Pass current environm>
#$ -cwd                                                 # Run file from current>
#$ -l h=c2                                              # Run array job on this>
#$ -o /data/hartfield/atsweeps/scripts/output/          # Folder for STDOUT print
#$ -e /data/hartfield/atsweeps/scripts/error/           # Folder for STDERR print

# Obtain sequencing depth for all sites in chromosome 1
vcftools --vcf /data/hartfield/atsweeps/analyses/Vcf_subset_AT.vcf --site-mean-depth --chr 1 --min-alleles 2 --max-alle>
 
# Copy file back to head folder
rsync -avz /scratch/mfellbaum/filtered_site_mean_depth_chr1_AT.ldepth.mean  /data/hartfield/atsweeps/analyses/
  
# EOF