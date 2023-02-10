#!/bin/sh
#Modified 10/02/2023
# Modification of Matthew Fellbaum script for obtaining sequencing depth in chromosome 1
# Obtain mean sequencing depth for all sites in the genome from the vcf subset file
# Need to first run 'sconda atsweeps'

# Grid Engine options (lines prefixed with #$)
 -N vcf_AT                                            # Name of job in 'wstat' list
 -V                                                   # Pass the variables in the current environment to the job
 -cwd                                                 # Run file from current directory
 -l h=c1                                              # Run array job on this node
 -o /data/hartfield/atsweeps/scripts/output/          # Folder for STDOUT print
 -e /data/hartfield/atsweeps/scripts/error/           # Folder for STDERR print

# Obtain sequencing depth for all sites in the genome
vcftools --vcf /data/hartfield/atsweeps/analyses/Vcf_subset_AT.vcf --site-mean-depth --min-alleles 2 --max-alleles 2 --remove-indels --max-missing 1
--out /scratch/byanez/filtered_genome_site_mean_depth

 
# Copy file back to head folder
rsync -avz /scratch/byanez/filtered_site_mean_depth.ldepth.mean  /data/hartfield/atsweeps/analyses/byanez/sitedepth
  
# EOF
