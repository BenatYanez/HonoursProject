#!/usr/bin/env bash

# Modified 20/02/2023
#
# Modification of Matthew Fellbaum's script
# Obtain number of unique haplotypes as defined by BED file for the simulations
# Need to first run 'sconda atsweeps'

# Grid Engine options (lines prefixed with #$)
#$ -N vcf_hap_count_sim                                                                     # Name of job in 'wstat' list
#$ -V                                                                                           # Pass current environment to job
#$ -cwd                                                                                 # Run file from current working directory
#$ -l h=c6                                                                      # Run array job on this sub-server
#$ -pe smp 5
#$ -t 1-5                   # Run command for each chromosome
#$ -o /data/hartfield/atsweeps/scripts/output/          # Folder for STDOUT print
#$ -e /data/hartfield/atsweeps/scripts/error/           # Folder for STDERR print

#Create a variable to carry out the script for every chromosome
chrmnm=${SGE_TASK_ID}

#Define size of BED file in kb
size=10
#Define the number of replicates
replicate=150
for (( n=1; n <= $replicate; ++n))
 do
# Obtain hapcount as defined by BED file
 vcftools --vcf /data/hartfield/atsweeps/analyses/byanez/VCF/Vcf_Chr"$chrmnm"_Neutral_Diploid_repeat"$n"_Real.vcf --hapcount /data/hartfield/atsweeps/analyses/byanez/bedfiles/hapcount_BED_file_chrom"$chrmnm"_"$size"kb_neutral.bed \
         --out /scratch/byanez/hapcount_"$size"kb_chromosome_"$chrmnm"_neutralsim_repeat"$n"_Real

# Copy file back to head folder
 rsync -avz /scratch/byanez/hapcount_"$size"kb_chromosome_"$chrmnm"_neutralsim_repeat"$n"_Real.hapcount /data/hartfield/atsweeps/analyses/byanez/hapcount
done
# EOF
