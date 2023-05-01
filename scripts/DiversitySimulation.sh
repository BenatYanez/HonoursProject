#!/usr/bin/env bash

# Modified 20/02/2023
#
# Modification of Matthew Fellbaum's script
# Obtain the diversity of the simulation between sites of defined site
# Need to first run 'sconda atsweeps'

# Grid Engine options (lines prefixed with #$)
#$ -N vcf_diveristy_sim                                                                     # Name of job in 'wstat' list
#$ -V                                                                                           # Pass current environm$#$ -cwd                                                                                 # $#$ -cwd
#$ -l h=c6
#$ -t 1-5                   # Run command for each chromosome
#$ -o /data/hartfield/atsweeps/scripts/output/          # Folder for STDOUT print
#$ -e /data/hartfield/atsweeps/scripts/error/           # Folder for STDERR print
#Create a variable to carry out the script for every chromosome
chrmnm=${SGE_TASK_ID}
#Define size of window
size=100000
#Define the number of replicates
replicate=150
#Change number in {} to change the number of simulations
for (( n=1; n <= $replicate; ++n))
 do
# Obtain hapcount as defined by BED file
 vcftools --vcf /data/hartfield/atsweeps/analyses/byanez/VCF/Vcf_Chr"$chrmnm"_Neutral_Diploid_repeat"$n"_Real.vcf --window-pi "$size"  --out /scratch/byanez/diversity_windowof_"$size"_chrm"$chrmnm"_Sim"$n"

# Copy file back to head folder
 rsync -avz /scratch/byanez/diversity_windowof_"$size"_chrm"$chrmnm"_Sim"$n".windowed.pi /data/hartfield/atsweeps/analyses/byanez/diversity
done
