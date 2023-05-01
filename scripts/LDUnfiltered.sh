
#!/bin/sh
#Calculate linakge disequilibrium at sites 20kb apart for data that has been filtered
# Need to first run 'sconda atsweeps'

# Grid Engine options (lines prefixed with #$)
#$ -N vcf_LD_NoFil                                                                    # Name of job in 'wstat' list
#$ -V                                                                                           # Pass current environment to job
#$ -cwd                                                                                 # Run file from current working directory
#$ -l h=c4                                                                      # Run array job on this sub-server
#$ -pe smp 10
#$ -t 1-5                   # Run command for each chromosome
#$ -o /data/hartfield/atsweeps/scripts/output/          # Folder for STDOUT print
#$ -e /data/hartfield/atsweeps/scripts/error/           # Folder for STDERR print

#Create a variable to carry out the script for every chromosome
chrmnm=${SGE_TASK_ID}
vcftools --vcf /data/hartfield/atsweeps/analyses/Vcf_subset_AT.vcf --hap-r2 --chr $chrmnm \
--maf 0.1 --min-meanDP 17 --max-meanDP 49 --thin 20000 --out /scratch/byanez/NoFilter_LD_chromosome_"$chrmnm"

#Copy file back to head folder
rsync -avz /scratch/byanez/NoFilter_LD_chromosome_"$chrmnm".hap.ld /data/hartfield/atsweeps/analyses/byanez/linkage
