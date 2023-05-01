#Run a simulation of ancestry and mutations based on seocndaryContact6 model from Huber et la., 2015 and obtain a vcf file for each repeat
import sys
import msprime
import tskit
import os.path
import numpy as np

demography = msprime.Demography()

#Add South Sweden
demography.add_population(name="South", initial_size=205840)

#Add North Sweden
demography.add_population(name="North", initial_size=21080)

#Add Ancestral Population
demography.add_population(name="Ancestral", initial_size=124000)

#Set initial migration rate, might need to change "South" to 0 and North to 1, transformed migration rate from 2N0 to m use 2N0 because the haploid individuals

demography.set_migration_rate(source="South", dest="North", rate=2.38e-5)

demography.set_migration_rate(source="North", dest="South", rate=4.96e-6)

#Change migration rate to zero at time 0.08- Number of generations should be 39000
demography.add_symmetric_migration_rate_change(time=39000, populations=["South","North"], rate=0)
#Merge the two populations into one ancestral one at time 0.31, Number of generations should be 153000
demography.add_population_split(time=153000, derived=["South","North"], ancestral="Ancestral")

#Translate the two variables specified outside the script into the code
chromosome = int(sys.argv[1])
repeat = int(sys.argv[2])

#Create a recombination rate map from the adjusted recombination rates
recomb = open (f'/data/hartfield/atsweeps/analyses/byanez/Recombination/Salome_Chr{chromosome}Map_Adjusted_Space.dat',)
rate_map = msprime.RateMap.read_hapmap(recomb, map_col=2)
recomb.close()
#Simulate the ancestry, haploid genome since it is highly selfing which in simulations is similar to outcrossing diploid genomes
ts = msprime.sim_ancestry({"South":45},ploidy=1, recombination_rate=rate_map, demography=demography)
ts = msprime.sim_mutations (ts, rate=7e-9,model=msprime.BinaryMutationModel())

#Save file
directory= '/scratch/byanez/'
vcf_file_name = f'Vcf_Chr{chromosome}_Neutral_repeat{repeat}_Real'
complete_file = os.path.join(directory,vcf_file_name)
with open(f'{complete_file}.vcf', "w") as vcf_files:
  ts.write_vcf(vcf_files)
