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

demography.set_migration_rate(source="North", dest="South", rate=2.38e-5)

demography.set_migration_rate(source="South", dest="North", rate=4.96e-6)

#Change migration rate to zero at time 0.08- Number of generations should be 39000
demography.add_symmetric_migration_rate_change(time=39000, populations=["South","North"], rate=0)
#Merge the two populations into one ancestral one at time 0.31, Number of generations should be 153000
demography.add_population_split(time=153000, derived=["South","North"], ancestral="Ancestral")
#print (demography.debug())

chromosome = int(sys.argv[1])
repeat = int(sys.argv[2])
#Create a recombination rate map, maybe mat need to half rates becuase of a haploid genome
recomb = open (f'/data/hartfield/atsweeps/analyses/byanez/Recombination/Salome_Chr{chromosome}Map_Adjusted_Space.dat',)
rate_map = msprime.RateMap.read_hapmap(recomb, map_col=2)
#print(rate_map)
recomb.close()
#Run the ancestry simulation for two individuals during a number of replicates
replicates =msprime.sim_ancestry({"South":2},ploidy=1, recombination_rate=rate_map, demography=demography, num_replicates=repeat)
X=np.zeros(repeat)
np.set_printoptions(threshold=sys.maxsize)
for replicate_index, ts in enumerate(replicates):
 breakpoints= ts.breakpoints(as_array=True)
 differences= np.diff(breakpoints)
 avg_diff= np.mean(differences)
 #print(avg_diff)
 X[replicate_index] =avg_diff
print("Distance between recombination(bp)(Average of Chromosome):")
print(X)
print("Average distance between recombinations (bp):")
print(np.mean(X))

file_name = f'Chr{chromosome}_RecombinationLength'
np.savetxt(f'{file_name}.txt',X,delimiter=",")
