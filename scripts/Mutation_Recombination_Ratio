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
recomb.close()
#Create a function to simulate the ancestry as a replicate
def sim_replicates(repeat):
 ancestry_reps= msprime.sim_ancestry({"South":45},ploidy=1, recombination_rate=rate_map, demography=demography, num_replicates=repeat)
 for ts in ancestry_reps:
  mutated_ts = msprime.sim_mutations(ts, rate=7e-9, model=msprime.BinaryMutationModel())
  yield mutated_ts
np.set_printoptions(suppress=True,formatter={'float_kind':'{:f}'.format})

#Obtain the number of mutations and trees in each replicate

S=np.zeros(repeat)
X=np.zeros(repeat)
Y=np.zeros(repeat)
for replicate_index, ts in enumerate(sim_replicates(repeat)):
 X[replicate_index] =ts.num_trees
 S[replicate_index] =ts.num_mutations
 Y[replicate_index]= ts.num_trees / ts.num_mutations
print("Number of Trees(Recombination Events)")
print(X)
print("Number of Mutations")
print(S)
print("Mean number of Mutations:"), print(np.mean(S))
print("Standard Deviation of Mutations:"), print(np.std(S))
print("Mean number of Trees(Recombiantion Events)"),print (np.mean(X))
print("Standard Deviation of Trees"), print (np.std(X))
print("Ratio of Recombination to Mutation")
print(Y)
print("Mean Ratio"), print (np.mean(Y)), print ("Standard Deviation:"), print(np.std(Y))
file_name = f'Chr{chromosome}_MutationRecombinationRatio'
np.savetxt(f'{file_name}.txt',Y,delimiter=",")
