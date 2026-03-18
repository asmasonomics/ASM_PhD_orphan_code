## python model_sampling.py [-h] insert_freq sample_size pop_size

from __future__ import division
import argparse
import math
from random import sample
import numpy

parser = argparse.ArgumentParser(
        description="Model for individual sequencing data",
        epilog="Author: Andrew Mason; Release: 06/11/2016; Contact: andrew.mason@roslin.ed.ac.uk")
parser.add_argument("insert_freq", help="decimal frequency for the ALVE insert")
parser.add_argument("sample_size", help="number of unique individuals to sample from the population")
parser.add_argument("pop_size", help="total population size available for sampling")
usr_args = parser.parse_args()

## S1 - define genotypes and pop structure
# use HW to get individuals of genotypes given the insert frequency 
hom_insert = float(usr_args.insert_freq) * float(usr_args.insert_freq)
hom_indiv = int(round(hom_insert * int(usr_args.pop_size)))
het_insert = 2 * (1 - float(usr_args.insert_freq)) * float(usr_args.insert_freq)
het_indiv = int(round(het_insert * int(usr_args.pop_size)))
wt_indiv = int(usr_args.pop_size) - (hom_indiv + het_indiv)

# create list of indiv genotypes, where 0=wt, 1=het, 2=hom_insert
pop_genotypes = []
pop_genotypes += [0] * wt_indiv
pop_genotypes += [1] * het_indiv
pop_genotypes += [2] * hom_indiv

#print(pop_genotypes)


## S2 - iterate through a sampling procedure
modeled_genotypes = []
i=0
#iterations = (math.factorial(int(usr_args.pop_size)))/((math.factorial(int(usr_args.sample_size)))*(math.factorial(int(usr_args.pop_size)-int(usr_args.sample_size))))
#print(iterations)
iterations = 1000000
while i<iterations:
	subset = sample(pop_genotypes, int(usr_args.sample_size))
	mod_freq = round((float(sum(subset) / (int(usr_args.sample_size)*2))), 3)
	modeled_genotypes += [mod_freq]
	i+=1

modeled_genotypes.sort()
#print(modeled_genotypes)
print("Min freq:\t" + str(modeled_genotypes[0]))
print("Max freq:\t" + str(modeled_genotypes[-1]))
print("Average freq:\t" + str(round((float(sum(modeled_genotypes)/len(modeled_genotypes))),4)))
print("\n% detecting:\t" + str(100-(((modeled_genotypes.count(0.0))/iterations)*100)))




