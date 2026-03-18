##python hexamer_model.py ref_genome_seq_only.fa

import sys
import random
import numpy
import scipy.stats

genome_str = open(sys.argv[1]).read().rstrip("\n")
print(len(genome_str))
genome_len = len(genome_str) - 6

cat0lt = []
cat1lt = []
cat2lt = []
cat3lt = []
cat4lt = []
cat5lt = []
cat6lt = []


rep = 0
while rep < 1000:

	nums = []
	x = 0
	while x < 312:
		st = random.randint(0,genome_len)
		hexamer = (genome_str[st:(st+6)]).upper().replace("G","C")
		nums += [hexamer.count("C")]
		x+=1
	
	
	cat0lt += [nums.count(0)]
	cat1lt += [nums.count(1)]
	cat2lt += [nums.count(2)]
	cat3lt += [nums.count(3)]
	cat4lt += [nums.count(4)]
	cat5lt += [nums.count(5)]
	cat6lt += [nums.count(6)]
	
	rep += 1
	

print(sum(cat0lt)/float(len(cat0lt)), numpy.std(cat0lt), scipy.stats.t.interval(0.95, len(cat0lt)-1,loc=numpy.mean(cat0lt),scale=scipy.stats.sem(cat0lt)))
print(sum(cat1lt)/float(len(cat1lt)), numpy.std(cat1lt), scipy.stats.t.interval(0.95, len(cat1lt)-1,loc=numpy.mean(cat1lt),scale=scipy.stats.sem(cat1lt)))
print(sum(cat2lt)/float(len(cat2lt)), numpy.std(cat2lt), scipy.stats.t.interval(0.95, len(cat2lt)-1,loc=numpy.mean(cat2lt),scale=scipy.stats.sem(cat2lt)))
print(sum(cat3lt)/float(len(cat3lt)), numpy.std(cat3lt), scipy.stats.t.interval(0.95, len(cat3lt)-1,loc=numpy.mean(cat3lt),scale=scipy.stats.sem(cat3lt)))
print(sum(cat4lt)/float(len(cat4lt)), numpy.std(cat4lt), scipy.stats.t.interval(0.95, len(cat4lt)-1,loc=numpy.mean(cat4lt),scale=scipy.stats.sem(cat4lt)))
print(sum(cat5lt)/float(len(cat5lt)), numpy.std(cat5lt), scipy.stats.t.interval(0.95, len(cat5lt)-1,loc=numpy.mean(cat5lt),scale=scipy.stats.sem(cat5lt)))
print(sum(cat6lt)/float(len(cat6lt)), numpy.std(cat6lt), scipy.stats.t.interval(0.95, len(cat6lt)-1,loc=numpy.mean(cat6lt),scale=scipy.stats.sem(cat6lt)))





