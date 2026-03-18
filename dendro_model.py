## python dendro_model.py

import sys
import subprocess
import random

# define number of occurences of ALVEs in each line (columns)
occur = [2, 3, 3, 6, 8, 4, 2, 11, 3, 5, 2, 8, 8, 5, 5, 4, 5, 6, 4, 6, 6, 4, 4, 1, 2, 2, 5, 5, 3, 0, 11, 11, 2, 3, 10, 12, 5, 3, 2, 2, 5, 5, 1, 4, 2, 4, 2, 7, 3, 4]
# total different ALVES
alves = 66
# number of datasets used
lines = len(occur)

# open output file
out = open("modelled_matrix.txt", "w")

# perform 10 iterations of random matrix construction
z = 0
while z < 10:

	model = []
	i = 0
	# for each column (dataset) produce a random order of presence/absence based on the observed occurences
	while i < lines:
		col = ([0] * (alves - occur[i])) + ([1] * occur[i])
		random.shuffle(col)
		model += [col]
		i+=1
	
	# turn the created list of columns into actual columns
	model_formatted = zip(*model)
	
	# format the data into a string ready for MATLAB model
	data_format = ""
	for j in model_formatted:
		data_format += str(j).replace("(","").replace(")",";").replace(" ","")
	df = data_format.rstrip(";")
	
	# print to file
	out.write("data" + str(z+1) + " = [" + df + "];\n")

	z+=1

out.close()
