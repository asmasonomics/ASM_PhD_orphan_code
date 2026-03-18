## python random_insert_generator_per_chromosome.py [-h] [--chromo CHROMO] [--replicates REPLICATES] [--iterations ITERATIONS] ref_seq_file

from __future__ import division
import argparse
import os
from random import randrange
import subprocess
import sys
sys.path.append("/nfs_netapp/v1amaso4/lib/python2.7/mason-code")
import static_functions

parser = argparse.ArgumentParser(
	description="random_insert_generator_per_chromosome creates a .txt positions file of randomly generated inserts from desired chromosomes.",
	epilog="Author: Andrew Mason; Release: 08/12/14; Contact: andrew.mason@roslin.ed.ac.uk")
parser.add_argument("seq_file", help="reference genome fasta file")
parser.add_argument("--chromo", help="specify name of chromosome to generate positions for.")
parser.add_argument("--iterations", help="specify number of iterations (usually equal to annotated positions on the same chromo")
parser.add_argument("--replicates", help="specify number of repeats for the simulation")
usr_args = parser.parse_args()

print("\nSplitting sequence file into individual query files and removing the contigs to leave sex chromosomes and assembled chromosomes.")

seq_file = open(usr_args.seq_file).read()
indiv_seq_files = static_functions.fasta_splitter(seq_file, "Y")
seq_files = list(indiv_seq_files)
i=0
for filename in indiv_seq_files:
		if usr_args.chromo:
			if not (filename.startswith(str(usr_args.chromo) + "_")):
				seq_files.pop(i)
				os.remove(filename)
			else:
				i+=1
		else:
			if not (filename[0].isdigit() or filename.startswith(("W", "Z"))):
				seq_files.pop(i)
				os.remove(filename)
			else:
				i+=1
print("Files split. Proceeding with random insert generation.\n")

for filename in seq_files:
		sequence = ""
		file_name = (filename.split("_"))[0]
		#print("\nProcessing chromosome " + file_name)
		for line in open(filename).read().rstrip("\n").split("\n"):
				if not line.startswith(">"):
						sequence += line
		seq_one_line = sequence.replace("\r", "").replace("\n", "")
		seq_len = len(seq_one_line)
		iterations = 10000
		base_steps = int(round((seq_len / 100), 0))
		if (base_steps > 10000):
			iterations = base_steps
		if usr_args.iterations:
                        iterations = int(usr_args.iterations)    

                reps = 1000
                if usr_args.replicates:
                        reps = int(usr_args.replicates)
                
		z = 1
		zero_pad = len(str(reps))
		while z <= reps:
			file_num = str(z).zfill(zero_pad)
		
			i=0
			insert_list = []
			while i < iterations:
					insert_list += [randrange(seq_len)]
					i+=1
			insert_list.sort(key=int)

			out_file = open("random_inserts_" + file_num + ".txt", "w")
			for insert in insert_list:
					out_file.write(str(file_name) + "\t" + str(insert) + "\t" + str(insert + 1) + "\tm\tRG\n")
			out_file.close()
                        if iterations > 200:
                                subprocess.call("sort -n -k1,1 -k2,2 random_inserts_" + file_num + ".txt -o random_inserts_" + file_num + ".txt", shell=True)
                                subprocess.call("zip -9 random_inserts_" + file_num + " random_inserts_" + file_num + ".txt 2>&1 > /dev/null", shell=True)
                                subprocess.call("rm random_inserts_" + file_num + ".txt", shell=True)

                        if (z%5000 == 0):
                                print(str(format(z, ",d")) + " random numbers sets have been generated of " + str(format(reps, ",d")))
			z+=1
		
		os.remove(filename)
		print("Processed chromosome " + file_name)

print("\nProcess complete.\n")
