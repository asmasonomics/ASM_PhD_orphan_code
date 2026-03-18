## python random_insert_generator.py [-h] [-z] [--iterations ITERATIONS] [--replicates REPLICATES] ref_seq_file

import argparse
from random import randrange
import subprocess

parser = argparse.ArgumentParser(
	description="random_instert_generator creates .txt positions file of randomly generated inserts with replicates.",
	epilog="Author: Andrew Mason; Release: 30/03/15; Contact: andrew.mason@roslin.ed.ac.uk")
parser.add_argument("seq_file", help="reference genome fasta file")
parser.add_argument("--iterations", type=int, default=10000000, help="(int) number of iterations to perform (default = 10 million)")
parser.add_argument("--replicates", type=int, default=1000, help="(int) number of replicates to perform (default = 1 thousand)")
parser.add_argument("-z", "--zip", action="store_true", help="select to zip each output positions file")
usr_args = parser.parse_args()

print("\nGenerating " + str(usr_args.iterations) + " random insertions from the \"" + ((usr_args.seq_file).split("/"))[-1] + "\" reference fasta.")
print("Process will be replicated " + str(usr_args.replicates) + " times.\n")

# prepare sequence and contig info from the reference sequence file
sequence = ""
chromo_info = []
genome = open(usr_args.seq_file).read().rstrip("\n").split("\n")
tick = 0
i=0
for line in genome:
        if not line.startswith(">"):
                sequence += line.replace("\r","")
		tick += 1
        else:
                # this section clearly requires a certain header format - adjust if required (based on Ensembl headers)
                chromo = ((line.split(" "))[0]).replace(">", "")
                length = len(genome[tick].replace("\r",""))
                start = i+1
                i+=length
                end = i
                chromo_info += [[chromo,start,end]]
		tick+=1
#seq_one_line = sequence.replace("\r", "").replace("\n", "")

# perform random number generation individually for the number of defined replicates

zero_pad = len(str(usr_args.replicates))           # define how many leading zeros should be used for file padding
#zero_pad = 4
z=1
print("Initiating random number generation.")
while z <= usr_args.replicates:
#while z < 1001:
        # perform the random number generation, add to list and then sort list
        i=0
        insert_list = []
        while i < int(usr_args.iterations):
                insert_list += [randrange(len(sequence))]
                i+=1
        insert_list.sort(key=int)

        # create appropriate output file format with correct coordinates
        file_name = "random_inserts_" + str(str(z).zfill(zero_pad)) + ".txt"
        out_file = open(file_name, "w")
        for insert in insert_list:
                for chromo in chromo_info:
                        if (insert <= chromo[2]):
                                start = int(insert + 1) - int(chromo[1])
                                out_file.write(str(chromo[0]) + "\t" + str(start) + "\t" + str(start + 1) + "\tm\tRG\n")
                                break
        out_file.close()
        subprocess.call(("sort -k1,1 -k2,2n " + file_name + " -o " + file_name), shell=True)
        
	if usr_args.zip:
		subprocess.call(("zip -9 random_inserts_" + str(str(z).zfill(zero_pad)) + " " + file_name + " 2>&1 > /dev/null"), shell=True)
        	subprocess.call(("rm " + file_name), shell=True)

        if (z%500 == 0):
                print(str(format(z, ",d")) + " random numbers sets have been generated of " + str(format(usr_args.replicates, ",d")))
        z+=1

print("\nProcess complete.\n")
