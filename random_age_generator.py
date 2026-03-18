## python random_age_generator.py [-h] [--iterations ITERATIONS] ori_pos_file

import argparse
import subprocess
import sys
sys.path.append("/nfs_netapp/v1amaso4/lib/python2.7/mason-code")
import static_functions
from random import shuffle

parser = argparse.ArgumentParser(
	description="random_age_generator takes a positions file with LTR identities and redistributes the identities at random.",
	epilog="Author: Andrew Mason; Release: 11/05/15; Contact: andrew.mason@roslin.ed.ac.uk")
parser.add_argument("ori_pos_file", help="original positions file, tab delimited, column order: chr, start, end, strand, source, identity")
parser.add_argument("--iterations", type=int, default=1000, help="set number of random age files to produce (default=1000)")
usr_args = parser.parse_args()

# extract the positional information and store in temp file
subprocess.call(("cut -f1-5 " + usr_args.ori_pos_file + " > pos.tmp"), shell=True)
# extract identities and store in temp file
subprocess.call(("cut -f6 " + usr_args.ori_pos_file + " > ages.tmp"), shell=True)

# format positional information and store as list of lists
element_pos = static_functions.list_initial_formatter("pos.tmp")
# store ages as list
age_list = open("ages.tmp").read().rstrip("\n").split("\n")
# define number of leading zeros required
fill = len(str(usr_args.iterations))

# for each iteration: create outfile, shuffle age list, write pos with new age, remove blank lines
i = 1
while i <= (usr_args.iterations):
	out_file_name = "random_age_" + str(i).zfill(fill) + ".txt"
	out_file = open(out_file_name, "w")
	shuffle(age_list)
	z=0
	for pos in element_pos:
		out_file.write(str(pos[0]) + "\t" + str(pos[1]) + "\t" + str(pos[2]) + "\t" + str(pos[3]) + "\t" + str(pos[4]) + "\t" + str(age_list[z]) + "\n")
		z+=1
	out_file.close()
	subprocess.call(("sed -i \'/^$/d\' " + out_file_name), shell=True)
	i+=1
	
# remove tmp files
subprocess.call("rm *.tmp", shell=True)
