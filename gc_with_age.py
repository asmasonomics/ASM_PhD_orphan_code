## python gc_with_age.py [-h] seq_file.fasta

from __future__ import division
import subprocess
import argparse
parser = argparse.ArgumentParser(
	description="gc_with_age takes a multi sequence fasta file with sequences of different ltr identity, splits it and reports the gc for each element with the age",
	epilog="Author: Andrew Mason; Release: 12/06/15; Contact: andrew.mason@roslin.ed.ac.uk")
parser.add_argument("seq_file", help="fastA format sequence file to be analysed")
usr_args = parser.parse_args()

import sys
sys.path.append("/nfs_netapp/v1amaso4/lib/python2.7/mason-code")
import static_functions

all_files = static_functions.fasta_splitter(open(usr_args.seq_file).read(), "Y")

outfile = open("gc_with_age.txt", "w")
dev_results = open("deviant_results.txt", "w")

for seq in all_files:
    sequence=""
    for line in (open(seq).read().rstrip("\n").split("\n")):
            if not line.startswith(">"):
                sequence+=line

    seq_one_line = ((sequence.replace("\r", "").replace("\n", "")).upper()).replace("U", "T")
    GC = seq_one_line.count("G") + seq_one_line.count("C")
    ATGC = GC + seq_one_line.count("A") + seq_one_line.count("T")
    GC_content = round(((float(100 / ATGC))*GC), 3)
    age = ((seq.split("_"))[-1]).replace(".fas", "")
    outfile.write(str(age) + "\t" + str(GC_content) + "\n")

    if (float(GC_content)>60.36) or (float(GC_content)<33.54):
        dev_results.write(str(seq.replace(".fas", "")) + "\t" + str(GC_content) + "\n")

outfile.close()
dev_results.close()
subprocess.call("rm *.fas", shell=True)
