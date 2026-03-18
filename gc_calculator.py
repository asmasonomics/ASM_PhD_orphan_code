## python gc_calculator.py [-h] seq_file.fasta

from __future__ import division
import numpy
import argparse
parser = argparse.ArgumentParser(
	description="gc_calculator reports the total number of sequences and total GC content of a single or multi fastA format sequence file.",
	epilog="Author: Andrew Mason; Release: 25/11/14; Contact: andrew.mason@roslin.ed.ac.uk")
parser.add_argument("seq_file", help="fastA format sequence file to be analysed")
usr_args = parser.parse_args()

seq_num = 0
seq = ""
for line in (open(usr_args.seq_file).read().rstrip("\n").split("\n")):
    if line.startswith(">"):
        seq_num += 1
        seq += "$"
    else:
        seq += line

print("\nPreparing sequence file: " + usr_args.seq_file)
print("All non \"AaTtGgCc\" characters will be ignored.\n")

total_seq = ((seq.replace("\r", "").replace("\n", "")).upper()).replace("U", "T")
indiv_seq = total_seq.split("$")
indiv_seq.pop(0)

gc_values = 0
GC_val = ""
n_values = 0
non_natgc = 0
total_len = 0
for seq in indiv_seq:
    seq_len = len(seq)
    n_count = seq.count("N")
    a_count = seq.count("A")
    t_count = seq.count("T")
    g_count = seq.count("G")
    c_count = seq.count("C")
    GC = g_count + c_count
    atgc_count = a_count + t_count + g_count + c_count

    gc_values += round(((float(100 / atgc_count)) * GC), 2)
    GC_val += str(round(((float(100 / atgc_count)) * GC), 2)) + "$"
    n_values += n_count
    non_natgc += (seq_len - (atgc_count + n_count))
    total_len += seq_len

avg_gc = str(round((gc_values / seq_num), 2)) + "%"
sd = round((numpy.std(map(float, (GC_val.rstrip("$")).split("$")))), 2)

print("Total sequences in file           : " + str(seq_num))
print("Total sequence length             : " + str(round(float(total_len / 1000000), 4)) + " Mb\n")
print("Total sequence characters in file : " + str(total_len))
print("ATGC characters analysed          : " + str(total_len - (n_values + non_natgc)) + "\n")
print("Percentage of characters analysed : " + str(round((float(100/total_len) * (total_len - (n_values + non_natgc))), 2)) + "%")
print("Total number of N                 : " + str(n_values) + " (" + str(round((float(100/total_len) * n_values), 2)) + "%)")
print("Total number of non ATGCN         : " + str(non_natgc) + " (" + str(round((float(100/total_len) * non_natgc), 2)) + "%)")
print("\nAverage sequence GC content       : " + avg_gc + " (SD: " + str(sd) + ")\n")
