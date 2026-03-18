## python ltr_identity_extractor.py [-h] [-u] ltr_pos_file ref_genome

from __future__ import division
import argparse
import subprocess
import re
import sys

sys.path.append("/nfs_netapp/v1amaso4/lib/python2.7/mason-code")
import static_functions

parser = argparse.ArgumentParser(
	description="processes LTR positions, extracts sequence, aligns and calculates identity",
	epilog="Author: Andrew Mason; Release: 14/04/15; Contact: andrew.mason@roslin.ed.ac.uk")
parser.add_argument("ltr_pos_file", help="positions file of ltrs - contig, LTR1_start, LTR1_end, LTR2_start, LTR2_end, source")
parser.add_argument("ref_gen", help="reference genome file for sequence extraction")
parser.add_argument("-u", "--use_full", action="store_true", help="if select seq-extract will use full contig name for extract, rather than default not")
usr_args = parser.parse_args()

# create genome_dict for sequence extract
print("\nCreating reference genome dictionary for sequence extraction.")
ref_dict = static_functions.seq_dict_creator(usr_args.ref_gen)

lines = (subprocess.check_output("wc -l " + usr_args.ltr_pos_file + " | awk \'{print $1}\'", shell=True)).rstrip("\n")
print("Complete.\n\nIterating through " + str(lines) + " LTR pairs.\nLTR pairs completed:")
i = 1
identity_file = open("identity.tmp", "w")
for ltrs in open(usr_args.ltr_pos_file).read().rstrip("\n").split("\n"):
    # extract positions
    ltr_pos = open("element_ltrs.tmp", "w")
    ltr_pos.write(str(ltrs.split("\t")[0]) + "\t" + str(ltrs.split("\t")[1]) + "\t" + str(ltrs.split("\t")[2]) + "\t+\t" + str(ltrs.split("\t")[5]) + "\n" + str(ltrs.split("\t")[0]) + "\t" + str(ltrs.split("\t")[3]) + "\t" + str(ltrs.split("\t")[4]) + "\t+\t" + str(ltrs.split("\t")[5]))
    ltr_pos.close()
    # extract sequence
    cont_full = "N"
    if usr_args.use_full:
        cont_full = "Y"
    static_functions.seq_extract("element_ltrs.fasta", (static_functions.list_initial_formatter("element_ltrs.tmp")), ref_dict, cont_full)

    # perform alignment
    subprocess.call(("muscle -diags -quiet -clw -in element_ltrs.fasta -out element_ltrs.clw"), shell=True)

    # get identity
    identity_count = (subprocess.check_output(("grep -o \"*\" element_ltrs.clw | wc -l"), shell=True)).rstrip("\n")
    aln_lines = (subprocess.check_output(("awk \'{print $2}\' element_ltrs.clw | tail -n +4 | grep -v \"\\*\" | sed \'/^\\s*$/d\' | awk \'{if (NR%2){top=top$0} else {bottom=bottom$0}}END{print top \"\\n\" bottom}\'"), shell=True)).rstrip("\n").split("\n")
    
    start_gap = 0
    end_gap = 0
    for seq in aln_lines:
        st_m = re.search(r"^(-+)", seq)
        if st_m:
            start_gap = len(st_m.group(1))
        en_m = re.search(r"(-+)$", seq)
        if en_m:
            end_gap = len(en_m.group(1))

    aln_length = len(aln_lines[0]) - (start_gap + end_gap)
    seq_identity = round(((int(identity_count) * 100) / aln_length), 2)
    identity_file.write(str(seq_identity) + "\n")

    # print progress
    if ((i%100 == 0) or (i+1 == lines)):
        print(str(i))
    i+=1
    
identity_file.close()

subprocess.call(("paste " + usr_args.ltr_pos_file + " identity.tmp > ltr_positions_with_seq_identity.txt"), shell=True)
subprocess.call("rm element_ltrs* identity.tmp", shell=True)
print("Complete.\n\nProcess complete.\n")
    
