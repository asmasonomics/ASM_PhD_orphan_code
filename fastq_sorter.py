## python fastq_sorter.py read1.fastq.gz read2.fastq.gz

import argparse
import subprocess

parser = argparse.ArgumentParser(
	description="sort paired fastq files",
	epilog="Author: Andrew Mason; Release: 03/11/2016; Contact: andrew.mason@roslin.ed.ac.uk")
parser.add_argument("read1", help="read1 of paired fastq files")
parser.add_argument("read2", help="read2 of paired fastq files")
usr_args = parser.parse_args()

read1_outname = (usr_args.read1).replace("_read","_sorted_read").replace(".gz","")
read2_outname = (usr_args.read2).replace("_read","_sorted_read").replace(".gz","")

def fastq_processor(read_file, pair_str):
	
	read = []
	read_names = []
	read_names2 = []
	
	i = 0
	tmp = []
	for line in ((subprocess.check_output("zcat " + read_file, shell=True)).rstrip("\n")).split("\n"):
		tmp += [line]
		i+=1
		if i==1:
			read_names += [line.replace(pair_str,"")]
			read_names2 += [line]
		elif i==4:
			read += [tmp]
			tmp = []
			i=0
	
	return read, read_names, read_names2



read1, read1_names, r1_names = fastq_processor(usr_args.read1, "/1")
read2, read2_names, r2_names = fastq_processor(usr_args.read2, "/2")

read1_names.sort()
r1_names.sort()
r1_set = set(read1_names)
read2_names.sort()
r2_names.sort()
r2_set = set(read2_names)

unmatched = r1_set.symmetric_difference(r2_set)
matched = r1_set.intersection(r2_set)


def write_fastq(outname,read_list):
	out = open(outname, "w")
	for i in read_list:
		for j in i:
			out.write(j + "\n")
	out.close()
	subprocess.call("gzip " + outname, shell=True)
	return


if len(unmatched) == 0:
	
	r1_order = dict((x[0], tuple(x[1:])) for x in read1)
	read1[:] = [((x,) + r1_order[x]) for x in r1_names]
	write_fastq(read1_outname, read1)
	
	r2_order = dict((x[0], tuple(x[1:])) for x in read2)
	read2[:] = [((x,) + r2_order[x]) for x in r2_names]
	write_fastq(read2_outname, read2)
	
else:
	print("There are unmatched pairs:")
	print(unmatched)
	
