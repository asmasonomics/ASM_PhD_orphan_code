## python pseudo_fastq.py read_cat.fastq.gz prefix

import argparse
import subprocess

parser = argparse.ArgumentParser(
	description="create pseudo single end fastq from concatenated paired data",
	epilog="Author: Andrew Mason; Release: 16/11/2016; Contact: andrew.mason@roslin.ed.ac.uk")
parser.add_argument("read_cat", help="concatenated fastq file")
parser.add_argument("prefix", help="prefix for each fastq line (which will be followed by a number), and prefix for out fastq")
usr_args = parser.parse_args()

read = []
i = 0
j = 0
tmp = []
for line in ((subprocess.check_output("zcat " + usr_args.read_cat, shell=True)).rstrip("\n")).split("\n"):
	tmp += [line]
	i+=1
	if i==4:
		read += [tmp]
		j+=1
		tmp = []
		i=0

zero_fill = len(str(j))
out = open(usr_args.prefix + ".fastq", "w")
y=1
for x in read:
	out.write("@" + usr_args.prefix + "-" + str(y).zfill(zero_fill) + "\n" + x[1] + "\n" + x[2] + "\n" + x[3] + "\n")
	y+=1
out.close()
subprocess.call("gzip " + usr_args.prefix + ".fastq", shell=True)
