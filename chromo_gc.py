## python chromo_gc.py ref_genome_single_line.fa alve_pos.txt
#alve_pos.txt is a file with the chromosome name and hexamer location

from __future__ import division
import sys
import subprocess

chromo = ["NC_006088.4", "NC_006089.4", "NC_006090.4", "NC_006091.4", "NC_006092.4", "NC_006093.4", "NC_006094.4", "NC_006095.4", "NC_006096.4", "NC_006097.4", "NC_006098.4", "NC_006099.4", "NC_006100.4", "NC_006101.4", "NC_006104.4", "NC_006105.4", "NC_006107.4", "NC_006110.4", "NC_006115.4", "NC_006127.4", "NC_008465.3"]

chromo_seq = []

for i in chromo:
	seq = ((subprocess.check_output(("grep -A1 " + i + " " + sys.argv[1] + " | grep -v \"^>\" "), shell=True)).rstrip("\n")).upper()
	#print(i + "\t" + str(seq.count("G") + seq.count("C")) + "\t" + str(len(seq)))
	gc_val = float((seq.count("G") + seq.count("C"))/len(seq))
	chromo_seq += [[i,len(seq),seq,gc_val]]
	

alve_locs = open(sys.argv[2]).read().rstrip("\n").split("\n")

out = open("alve_gc_20bp_window.txt", "w")

for alve in alve_locs:
	ch = alve.split("\t")[0]
	pos = int(alve.split("\t")[1])
	
	ch_slice = 0
	i = 0
	for c in chromo_seq:
		if c[0] == ch:
			ch_slice = i
			break
		i += 1
	
	pre_pos = 10
	if (pos-10) < 0:
		pre_pos = pos
	
	post_pos = 10
	if ((pos+10) > chromo_seq[ch_slice][1]):
		post_pos = chromo_seq[ch_slice][1] - pos
	
	print(chromo_seq[ch_slice][1], pos, pre_pos, post_pos)
	
	window_size = pre_pos + post_pos
	window_st = pos - pre_pos
	window_end = pos + post_pos
	
	window_gc = ((chromo_seq[ch_slice][2][window_st:window_end]).count("G") + (chromo_seq[ch_slice][2][window_st:window_end]).count("C")) / window_size
	#print(str(window_gc) + "\t" + str(chromo_seq[ch_slice][3]))
	out.write(str(window_gc) + "\t" + str(chromo_seq[ch_slice][3]) + "\n")
	
out.close()

