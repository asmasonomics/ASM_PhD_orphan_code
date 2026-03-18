### python cluster_counter.py [-h] [-i] [-r ref_genome] pos_file 

from __future__ import division
import argparse
import sys
import subprocess
from collections import OrderedDict
import itertools

parser = argparse.ArgumentParser(
	description="cluster_counter takes a file of positions and calculates whether elements are found within clusters.",
	epilog="Author: Andrew Mason; Release: 29/04/15; Contact: andrew.mason@roslin.ed.ac.uk"
	)
parser.add_argument("pos_file", help="Specify the positions file for analysis")
parser.add_argument("-i", "--identity", action="store_true", help="select to report the identity - always in the final column of the pos file")
parser.add_argument("-r", "--ref_genome", help="reference genome for sequence extraction, only needed if the -i flag is used")
parser.add_argument("--element_num", type=int, default=5, help="(int) select min number of elements needed for a cluster (default = 5)")
parser.add_argument("--region_size", type=int, default=1000000, help="(int) select max range for a cluster (default = 1000000, 1Mbp)")
usr_args = parser.parse_args()

sys.path.append("/nfs_netapp/v1amaso4/lib/python2.7/mason-code")
import static_functions
if usr_args.identity:
    ref_dict = static_functions.seq_dict_creator(usr_args.ref_genome)

total_hits = (subprocess.check_output("wc -l " + usr_args.pos_file + " | awk \'{print $1}\'", shell=True)).rstrip("\n")
elements_by_chro = static_functions.list_formatter((open(usr_args.pos_file).read()), 0)
elements = list(elements_by_chro)
i=0
for chro in elements_by_chro:
    # default cluster is defined as 5 elements within a 1 million bp window, therefore if a contig has fewer than 5 hits it can be discarded
    if (len(chro) < int(usr_args.element_num)):
        elements.pop(i)
    else:
        i+=1

refreshed_elements = list(elements)
k=0
for chro in elements:
    elements_to_remove = []
    i=0
    for element in chro:
        # remove individual positions more than 1Mbp away from nearest other position
        if i==0:
            if ((int(element[2]) + int(usr_args.region_size)) < int(chro[i+1][1])):
                elements_to_remove += [i]
            i+=1
        elif i==(len(chro)-1):
            if ((int(element[1]) - int(usr_args.region_size)) > int(chro[i-1][2])):
                elements_to_remove += [i]
            i+=1
        else:
            if not (((int(element[2]) + int(usr_args.region_size)) > int(chro[i+1][1])) or ((int(element[1]) - int(usr_args.region_size)) < int(chro[i-1][2]))):
                elements_to_remove += [i]
            i+=1

    remove = elements_to_remove[::-1]
    for num in remove:
        chro.pop(num)

    # again, remove conitgs with fewer than 5 possible cluster positions
    if (len(chro) < int(usr_args.element_num)):
        refreshed_elements.pop(k)
    else:
        k+=1

## calculate clustered elements
clustered_elements = 0
cluster_number = 0
for chro in refreshed_elements:
    i=0
    chro_clust = []
    for element in chro:
        start = int(element[1]) - int(usr_args.region_size)
        if start < 0:
            start=0
        end = int(element[2]) + int(usr_args.region_size)

        nearby_elements = []
        finish = len(chro) - 1
        k=0
        while k < finish:
            if k < i:
                if (start < int(chro[k][2])):
                    nearby_elements += [[k, "B",(int(element[1]) - int(chro[k][2]))]]
                k+=1
            elif k == i:
                k=i+1
            elif k > i:
                if (end > int(chro[k][1])):
                    nearby_elements += [[k, "A",(int(chro[k][1]) - int(element[2]))]]
                else:
                    break
                k+=1
        i+=1

        if len(nearby_elements) >= (int(usr_args.element_num) - 1):
            
            # all nearby elements are in the 1 million bases after the position
            if nearby_elements[0][1] == "A":
                for e in nearby_elements:
                    chro_clust += [chro[(e[0])]]
                    chro_clust += [element]
            # all nearby elements are in the 1 million bases before the position
            elif nearby_elements[-1][1] == "B":
                for e in nearby_elements:
                    chro_clust += [chro[(e[0])]]
                    chro_clust += [element]
            else:
                # the range between first and last elements in the list is less than or equal to 1 million
                if ((nearby_elements[0][2] + nearby_elements[-1][2]) <= int(usr_args.region_size)):
                    for e in nearby_elements:
                        chro_clust += [chro[(e[0])]]
                        chro_clust += [element]
                # now consider the list in groups of 5 to see whether or not there is 1 million base pairs covering the group
                else:
                    m=0
                    while (m+(int(usr_args.element_num)-1) < len(nearby_elements)):
                        # the sub group of 5 are all either before or after the element
                        if (nearby_elements[m][1] == nearby_elements[m+(int(usr_args.element_num)-1)][1]):
                            for e in nearby_elements[m:m+(int(usr_args.element_num)-1)]:
                                chro_clust += [chro[(e[0])]]
                                chro_clust += [element]
                        # remaining subgroups are combinations of B and A, therefore a simple range can be used
                        elif ((nearby_elements[m][2] + nearby_elements[m+(int(usr_args.element_num)-1)][2]) <= int(usr_args.region_size)):
                            for e in nearby_elements[m:m+(int(usr_args.element_num)-1)]:
                                chro_clust += [chro[(e[0])]]
                                chro_clust += [element]
                        m+=1

    print("\nChromosome " + str(chro[0][0]))
    inserts = (subprocess.check_output("awk \'$1==\"" + str(chro[0][0]) + "\"{print $0}\' " + usr_args.pos_file + " | wc -l", shell=True)).rstrip("\n")
    print(str(inserts) + " positions analysed")

    if len(chro_clust) > 0:
        # sort the cluster matches and remove duplicates
        chro_clust.sort(lambda l1, l2: cmp(l1[1], l2[1]) or (l1[1] == l2[1] and (cmp(l1[2], l2[2]))))
        sorted_chro_clust = list(chro_clust for chro_clust,_ in itertools.groupby(chro_clust))

        clusters=[]
        q=0
        k=1
        for match in sorted_chro_clust:
            if (q==(len(sorted_chro_clust)-1)):
                coord = "(" + str(sorted_chro_clust[(q+1)-k][1]) + "-" + str(str(match[2])) + " : " + str(format(((match[2] + 1) - (sorted_chro_clust[(q+1)-k][1])), ",d")) + ")" 
                clusters += [k,coord]
                k=0
            elif (((sorted_chro_clust[q+1][1])-int(usr_args.region_size)) > match[2]):
                coord = "(" + str(sorted_chro_clust[(q+1)-k][1]) + "-" + str(str(match[2])) + " : " + str(format(((match[2] + 1) - (sorted_chro_clust[(q+1)-k][1])), ",d")) + ")" 
                clusters += [k,coord]
                k=0
            k+=1
            q+=1

        cluster_num = int(len(clusters) / 2)
        cluster_sizes = str(clusters).replace(" ", "").replace("[", "").replace("]", "").replace("\'", "").replace("\"", "").replace(",(", " (").replace("),", "); ").replace(":", " : ")
        
        perc_per_chro = round(((float(100/int(inserts)))*len(sorted_chro_clust)), 2)
        print(str(len(sorted_chro_clust)) + " elements within " + str(cluster_num) + " cluster(s) --> " + str(perc_per_chro) + "% ; cluster size(s): " + str(cluster_sizes))

        if usr_args.identity:
            identity_scores = []
            for hit in sorted_chro_clust:
                identity_scores += [float(hit[-1])]
            tot_identity = len(sorted_chro_clust)
            cluster_lengths = []
            for item in str(clusters).replace("[", "").replace("]", "").split(", "):
                if item.isdigit():
                    cluster_lengths += [int(item)]
                    
            cluster_identities = []
            cluster_members = []
            z=0
            tot = 0
            while z < len(cluster_lengths):
                start = tot
                tot += cluster_lengths[z]
                cluster_identities += [identity_scores[start:tot]]
                cluster_members += [sorted_chro_clust[start:tot]]
                z+=1

            z=1
            for cluster in cluster_identities:
                print("Cluster " + str(z) + ": " + str(cluster).replace("[", "").replace("]", ""))
                file_name_core = str(chro[0][0]) + "_cluster" + str(z).zfill(2)
                clust_pos = open((file_name_core + "_pos.txt"), "w")
                for hit in cluster_members[z-1]:
                    clust_pos.write(str(hit[0]) + "\t" + str(hit[1]) + "\t" + str(hit[2]) + "\t" + str(hit[3]) + "\t" + str(hit[4]) + "_" + str(hit[5]) + "\n")
                clust_pos.close()
                static_functions.seq_extract((file_name_core + "_seq.fasta"), (static_functions.list_initial_formatter(file_name_core + "_pos.txt")), ref_dict, "N")
                z+=1

        clustered_elements += len(sorted_chro_clust)
        cluster_number += cluster_num
    else:
        print("No clusters found")

# print overall output messages
print("\nAnalysed a total of " + str(total_hits) + " insertion sites")
print(str(clustered_elements) + " elements identified within " + str(cluster_number) + " cluster(s)")
perc = round((float(100 / int(total_hits))*clustered_elements), 2)
print(str(perc) + "% of elements found within clusters\n")
                




