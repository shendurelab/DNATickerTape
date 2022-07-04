# Update in 2022.07.04: This script was used only in the preprint, and not in the revised and peer-reviewed manuscript, 
# because we substituted previous bulk lineage tracing experiment with single-cell RNA-seq assay.




# 2021.09.15 Python script to generate bigram matrix, unigram order, and bigram order for 5xTAPE-1 experiments
# Start with 16 PCR replicates of paired-end sequencing files (post UMI-amplicon collapsing)


# importing tools
import gzip
import subprocess
import os,sys,csv,re
import itertools
from collections import Counter
from collections import OrderedDict
from operator import itemgetter
from datetime import datetime
from optparse import OptionParser,OptionGroup
import numpy as np

# Additional functions
def reverse_complement(seq):
    return complement(seq[::-1])

def complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'} 
    bases = list(seq) 
    bases = [complement[base] for base in bases] 
    return ''.join(bases)




# Setting input file names
read_file = 'iPE2H3-SRv1x5-D32-rep'

        
# Pattern match for TargetBC-5xTAPE-1 GCACG(.*?)TGATGG.*?AGCACG(.*?)TGATGG
regex = r'.*?TCCATGGTCAAT(.*?)CAAGGAAAGGAA.*?'
pattern_LT = re.compile(regex)

read_table = {}
Total_reads = 0
Clean_reads = 0


# Going through each fastq.gz files == REP1
for ii in range(1,17):
    read_file_rep = read_file + str(ii) + '.assembled.fastq.gz_R1.fastq.gz'

    print ("===== Reading Forward reads for REP %i ===== \n" % ii)
    cmd = "zcat "+read_file_rep
    with os.popen(cmd) as handle: 
        handle = iter(handle)
        for line in handle:
            if line[0] == '@':
                read = next(handle).rstrip('\n')
                Total_reads +=1
                inserts = pattern_LT.match(read)
                if inserts:
                    Target_TAPE = inserts.groups()[0][0:8] + ',' + inserts.groups()[0][8:]
                    Clean_reads +=1
                    try:
                        read_table[Target_TAPE] += 1
                    except KeyError:
                        read_table[Target_TAPE] = 1
    handle.close()


read_table_ordered = OrderedDict(sorted(read_table.items(), key=lambda t: t[0], reverse=False))


# Writing to a csv file
out_file = 'iPE2H3-SRv1x5-D32-repAll_bcTAPE_all.csv'

with open(out_file, 'w', newline='') as f0:
    f0.write(read_file+','+str(Total_reads)+','+str(Clean_reads)+'\n')
    f0.write('TargetBC,TAPEsequence,rep1,rep2,rep3\n')        
    for key,value in read_table_ordered.items():
        if value > 0:
            f0.write(key+','+str(value)+'\n')
f0.close()

