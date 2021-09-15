# 2021.09.15 Python script to generate bigram matrix, unigram order, and bigram order for 5xTAPE-1 experiments



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
from scipy.cluster.vq import vq, kmeans, whiten
from scipy.cluster.vq import kmeans2


# Additional functions
def reverse_complement(seq):
    return complement(seq[::-1])

def complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'} 
    bases = list(seq) 
    bases = [complement[base] for base in bases] 
    return ''.join(bases)


# Processes 3 transfection replicates to make Unigram and Bigram dictionaries
def TAPEx5_unigram_bigram_reps(read_file, insert_size):
    
    # Initialize
    Site1_insert_table = {}
    Site2_insert_table = {}
    Site3_insert_table = {}
    Site4_insert_table = {}
    Site5_insert_table = {}
    Site12_insert_table = {}
    Total_reads = 0

    # Expected
    Site_XO = 0
    Site_XXO = 0
    Site_XXXO = 0
    Site_XXXXO = 0
    Site_XXXXX = 0
    
    # Pattern match for TAPE-1
    regex = r'.*?GCACG(.*?)TGATGG.*?AGCACG(.*?)TGATGG.*?AGCACG(.*?)TGATGG.*?AGCACG(.*?)TGATGG.*?AGCACG(.*?)TGATGG.*?'
    pattern_TAPEv1 = re.compile(regex)

    # Reading in the fastq.gz file for rep1
    print ("===== Reading Forward reads for REP1 ===== \n")
    cmd = "zcat "+read_file
    with os.popen(cmd) as handle: 
        handle = iter(handle)
        for line in handle:
            if line[0] == '@':
                read = next(handle).rstrip('\n')
                inserts = pattern_TAPEv1.match(read)
                if inserts:
                    Total_reads +=1
                    # First is filled
                    if len(inserts.groups()[0]) == insert_size:
                        try:
                            Site1_insert_table[inserts.groups()[0]] += 1
                        except KeyError:
                            Site1_insert_table[inserts.groups()[0]] = 1
                        
                        # First-Second are filled
                        if len(inserts.groups()[1]) == insert_size:
                            
                            insert_order12 = inserts.groups()[0] + '-' + inserts.groups()[1]
                            try:
                                Site12_insert_table[insert_order12] += 1
                            except KeyError:
                                Site12_insert_table[insert_order12] = 1
                            
                            try:
                                Site2_insert_table[inserts.groups()[1]] += 1
                            except KeyError:
                                Site2_insert_table[inserts.groups()[1]] = 1
                            
                            # First-Second-Third are filled
                            if len(inserts.groups()[2]) == insert_size:
                                
                                insert_order23 = inserts.groups()[1]+ '-' + inserts.groups()[2]
                                try:
                                    Site12_insert_table[insert_order23] += 1
                                except KeyError:
                                    Site12_insert_table[insert_order23] = 1

                                try:
                                    Site3_insert_table[inserts.groups()[2]] += 1
                                except KeyError:
                                    Site3_insert_table[inserts.groups()[2]] = 1
                                
                                # First-Second-Third-Fourth
                                if len(inserts.groups()[3]) == insert_size:
                                    
                                    insert_order34 = inserts.groups()[2]+'-'+inserts.groups()[3]
                                    try:
                                        Site12_insert_table[insert_order34] += 1
                                    except KeyError:
                                        Site12_insert_table[insert_order34] = 1
                                    
                                    try:
                                        Site4_insert_table[inserts.groups()[3]] += 1
                                    except KeyError:
                                        Site4_insert_table[inserts.groups()[3]] = 1
                                        
                                    # First-Second-Third-Fourth-Fifth
                                    if len(inserts.groups()[4]) == insert_size:
                                        insert_order45 = inserts.groups()[3]+'-'+inserts.groups()[4]
                                        try:
                                            Site12_insert_table[insert_order45] += 1
                                        except KeyError:
                                            Site12_insert_table[insert_order45] = 1
                                    
                                        try:
                                            Site5_insert_table[inserts.groups()[4]] += 1
                                        except KeyError:
                                            Site5_insert_table[inserts.groups()[4]] = 1
                                        Site_XXXXX +=1
                                    elif len(inserts.groups()[4]) == 0:
                                        Site_XXXXO +=1
                                elif len(inserts.groups()[3]) == 0:
                                    Site_XXXO +=1
                            elif len(inserts.groups()[2]) == 0:
                                Site_XXO +=1
                        elif len(inserts.groups()[1]) == 0:
                            Site_XO +=1                            
    handle.close()

    # Reading in the fastq.gz file for rep2
    print ("===== Reading Forward reads for REP2 ===== \n")
    cmd = "zcat "+read_file.replace('rep1','rep2')
    with os.popen(cmd) as handle: 
        handle = iter(handle)
        for line in handle:
            if line[0] == '@':
                read = next(handle).rstrip('\n')
                inserts = pattern_TAPEv1.match(read)
                if inserts:
                    Total_reads +=1
                    # First is filled
                    if len(inserts.groups()[0]) == insert_size:
                        try:
                            Site1_insert_table[inserts.groups()[0]] += 1
                        except KeyError:
                            Site1_insert_table[inserts.groups()[0]] = 1
                        
                        # First-Second are filled
                        if len(inserts.groups()[1]) == insert_size:
                            
                            insert_order12 = inserts.groups()[0] + '-' + inserts.groups()[1]
                            try:
                                Site12_insert_table[insert_order12] += 1
                            except KeyError:
                                Site12_insert_table[insert_order12] = 1
                            
                            try:
                                Site2_insert_table[inserts.groups()[1]] += 1
                            except KeyError:
                                Site2_insert_table[inserts.groups()[1]] = 1
                            
                            # First-Second-Third are filled
                            if len(inserts.groups()[2]) == insert_size:
                                
                                insert_order23 = inserts.groups()[1]+ '-' + inserts.groups()[2]
                                try:
                                    Site12_insert_table[insert_order23] += 1
                                except KeyError:
                                    Site12_insert_table[insert_order23] = 1

                                try:
                                    Site3_insert_table[inserts.groups()[2]] += 1
                                except KeyError:
                                    Site3_insert_table[inserts.groups()[2]] = 1
                                
                                # First-Second-Third-Fourth
                                if len(inserts.groups()[3]) == insert_size:
                                    
                                    insert_order34 = inserts.groups()[2]+'-'+inserts.groups()[3]
                                    try:
                                        Site12_insert_table[insert_order34] += 1
                                    except KeyError:
                                        Site12_insert_table[insert_order34] = 1
                                    
                                    try:
                                        Site4_insert_table[inserts.groups()[3]] += 1
                                    except KeyError:
                                        Site4_insert_table[inserts.groups()[3]] = 1
                                        
                                    # First-Second-Third-Fourth-Fifth
                                    if len(inserts.groups()[4]) == insert_size:
                                        insert_order45 = inserts.groups()[3]+'-'+inserts.groups()[4]
                                        try:
                                            Site12_insert_table[insert_order45] += 1
                                        except KeyError:
                                            Site12_insert_table[insert_order45] = 1
                                    
                                        try:
                                            Site5_insert_table[inserts.groups()[4]] += 1
                                        except KeyError:
                                            Site5_insert_table[inserts.groups()[4]] = 1
                                        Site_XXXXX +=1
                                    elif len(inserts.groups()[4]) == 0:
                                        Site_XXXXO +=1
                                elif len(inserts.groups()[3]) == 0:
                                    Site_XXXO +=1
                            elif len(inserts.groups()[2]) == 0:
                                Site_XXO +=1
                        elif len(inserts.groups()[1]) == 0:
                            Site_XO +=1                            
    handle.close()
    
    # Reading in the fastq.gz file for rep3
    print ("===== Reading Forward reads for REP3 ===== \n")
    cmd = "zcat "+read_file_x5.replace('rep1','rep3')
    with os.popen(cmd) as handle: 
        handle = iter(handle)
        for line in handle:
            if line[0] == '@':
                read = next(handle).rstrip('\n')
                inserts = pattern_TAPEv1.match(read)
                if inserts:
                    Total_reads +=1
                    # First is filled
                    if len(inserts.groups()[0]) == insert_size:
                        try:
                            Site1_insert_table[inserts.groups()[0]] += 1
                        except KeyError:
                            Site1_insert_table[inserts.groups()[0]] = 1
                        
                        # First-Second are filled
                        if len(inserts.groups()[1]) == insert_size:
                            
                            insert_order12 = inserts.groups()[0] + '-' + inserts.groups()[1]
                            try:
                                Site12_insert_table[insert_order12] += 1
                            except KeyError:
                                Site12_insert_table[insert_order12] = 1
                            
                            try:
                                Site2_insert_table[inserts.groups()[1]] += 1
                            except KeyError:
                                Site2_insert_table[inserts.groups()[1]] = 1
                            
                            # First-Second-Third are filled
                            if len(inserts.groups()[2]) == insert_size:
                                
                                insert_order23 = inserts.groups()[1]+ '-' + inserts.groups()[2]
                                try:
                                    Site12_insert_table[insert_order23] += 1
                                except KeyError:
                                    Site12_insert_table[insert_order23] = 1

                                try:
                                    Site3_insert_table[inserts.groups()[2]] += 1
                                except KeyError:
                                    Site3_insert_table[inserts.groups()[2]] = 1
                                
                                # First-Second-Third-Fourth
                                if len(inserts.groups()[3]) == insert_size:
                                    
                                    insert_order34 = inserts.groups()[2]+'-'+inserts.groups()[3]
                                    try:
                                        Site12_insert_table[insert_order34] += 1
                                    except KeyError:
                                        Site12_insert_table[insert_order34] = 1
                                    
                                    try:
                                        Site4_insert_table[inserts.groups()[3]] += 1
                                    except KeyError:
                                        Site4_insert_table[inserts.groups()[3]] = 1
                                        
                                    # First-Second-Third-Fourth-Fifth
                                    if len(inserts.groups()[4]) == insert_size:
                                        insert_order45 = inserts.groups()[3]+'-'+inserts.groups()[4]
                                        try:
                                            Site12_insert_table[insert_order45] += 1
                                        except KeyError:
                                            Site12_insert_table[insert_order45] = 1
                                    
                                        try:
                                            Site5_insert_table[inserts.groups()[4]] += 1
                                        except KeyError:
                                            Site5_insert_table[inserts.groups()[4]] = 1
                                        Site_XXXXX +=1
                                    elif len(inserts.groups()[4]) == 0:
                                        Site_XXXXO +=1
                                elif len(inserts.groups()[3]) == 0:
                                    Site_XXXO +=1
                            elif len(inserts.groups()[2]) == 0:
                                Site_XXO +=1
                        elif len(inserts.groups()[1]) == 0:
                            Site_XO +=1                            
    handle.close()

    print("===== Forward reads have been read ===== \n")



    unigram1_dict = OrderedDict(sorted(Site1_insert_table.items(), key=lambda t: t[1], reverse=True)[:120])
    unigram2_dict = OrderedDict(sorted(Site2_insert_table.items(), key=lambda t: t[1], reverse=True)[:120])
    unigram3_dict = OrderedDict(sorted(Site3_insert_table.items(), key=lambda t: t[1], reverse=True)[:120])
    unigram4_dict = OrderedDict(sorted(Site4_insert_table.items(), key=lambda t: t[1], reverse=True)[:120])
    unigram5_dict = OrderedDict(sorted(Site5_insert_table.items(), key=lambda t: t[1], reverse=True)[:120])
    
    bigram12_dict = OrderedDict(sorted(Site12_insert_table.items(), key=lambda t: t[1], reverse=True))


    
    return unigram1_dict,unigram2_dict,unigram3_dict,unigram4_dict,unigram5_dict,bigram12_dict



# Outputs unigram-orders
def TAPEx5_positional_unigram(unigram1_dict,unigram2_dict,unigram3_dict,unigram4_dict,unigram5_dict):
    insert_BC_list =[]
    bases = ['A','C','G','T']
    for i1 in bases:
        for i2 in bases:
            for i3 in bases:
                insert_BC_list.append(i1+i2+i3)
            
    #np.random.seed(1)


    # Correct unigram with editing score
    TAPE_3N3 = {}
    with open('TAPE-3N3-editScore.csv', newline='') as csvfile:
        TAPE_3N3_rows = csv.reader(csvfile, delimiter=',')
        for row in TAPE_3N3_rows:
            TAPE_3N3[row[1]] = row[2]

    unigram1_dict_corr = {}
    for key,value in unigram1_dict.items():
        ES = 2**float(TAPE_3N3[key[0:3]])
        unigram1_dict_corr[key] = value / ES


    # Using k-means to recover introduced barcodes
    unigram_counts = np.zeros((64,1))

    unigram_key = np.array(list(unigram1_dict_corr.keys()))
    unigram_value = np.array(list(unigram1_dict_corr.values()))





    # Convert read-counts to log2
    unigram_log = np.log2(unigram_value)
    uni_centroid, uni_label = kmeans2(unigram_log, 2, minit='points')
    kmean_indexes = np.nonzero(uni_label == uni_centroid.argmax())[0]

    for ii in kmean_indexes:
        if unigram_key[ii][3:] == 'GGA':
            first = unigram_key[ii][0:3]
            try:
                unigram_counts[insert_BC_list.index(first)] += unigram_value[ii]
            except ValueError:
                print('Unknown unigram barcode is: ' + first)
            
    unigram_set = np.nonzero(unigram_counts)[0]
    
    
    
    unigram_array = np.zeros((64,5))
    array_count = 0
    for key in insert_BC_list:
        current_key = key + 'GGA'
        try:
            current_U1 = unigram1_dict[current_key]
        except KeyError:
            current_U1 = '0'
        try:
            current_U2 = unigram2_dict[current_key]
        except KeyError:
            current_U2 = '0'    
        try:
            current_U3 = unigram3_dict[current_key]
        except KeyError:
            current_U3 = '0'
        try:
            current_U4 = unigram4_dict[current_key]
        except KeyError:
            current_U4 = '0'
        try:
            current_U5 = unigram5_dict[current_key]
        except KeyError:
            current_U5 = '0'
        current_U = np.array((current_U1,current_U2,current_U3,current_U4,current_U5))
        unigram_array[array_count,:] = current_U
        array_count += 1
    
    unigram_array = unigram_array[unigram_set]
    unigram_freq = unigram_array.astype('int')[:,0]/np.sum(unigram_array, axis=1)
    unigram_order = unigram_set[np.argsort(unigram_freq)[::-1][:len(unigram_set)]]
    
    return unigram_order


# Outputs the order of barcodes
def TAPEx5_bubble_sort(unigram_order, bigram_dict):
    insert_BC_list =[]
    bases = ['A','C','G','T']
    for i1 in bases:
        for i2 in bases:
            for i3 in bases:
                insert_BC_list.append(i1+i2+i3)
    
    bigram_matrix = np.zeros((64, 64))
    for key,value in bigram_dict.items():
        first =key[0:3]
        second=key[7:10]
        try:
            bigram_matrix[insert_BC_list.index(first),insert_BC_list.index(second)] += value
        except ValueError:
            print('Unknown bigram barcode is: ' + first + '-' + second)
    
    # Normalizing the matrix
    bigram_matrix_colnorm = bigram_matrix + 1
    sum_of_cols = bigram_matrix_colnorm.sum(axis=0)
    bigram_matrix_colnorm = bigram_matrix_colnorm / sum_of_cols[np.newaxis, :]
    sum_of_rows = bigram_matrix_colnorm.sum(axis=1)
    bigram_matrix_rowcolnorm = bigram_matrix_colnorm / sum_of_rows[:, np.newaxis]



    sorted_order = np.copy(unigram_order)
    print('Initial order is: \n' )
    print(sorted_order)
    
    total_swap = 0
    for jj in range(20):
        swap_count = 0
        for ii in range(unigram_order.size-1):
            if (bigram_matrix_rowcolnorm[sorted_order[ii],sorted_order[ii+1]]) / (bigram_matrix_rowcolnorm[sorted_order[ii+1],sorted_order[ii]]) < 1:
                sorted_order[[ii,ii+1]] = sorted_order[[ii+1,ii]]
                swap_count += 1
            #print(sorted_order)
    
        #print('Iteration count %d and swap count %d' % (jj,swap_count))
        total_swap += swap_count
        if swap_count == 0:
            break
        continue
    print('Final iteration count %d and swap count %d' % (jj,total_swap))
    print('Final sorted order is: \n' )
    print(sorted_order)
    return sorted_order


# Saves bigram_matrix in a long-format for plotting in R
def TAPE_bigram_long(bigram_dict,file_name):
    insert_BC_list =[]
    bases = ['A','C','G','T']
    for i1 in bases:
        for i2 in bases:
            for i3 in bases:
                insert_BC_list.append(i1+i2+i3)

    with open(file_name, 'w', newline='') as f0:
        f0.write('Position1,Position2,ReadCounts\n')
        for ii in range(0,64):
            for jj in range(0,64):
                current_barcode_pair = insert_BC_list[ii]+'GGA-'+insert_BC_list[jj]+'GGA'
                try:
                    current_RC = bigram_dict[current_barcode_pair]
                except KeyError:
                    current_RC = 0
                current_row = insert_BC_list[ii] + ',' + insert_BC_list[jj] + ',' + str(current_RC) + '\n'
                f0.write(current_row)
    f0.close()
    return()



# Command to process each data (paired-end sequencing data for 'TAPE-3NTE-M1N1s1-E4' condition is used here)
# 6-bp or NNNGGA insertion is expected
unigram1_dict,unigram2_dict,unigram3_dict,unigram4_dict,unigram5_dict,bigram12_dict = TAPEx5_unigram_bigram_reps('TAPE-3NTE-M1N1s1-E4-rep1.assembled.fastq.gz',6)

# Command to export the bigram matrix for plotting in R
TAPE_bigram_long(bigram12_dict,'TAPE-3NTE-M1N1s1-E6.csv')

# Command to get unigram and bigram orders
unigram_order_M1 = TAPEx5_positional_unigram(unigram1_dict,unigram2_dict,unigram3_dict,unigram4_dict,unigram5_dict)
bigram_order_M1 = TAPEx5_bubble_sort(unigram_order_M1, bigram12_dict)

