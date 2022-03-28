# -*- coding: utf-8 -*-
"""
Created on Sun Mar 14 15:55:13 2021

@author: guest1
"""
import sys
import math
import random

####
def complement(s):
    basecomplement = {
         "A":"T",
          "T":"A",
          "G":"C",
          "C":"G",
          "a":"t",
          "t":"a",
          "g":"c",
          "c":"g",
          "N":"N",
          "n":"n"
          }
    letters = list(s)
    letters = [basecomplement[base] for base in letters]
    return ''.join(letters)

def reverse(seq):
    seq_r = seq[::-1]
    sq = complement(seq_r)
    return sq

####input the mitochondria genome####

seq=''
with open(sys.argv[1],'r') as ip:
    for line in ip:
        if line.startswith('>'):
            line=line.split(' ')
            name=line[0][1:]+'_reads'
        else:
            line=line.strip()
            seq += line
            
#### calculate parameters for generating reads###

### depth： depth=(reads length * reads number)/genome size
len_seq = len(seq)
depth = int(sys.argv[2])
reads_l = int(sys.argv[3])

r_num=math.ceil((depth * len_seq)/reads_l)  # the number of reads we want to generate.

sites = list(range((len(seq) - reads_l +1))) # represents all the possible sites of the beginning of reads
fr = [1,2]      # represent the direction of the reads   

####
count = 0
cnt_f =0
cnt_r =0
with open(name + '_f.fasta','w') as op1:
    with open(name + '_r.fasta','w') as op2:
        while count < r_num:
            reads_s= random.choice(sites)
            reads = seq[reads_s:(reads_s + read_l)]
            flag = random.choice(fr)
            if flag == 2:
                cnt_r +=1
                title = '>reads_' + str(cnt_r) + '\n'
                reads = reverse(reads)
                read = reads + '\n'
                op2.write(title)
                op2.write(read)
            else:
                cnt_f +=1
                title = '>reads_' + str(cnt_f) + '\n'
                read = reads + '\n'
                op1.write(title)
                op1.write(read)
            count += 1
            




