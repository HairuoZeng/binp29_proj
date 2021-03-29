# -*- coding: utf-8 -*-
"""
Created on Sun Mar 14 16:03:03 2021

@author: guest1
"""
import sys
import datetime
import time


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

##################################
def get_weight(s1,s2):
    match = 1
    l = min(len(s1),len(s2))
    i=1
    while i < l:
        if s2[:i] == s1[-i:]:
            match = i
        else:
            match = i-1
            break
        i += 1
    
    return match

def print_result(s1,s2):   
    weight = get_weight(s1,s2)
    s = s1 + s2[weight:]
    #print(s)
    return s
#### import reads from selected file:
seq_list = []
with open(sys.argv[1],'r') as ip1:
    for line in ip1:
        line =line.strip()
        if line.startswith('>'):
            continue
        else:
            seq_list.append(line)
with open(sys.argv[2],'r') as ip2:
    for line in ip2:
        line =line.strip()
        if line.startswith('>'):
            continue
        else:
            line=reverse(line)
            seq_list.append(line)
            
seq_list = list(set(seq_list))
           
score = 0
kmer_list=[]
for i in seq_list:
    for n in range(31,len(i)):
        kmer= i[(n-31):n]
        kmer_list.append(kmer)
lth= len(kmer_list)
while len(kmer_list) > 1:
    t1 = datetime.datetime.now().microsecond
    t3 = time.mktime(datetime.datetime.now().timetuple())
    for i in range(1,len(kmer_list)):
        s1 = kmer_list[0]
        s2 = kmer_list[i]
        results1 = get_weight(s1,s2)
        results2 = get_weight(s2,s1)
        results = max([results1,results2])
        if score < results:
            score = results 
            t_align = kmer_list[0]
            m_align = kmer_list[i]
            if results == results1:
                cb_seq=print_result(t_align,m_align)
            else:
                cb_seq=print_result(m_align,t_align)
            loc = i
            print(loc)
    print(score)
    if score == 0:
        print( 'unmaped sequences exist')
        break
    kmer_list[0] = cb_seq
    kmer_list.pop(loc)
    score = 0
    print('done:'+ str( lth - len(kmer_list)) + '/' + str(lth))
    t2 = datetime.datetime.now().microsecond
    t4 = time.mktime(datetime.datetime.now().timetuple())
    strTime = 'funtion time use:%dms' % ((t4 - t3) * 1000 + (t2 - t1) / 1000) 
    print(strTime)
with open('results.fasta','w') as op:
    print(seq_list[1:],cb_seq)
    op.write(cb_seq + '\n')
    










