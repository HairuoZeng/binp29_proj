# -*- coding: utf-8 -*-
"""
Created on Tue Mar 16 09:30:06 2021

@author: guest1
"""

from functools import reduce
import sys


def get_weight(s1,s2):    
    l = min(len(s1),len(s2))
    while l>0:
        if s2[:l] == s1[-l:]:
            return l
        else:
            l-=1
    return 0

def print_result(s1,s2): 
    weight = get_weight(s1,s2)
    s = s1 + s2[weight:]
    print(s)
    return s

def dir_graph(l,t=3):    
    graph = {}
    for i in l:
        VW = []
        for j in l:
            if i!=j:
                weight = get_weight(i,j)
                if weight >= t:
                    VW.append(j)
        graph[i] = VW
    #print(graph)
    for i in graph.keys():  
        if not graph[i]:
            count = get_in_V(graph,i)
            if count ==0:
                graph.clear()
                print('The sequence:\n"{0}"\n can\'t align with others!'.format(i))
                break
    return graph

def get_in_V(graph,v): 
    count = 0
    all_in = reduce(lambda x,y:x+y,graph.values())
    for i in all_in:
        if i == v:
            count+=1
    return count

def aligner(graph,topo=[]):  
    while graph:
        V = graph.keys()
        for i in V:
            flag = 1
            in_num = get_in_V(graph,i)
            if in_num ==0:
                topo.append(i)
                graph.pop(i)
                flag = 0
                break
        if flag:                       
            print('The t score is too small!')
            return None
        else:
            aligner(graph,topo)
    return topo

####
def complement(s):
    basecomplement = {
         "A":"T",
          "T":"A",
          "G":"C",
          "C":"G",
          "a":"t",
          "t":"a",
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

#### import reads from selected file:
seq_list = []
with open('test_reads_f.fasta','r') as ip1:
    for line in ip1:
        line =line.strip()
        if line.startswith('>'):
            continue
        else:
            seq_list.append(line)
with open('test_reads_r.fasta','r') as ip2:
    for line in ip2:
        line =line.strip()
        if line.startswith('>'):
            continue
        else:
            line=reverse(line)
            seq_list.append(line)


graph = dir_graph(seq_list,t=30)
topo = aligner(graph)
print(graph)
if topo:
    result = reduce(print_result,topo)
    print(result)
else:
    result = topo
with open('lol1.txt','w') as op:
    op.write(str(result))
    
    
    