# -*- coding: utf-8 -*-
"""
Created on Tue Mar 16 09:30:06 2021
@author: guest1
"""

from functools import reduce
import sys

def Mapping_score(s1,s2):    
    length = min(len(s1),len(s2))
    while length>0:
        if s2[:length] == s1[-length:]:
            return length
        else:
            length -= 1
    return 0

def Num_indegree(graph,v): 
    nodes_in = list(graph.values())
    nodes_list=[]
    for i in nodes_in:
        nodes_list += i
    count = nodes_list.count(v)
    return count

def Dir_graph(seq_list,length,t):    
    graph = {}
    cost_graph = {}
    for i in seq_list:
        i_node = []
        count_i_node={}
        for j in seq_list:
            if i != j:
                score = Mapping_score(i,j)
                if score > t:
                    i_node.append(j)
                    count_i_node[j] = length - score
        graph[i] = i_node
        cost_graph[i] = count_i_node
    keys_graph=list(graph.keys())
    for i in keys_graph:  
        if graph[i] == []: #outdegree='0'
            count = Num_indegree(graph,i)
            if count ==0:
                graph.pop(i)
                cost_graph.pop(i)
                print('The sequence:\n"{0}"\n can\'t align with others!'.format(i))
    return graph,cost_graph

def BFS_topo(graph):  
    topo=[]
    while len(graph) > 0:
        All_nodes = graph.keys()
        for i in All_nodes:
            flag = True
            in_num = Num_indegree(graph,i)
            if in_num ==0:
                topo.append(i)
                graph.pop(i)
                flag = False
                break
        if flag == True:                       
            print('The t score is too small!')
            return None
    return topo


def Reads_alignment(s1,s2): 
    weight = Mapping_score(s1,s2)
    s = s1 + s2[weight:]
    return s

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
def Find_startpoint(graph):
    nodes = graph.keys()
    for i in nodes:
        if Num_indegree(graph,i) == 0 and graph[i]:
            start_point=i
    return start_point

def Dijkstra_algorithm(graph):
    start_point=Find_startpoint(graph)
    distance = graph[start_point]
    trace=[start_point]
    min_node = ''
    nodes= list(graph.keys())
    for i in nodes:
        if i in distance or i == start_point:
            continue
        else:
            distance[i]= 9999999
    count = len(cost_graph)
    while count > 1:
        sort_dis = sorted(distance.items(), key=lambda item: item[1])
        for p in sort_dis:
            if p[0] not in trace:
                min_node=p[0]
                min_distance=distance[min_node]
                trace.append(min_node)
        for j in cost_graph[min_node]:
            new_path = min_distance+ cost_graph[min_node][j]
            if distance[j] > new_path:
                distance[j] = new_path
    count -=1
    
    return trace
    
#%%
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
            
threshold=int(sys.argv[3])


length = len(seq_list[0])

graph,cost_graph= Dir_graph(seq_list,length,threshold)

checkpoint=str(input("Directed Graph Assemble methods: Dijkstra or BFS?"))

if checkpoint == 'Dijkstra':
    topo_list= Dijkstra_algorithm(cost_graph)
    method = 'Dijkstra_algorithm'
elif checkpoint == 'BFS':
    topo_list= BFS_topo(graph)
    method = 'BFS'
else:
    print('No matched algorithm. Using BFS instead')
    topo_list= BFS_topo(graph)
    method = 'BFS'
    
result = reduce(Reads_alignment,topo_list)

name= method+'_'+sys.argv[4]+'.txt'
with open(name,'w') as op:
    op.write(str(result))




