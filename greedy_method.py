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
def theta(a, b):
    if a == '-' or b == '-':   # gap or mismatch
        return -3.0
    elif a == b:                         # match
        return 1.0
    elif a != b:
        return -1.0
def make_score_matrix(seq1, seq2):
    """
    return score matrix and map(each score from which direction)
    0: diagnosis
    1: up
    2: left
    """
    seq1 = '-' + seq1
    seq2 = '-' + seq2
    score_mat = {}
    trace_mat = {}
    picked=0
    for i,p in enumerate(seq1):
        score_mat[i] = {}
        trace_mat[i] = {}
        for j,q in enumerate(seq2):
            score_mat[0][0] = 0
            if i == 0 :                    # first row, gap in seq1
                score_mat[i][j] = -j*3
                trace_mat[i][j] = 1
                continue
            if j == 0 :                    # first column, gap in seq2
                score_mat[i][j] = -i*3
                trace_mat[i][j] = 2
                continue
            ul = score_mat[i-1][j-1] + theta(p, q)     # from up-left, mark 0
            if trace_mat[i][j-1] == 1:
                l = score_mat[i][j-1] -1
            else:
                l  = score_mat[i][j-1]   + theta('-', q)   # from left, mark 1, gap in seq1
            if trace_mat[i-1][j] == 2:
                u = score_mat[i-1][j] -1
            else:
                u  = score_mat[i-1][j]   + theta(p, '-')   # from up, mark 2, gap in seq2
            picked = max([ul,l,u])
            score_mat[i][j] = picked
            trace_mat[i][j] = [ul, l, u].index(picked)   # record which direction
    score_align=picked
    return score_mat, trace_mat, score_align

def traceback(seq1, seq2, trace_mat):
    '''
    find one optimal traceback path from trace matrix, return path code
    -!- CAUTIOUS: if multiple equally possible path exits, only return one of them -!-
    '''
    seq1, seq2 = '-' + seq1, '-' + seq2
    i, j = len(seq1) - 1, len(seq2) - 1
    path_code = ''
    while i > 0 or j > 0:
        direction = trace_mat[i][j]
        if direction == 0:                    # from up-left direction
            i = i-1
            j = j-1
            path_code = '0' + path_code
        elif direction == 1:                  # from left
            j = j-1
            path_code = '1' + path_code
        elif direction == 2:                  # from up
            i = i-1
            path_code = '2' + path_code
    return path_code

def print_m(seq1, seq2, m):
    """print score matrix or trace matrix"""
    seq1 = '-' + seq1; seq2 = '-' + seq2
    print()
    print(' '.join(['%3s' % i for i in ' '+seq2]))
    for i, p in enumerate(seq1):
        line = [p] + [m[i][j] for j in range(len(seq2))]
        print(' '.join(['%3s' % i for i in line]))
    print()
    return

def pretty_print_align(seq1, seq2, path_code):
    align1 = ''
    middle = ''
    align2 = ''
    for p in path_code:
        if p == '0':
            align1 = align1 + seq1[0]
            align2 = align2 + seq2[0]
            if seq1[0] == seq2[0]:
                middle = middle + '|'
            else:
                middle = middle + ' '
            seq1 = seq1[1:]
            seq2 = seq2[1:]
        elif p == '1':
            align1 = align1 + '-'
            align2 = align2 + seq2[0]
            middle = middle + ' '
            seq2 = seq2[1:]
        elif p == '2':
            align1 = align1 + seq1[0]
            align2 = align2 + '-'
            middle = middle + ' '
            seq1 = seq1[1:]

 #   print('Alignment:\n\n   ' + align1 + '\n   ' + middle + '\n   ' + align2 + '\n')
    lenth=len(align2)
    return lenth,align1,align2

def NW_match(seq1,seq2):
    score_mat, trace_mat, score_align= make_score_matrix(seq1, seq2)
    #print_m(seq1, seq2, score_mat)
    #print_m(seq1, seq2, trace_mat)

    path_code = traceback(seq1, seq2, trace_mat)
    result=pretty_print_align(seq1, seq2, path_code)
    lenth = result[0]
    align1 =result[1]
    align2 = result[2]
    #print('   '+path_code)
    return score_align,align1,align2

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
            
score = (-2)*(len(seq_list[0]))
lth= len(seq_list)
while len(seq_list) >1:
    t1 = datetime.datetime.now().microsecond
    t3 = time.mktime(datetime.datetime.now().timetuple())
    for i in range(1,len(seq_list)):
        results=NW_match(seq_list[0],seq_list[i])
       #score=results=[0]
       #align1=results[1]
       #align2=results[2]
        if score <= results[0]:
            score = results[0] 
            t_align = results[1]
            m_align = results[2]
            loc = i
    combine_list=[]
    for j in range(0,len(t_align)):
        if t_align[j] == '-':
            combine_list.append(m_align[j])
        elif m_align[j] == '-':
            combine_list.append(t_align[j])
        else:
            combine_list.append(t_align[j])
    cb_seq = ''.join(combine_list)
    seq_list[0] = cb_seq
    seq_list.pop(loc)
    score = (-2)*(len(seq_list[0]))
    print('done:'+ str( lth - len(seq_list)) + '/' + str(lth),cb_seq)
    t2 = datetime.datetime.now().microsecond
    t4 = time.mktime(datetime.datetime.now().timetuple())
    strTime = 'funtion time use:%dms' % ((t4 - t3) * 1000 + (t2 - t1) / 1000) 
    print(strTime)
with open(sys.argv[3],'w') as op:
    print(seq_list,cb_seq)
    op.write(cb_seq + '\n')
    










