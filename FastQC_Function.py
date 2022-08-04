# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 02:27:24 2022

@author: Mohamad Idrees
"""
import numpy as np
from itertools import chain
import matplotlib.pyplot as plot
import statistics


def Phred33toQ(Ascii):
    Ascii=ord(Ascii)
    Q=Ascii-33
    return Q


"""
Per_Base_Sequence_Quality

"""
def Per_Base_Sequence_Quality(File):
    file=open(File)
    Q_Total=[]
    for Line in file:
        Q_Line=[]
        if Line.startswith('+'):
            L=file.readline()
            L=L.rstrip()
            for base in L:
                Q=Phred33toQ(base)
                Q_Line.append(Q)
            Q_Total.append(Q_Line)
    return Q_Total,Q_Line

Q_Total=[]
Q_Line=[]
Q_Total,Q_Line=Per_Base_Sequence_Quality('data(lecture 3).fastq')
Q_Total=np.array(Q_Total)
median_Q_of_all_reads=[]
mean_Q_of_all_reads=[]
for i in range (len(Q_Total[0])):
    Q_List=list(chain.from_iterable(Q_Total[0:,i:i+1]))
    mean=statistics.mean(Q_List)
    median=statistics.median(sorted(Q_List))
    median_Q_of_all_reads.append(median)
    mean_Q_of_all_reads.append(mean)
    
    
#plot.xlabel('Base index')
#plot.ylabel('Quality Score')
#plot.title('Per Base Sequence Histogram')
#plot.bar(range(len(mean_Q_of_all_reads)),mean_Q_of_all_reads)
#plot.show()                               
 
#---------------------------------------------------------------


"""
Per_Sequence_Quality_Scores
"""             

def Per_Sequence_Quality_Scores(File):
    file=open(File)
    Q_Total=[]
    for Line in file:
        Q_Line=[]
        if Line.startswith('+'):
            L=file.readline()
            L=L.rstrip()
            for base in L:
                Q=Phred33toQ(base)
                Q_Line.append(Q)
            mean=statistics.mean(Q_Line)
            mean=round(mean)
            Q_Total.append(mean)
    return Q_Total


qualities=Per_Sequence_Quality_Scores('data(lecture 3).fastq')
#print(qualities)
def createHist(qualities):
    range_Q=[0]*42    
    for quality in qualities:
        range_Q[quality]+=1
            
    return range_Q

Result=createHist(qualities)

#plot.xlabel('Q index')
#plot.ylabel('Number of reads')
#plot.title('Per_Sequence_Quality_Scores Histogram')
#plot.bar(range(len(Result)),Result)
#plot.show() 

#------------------------------------------------
"""
    Per Base Sequence Content
    
"""
"""
def Per_Base_Sequence_Content(File):
    file=open(File)
    Seq_Total=[]
    for Line in file:
        if Line.startswith('@'):
            L=file.readline()
            L=L.rstrip()
            Seq_Total.append(L)
    return Seq_Total

all_Seq=Per_Base_Sequence_Content('data(lecture 3).fastq')

for i in range (len(all_Seq[0])):
    Char_list=list(chain.from_iterable(all_Seq[0:,i:i+1]))
    print(Char_list)
    break;
    #mean=statistics.mean(Q_List)
    #median=statistics.median(sorted(Q_List))
    #median_Q_of_all_reads.append(median)
    #mean_Q_of_all_reads.append(mean)


"""      

#-----------------------------------------------------------
"""
   Per Sequence GC Content
"""  
from Bio.SeqUtils import GC
       
def Per_Sequence_GC_Content(File):
    file=open(File)
    Seq_Total=[]
    for Line in file:
        if Line.startswith('@'):
            L=file.readline()
            L=L.rstrip()
            Seq_Total.append(L)
    return Seq_Total

all_Seq=Per_Sequence_GC_Content('data(lecture 3).fastq')

def GC_calc(all_Seq):
    GC_seq=[]
    for seq in all_Seq:
        GC_seq.append(round(GC(seq)))
    return GC_seq
GC_seq=GC_calc(all_Seq)
range_GC=[0]*len(GC_seq)


for base in GC_seq:
    range_GC[base]+=1    
        

#plot.xlabel('GC Contnt')
#plot.ylabel('Number of Seq')
#plot.title('Per_Sequence_GC_Content Histogram')
#plot.bar(GC_seq,range(len(range_GC)))
#plot.show() 

#-------------------------------------------

"""
   Per Base N Content
""" 
"""
def Per_Sequence_N_Content(File):
    file=open(File)
    Seq_Total=[]
    for Line in file:
        if Line.startswith('@'):
            L=file.readline()
            L=L.rstrip()
            Seq_Total.append(L)
    return Seq_Total

all_Seq=Per_Sequence_N_Content('data(lecture 3).fastq')

def Calc_N(all_Seq):
    N_arr=[0]*len(all_Seq[0])
    for index in range(len(all_Seq[0])):
        for seq in all_Seq:
            if seq[index] == 'N':
                N_arr[index]+=1
                
    return N_arr

N_arr=Calc_N(all_Seq)

plot.xlabel('Base index')
plot.ylabel('Number of N')
plot.title('Per_Base_N_Content Histogram')
plot.bar(range(len(N_arr),N_arr))
plot.show() 
"""

#--------------------------------------------------
"""
   Sequence Length Distribution
""" 
def Per_Sequence_Length_Content(File):
    file=open(File)
    Seq_Total=[]
    for Line in file:
        if Line.startswith('@'):
            L=file.readline()
            L=L.rstrip()
            Seq_Total.append(L)
    return Seq_Total

all_Seq=Per_Sequence_Length_Content('data(lecture 3).fastq')

def Calc_length(all_Seq):
    seq_length=[0]*len(all_Seq[0])
    for seq in all_Seq:
       L=len(seq)
       seq_length[L]+=1
                
    return seq_length

seq_length=Calc_length(all_Seq)

#----------------------------------------