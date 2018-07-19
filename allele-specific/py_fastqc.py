#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 12 14:39:44 2018
FastQC reader in Python: Takes in fastq file and outputs graph of
                        -Quality
                        -Read Length
                        -%GC just above 50%
@author: nickgarza
"""


import matplotlib.pyplot as plt
import numpy as np

def main():
    fin = input("Please enter the fastQ file: ")
    file = open(fin,'r')
    lines = file.readlines()
    
    
    reads = []
    avgqual = []
    lengths = []
    gcArray = []
    for i in range(int(len(lines)/4)):
        read = lines[((i+1)*4 -4):(i+1)*4]
        reads.append(read)
    
    
    
    
    for i in range(len(reads)):
        values = (fastq(reads[i]))
        avgqual.append(values.getQual())
        lengths.append(values.getLength())
        gcArray.append(values.getGC())
    
   # ----------------------------------
    fig1 = plt.figure(1)
    plt.hist(avgqual)
    plt.xlabel("Average Qualities")
    plt.ylabel("Number of Reads")
    plt.title("Quality Across Sequences")
    fig1.savefig("ngarza1_quality.png")
    
    
    
    fig2 = plt.figure(2)
    plt.hist(lengths)
    plt.xlabel("Read Length")
    plt.ylabel("Number of Reads")
    plt.title("Read Lengths for FASTQ")
    fig2.savefig("ngarza1_lengths.png")
    
    
    fig3 = plt.figure(3)
    plt.hist(gcArray)
    plt.xlabel("%GC")
    plt.ylabel("Number of Reads")
    plt.title("%GC of Reads")
    fig3.savefig("ngarza1_gc.png")
    
   # --------------------------------
   
    outfile = open("ngarza1_filtered.fq",'w')
   
    
    
   
    for i in range(int(len(lines)/4)):
        values = fastq(reads[i])
        if values.getGC()>50:
            for item in reads[i]:
                print(item,file = outfile)
    
    file.close()
    
class fastq:
    
    def __init__(self, reads):
         
        self.names = reads[0] # name of line
        self.seq = reads[1]   # DNA sequence
        self.qual = reads[3].rstrip()   # Description of Quality
        
        

    def getQual(self):
        num = 0
        qual = 0
        avgList = []
        for s in self.qual:
            qual += ord(s)-33   # Gets quality from ASCII
            num+=1
            avgList.append(qual) # adds quality to list
            
        return np.mean(avgList) #returns quality averages
    
    def getGC(self):
        
        letters = self.seq
        gCount = letters.count('G') # counts the number of Gs
        cCount = letters.count('C') # counts the number of Cs
        
        return ((gCount + cCount)/len(letters)*100)
    
    
    
    def getLength(self):
        return len(self.seq)    # returns length of sequence from reads
    

       
        
        
main()
