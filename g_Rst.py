#!/usr/bin/env python
# encoding: utf-8

"""

Rst for STR (M. Slatkin 1995, A Measure of Population Subdivision Based on Microsatellite Allele Frequencies)
< Rst = (St - Sw) / St >
Sw = average sum of squares of the differences in allele size within each population
     is twice the average of the estimated variances of allele size within each population
St = is twice the estimated variance in allele size in the collection of populations together
Sb = average squared difference between all pairs of copies we define the between-population 

"""

import os,sys,string,gzip,argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-i','--input',help="Input vcf  file name ,copy number genotypes, missing value are code as -10000",required=False)
parser.add_argument('-o','--out',help='Output file name', required=True)
parser.add_argument('-p1','--pop1',help='population1 file name', required=True)
parser.add_argument('-p2','--pop2',help='population2 file name', required=True)

args = parser.parse_args() 
p1=args.pop1
p2=args.pop2
popname1=p1.split("/")[len(p1.split("/"))-1]
popname2=p2.split("/")[len(p2.split("/"))-1]


f = open(args.input, 'r')
mf = f.readlines()


allRstOut = open(args.out+"_"+popname1+"_"+popname2+'_all.Rst.out','w')

for line in mf:
    if line.split("\t")[0] == 'CHROM':
        header = line.replace("\n","").split("\t")
    idnames = header[0:]

if p1 != None:
    ids1 = [i.replace("\n","") for i in open(p1,"r")]  
    idx2keep1=[idnames.index(i) for i in ids1]

if p2 != None:
    ids2 = [i.replace("\n","") for i in open(p2,"r")]   
    idx2keep2=[idnames.index(i) for i in ids2]

for line in mf:
    if line.split("\t")[0] != 'CHROM':
        a=line.replace("\n","").replace("/","\t").split("\t")
        chr=a[0].replace("chr","")
        pos=a[1]
        a1=[]
        for i in idx2keep1:
            a1.append(a[2*i-2])
            a1.append(a[2*i-1])
        b1=a1
        a2=[]
        for i in idx2keep2:
            a2.append(a[2*i-2])
            a2.append(a[2*i-1])
        b2=a2
        ab1 = []
        for num in b1:
            if(num == "-10000"):
                continue
            else:
                ab1.append(round(float(num),3))
                
        ab2 = []
        for num in b2:
            if(num == "-10000"):
                continue
            else:
                ab2.append(round(float(num),3))
        ab3=ab1+ab2
        d1 = 2*round(np.var(ab1),6)
        d2 = 2*round(np.var(ab2),6)
        d3 = 2*round(np.var(ab3),6)
        if d3 != 0 and len(ab1) >= 1.2*len(ids1) and len(ab2) >= 1.2*len(ids2):
            dd = (d1*len(ab1)+d2*len(ab2))/(len(ab1)+len(ab2))
            Rst = (d3-dd)/d3
            Rst = round(Rst,6)
        else :
            Rst = "inf"
        allRstOut.write(chr+'\t'+pos+'\t'+str(Rst)+'\n') 

allRstOut.close()

