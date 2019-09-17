#!/usr/bin/env python
# encoding: utf-8

"""
show max values in list in population (major allele)
"""

import os,sys,string,gzip,argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-i','--input',help="Input vcf file name , missing value are code as -10000",required=False)
parser.add_argument('-o','--out',help='Output file name', required=True)
parser.add_argument('-p1','--pop1',help='population1 file name', required=True)

args = parser.parse_args() 
p1=args.pop1
popname1=p1.split("/")[len(p1.split("/"))-1]



f = open(args.input, 'r')
mf = f.readlines()


Specific = open(args.out+"_"+popname1+'_all.showmax.out','w')

for line in mf:
    if line.split("\t")[0] == 'CHROM':
        header = line.replace("\n","").split("\t")
    idnames = header[0:]

if p1 != None:
    ids1 = [i.replace("\n","") for i in open(p1,"r")]  
    idx2keep1=[idnames.index(i) for i in ids1]

  
def showmax(lt):
    index1 = 0
    max = 0
    for i in range(len(lt)):
        flag = 0
        for j in range(i+1,len(lt)):
            if lt[j] == lt[i]:
                flag += 1
        if flag > max:
            max = flag
            index1 = i
    return lt[index1]

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
        ab1 = []
        for num in b1:
            if(num == "-10000"):
                continue
            else:
                ab1.append(int(num))
        bb = showmax(ab1)
        Specific.write(chr+'\t'+pos+'\t'+str(bb)+'\n')

Specific.close()

