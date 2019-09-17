#!/usr/bin/env python
# encoding: utf-8


import os,sys,string,gzip,argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-i','--inp',help="Input vcf file name,GB genotypes, missing value is code -10000",required=False)
parser.add_argument('-o','--out',help='Output file name', required=True)

args = parser.parse_args() 

f = open(args.inp, 'r')

mf = f.readlines()

allGTOut = open(args.out+'_all.gb2gt.out','w')

for line in mf:
    a=line.replace("\n","").replace("/","\t").split("\t")
    if a[0] == "CHROM":
        continue
    else:
        chr=a[0]
        pos=a[1]
        b=a[2:]
        bb = []
        for num in b:
            bb.append(int(num))
        ab = []
        for num in bb:
            if(num == -10000 or num == 0):
                continue
            else:
                ab.append(int(num))
        ab2=list(set(ab))
        ab2.sort(reverse=False)
        aa = []
        for i in bb:
            if(i == -10000 or i == 0):
                aa.append(str(i))
            else:
                aa.append(str(int(ab2.index(i))+1))
        allGTOut.write(chr+'\t'+pos+'\t'+ '\t'.join(aa) +'\n') 


allGTOut.close()

