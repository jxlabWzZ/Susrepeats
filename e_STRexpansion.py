#!/usr/bin/env python
# encoding: utf-8


import os,sys,string,gzip,argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-i','--inp',help="Input vcf (copy number genotypes, missing value are code as "-10000") file name",required=False)
parser.add_argument('-o','--out',help='Output file name', required=True)

args = parser.parse_args() 

f = open(args.inp, 'r')

mf = f.readlines()

allexpanOut = open(args.out+'_all.expansion.out','w')

for line in mf:
    a=line.replace("\n","").split("\t")
    chr=a[0].replace("chr","")
    pos=a[1]
    b=a[2:]
    bb = []
    for num in b:
        bb.append(float(num))
    ab = []
    for num in bb:
        if(num == -10000):
            continue
        else:
            ab.append(float(num))
    ab.sort(reverse=False)
    d=len(ab)/2
    e=int(len(ab)*0.95)
    e0=int(len(ab)*0.05)
    median=ab[d]
    ssdd=round(np.std(ab),5)
    allexpanOut.write(chr+'\t'+pos+'\t'+str(ab[0])+'\t'+str(ab[e0])+'\t'+str(median)+'\t'+str(ab[e])+'\t'+str(ab[len(ab)-1])+'\t'+str(ssdd)+'\n') 


allexpanOut.close()

