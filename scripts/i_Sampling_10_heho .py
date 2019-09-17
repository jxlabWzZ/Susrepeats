#!/usr/bin/env python
# encoding: utf-8

"""
for sample size >=10, random sample 10 individuals for each breed.
then compute the count of homozygotes and heterozygotes in this breed.
repeats for 10 times
"""

import os,sys,string,gzip,argparse
import numpy as np
import random

parser = argparse.ArgumentParser()
parser.add_argument('-i','--input',help="Input vcf file name, GB genotypes",required=False)
parser.add_argument('-o','--out',help='Output file name', required=True)
parser.add_argument('-p1','--pop1',help='population1 file name', required=True)

args = parser.parse_args() 
p1=args.pop1
popname1=p1.split("/")[len(p1.split("/"))-1]



f = open(args.input, 'r')
mf = f.readlines()


Specific = open(args.out+"_"+popname1+'_all.specific.out','w')

for line in mf:
    if line.split("\t")[0] == 'CHROM':
        header = line.replace("\n","").split("\t")
    idnames = header[0:]

for i in [1,2,3,4,5,6,7,8,9,10]:
    if p1 != None:
        ids1 = [i.replace("\n","") for i in open(p1,"r")]
        ids2 = random.sample(ids1,10)
        idx2keep1=[idnames.index(i) for i in ids2]
    indxo=0
    indxe=0
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
                if(num == "m"):
                    continue
                else:
                    ab1.append(int(num))
            d1 = list(set(ab1))
            if len(d1)>=2 :
                indxe += 1
            else :
                indxo += 1
    Specific.write(str(indxe)+'\t'+str(indxo)+'\n')

Specific.close()

