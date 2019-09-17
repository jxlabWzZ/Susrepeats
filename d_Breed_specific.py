#!/usr/bin/env python
# encoding: utf-8

"""
breed_specific
"""

import os,sys,string,gzip,argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-i','--input',help="Input vcf file name,GB genotypes,missing value is -10000",required=False)
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

if p1 != None:
    ids1 = [i.replace("\n","") for i in open(p1,"r")]  
    idx2keep1=[idnames.index(i) for i in ids1]

idname2 = header[2:]
ids2 = [x for x in idname2 if x not in ids1]
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
        ab1 = []
        for num in b1:
            if(num == "-10000"):
                continue
            else:
                ab1.append(int(num))
        a2=[]
        for i in idx2keep2:
            a2.append(a[2*i-2])
            a2.append(a[2*i-1])
        b2=a2
        ab2 = []
        for num in b2:
            if(num == "-10000"):
                continue
            else:
                ab2.append(int(num))
        b3=a[2:]
        ab3 = []
        for num in b3:
            if(num == "-10000"):
                continue
            else:
                ab3.append(int(num))
        total=len(ab3)
        total0=len(ab1)
        d1 = list(set(ab1))
        d2 = list(set(ab2))
        d3 = list(set(ab3))
        for num in d3:
            if num not in d1:
                itmp = "0"
                itmp2 = "0"
                he1 = "0"
                he0 = "0"
                Specific.write(chr+'\t'+pos+'\t'+str(num)+'\t'+str(itmp)+'\t'+str(itmp2)+'\t'+str(he1)+'\t'+str(he0)+'\n') 
            else:
                if (num in d2) and (num in d1):
                    itmp = "1"
                    itmp2 = "0"
                    he1 = float(ab1.count(num)/float(total))
                    he0 = float(ab1.count(num)/float(total0))
                    he1 = round(he1,6)
                    he0 = round(he0,6)
                    Specific.write(chr+'\t'+pos+'\t'+str(num)+'\t'+str(itmp)+'\t'+str(itmp2)+'\t'+str(he1)+'\t'+str(he0)+'\t'+str(ab1.count(num))+'\n')
                else :
                    itmp = "1"
                    itmp2 = "1"
                    he1 = float(ab1.count(num)/float(total))
                    he0 = float(ab1.count(num)/float(total0))
                    he1 = round(he1,6)
                    he0 = round(he0,6)
                    Specific.write(chr+'\t'+pos+'\t'+str(num)+'\t'+str(itmp)+'\t'+str(itmp2)+'\t'+str(he1)+'\t'+str(he0)+'\t'+str(ab1.count(num))+'\n')

Specific.close()

