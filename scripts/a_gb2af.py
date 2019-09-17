#!/usr/bin/env python
# encoding: utf-8

"""
from GB file estimate alle count, alle frequency, heterozygous ,include minor(second major) allele frequency
"""

import os,sys,string,gzip,argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-i','--inp',help="Input vcf file name,GB genotypes, missing value is code -10000",required=False)
parser.add_argument('-o','--out',help='Output file name', required=True)

args = parser.parse_args() 

f = open(args.inp, 'r')

mf = f.readlines()

allAFOut = open(args.out+'_all.gb2.af.he.out','w')

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
            bb.append(num)
        ab = []
        for num in bb:
            if(num == "-10000"):
                continue
            else:
                ab.append(int(num))
        total_count = len(ab)
        ab2=list(set(ab))
        ab2.sort(reverse=False)
        alle_count = len(ab2)
        aa = []
        aa2 = []
        Hexp_0 = 0
        for i in ab2:
            aa.append(str(i))
            aa2.append(str(round(float(ab.count(i))/float(total_count),6)))
            aa.append(str(round(float(ab.count(i))/float(total_count),6)))
            Hexp_0 += round((float(ab.count(i))/float(total_count))*(float(ab.count(i))/float(total_count)),6)
        Hexp = 1-Hexp_0
        aa2.sort(reverse=False)
        min = aa2[0]
        max = aa2[-1]
        if alle_count >= 2 :
            sec = aa2[-2]
        else:
            sec = "mono"
        allAFOut.write(chr+'\t'+pos+'\t'+str(min)+'\t'+str(sec)+'\t'+str(max)+'\t'+str(Hexp)+'\t'+str(alle_count)+'\t'+str(total_count)+'\t'+ '\t'.join(aa) +'\n') 

allAFOut.close()
