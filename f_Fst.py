#!/usr/bin/env python
# encoding: utf-8


import os,sys,string,gzip,argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-i','--input',help="Input vcf file name, genotypes GB, missing value are code as -10000 ",required=False)
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


allFstOut = open(args.out+"_"+popname1+"_"+popname2+'_all.Fst.out','w')

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
                ab1.append(int(num))
                
        ab2 = []
        for num in b2:
            if(num == "-10000"):
                continue
            else:
                ab2.append(int(num))
        ab3=ab1+ab2
        d1=set(ab1)
        d2=set(ab2)
        d3=set(ab3)
        dd1=float(0)
        dd2=float(0)
        dd3=float(0)
        for item in d1:
            dd1=dd1+float(ab1.count(item)*ab1.count(item))/float(len(ab1)*len(ab1))
        
        for item in d2:
            dd2=dd2+float(ab2.count(item)*ab2.count(item))/float(len(ab2)*len(ab2))
            
        for item in d3:
            dd3=dd3+float(ab3.count(item)*ab3.count(item))/float(len(ab3)*len(ab3))
        if dd3 != 1 and len(ab1) >= 1.2*len(ids1) and len(ab2) >= 1.2*len(ids2):
            fst = ((1-dd3)-(((1-dd2)*len(ab2)+(1-dd1)*len(ab1))/(len(ab1)+len(ab2))))/(1-dd3)
            fst = round(fst,6)
        else :
            fst = "inf"
        allFstOut.write(chr+'\t'+pos+'\t'+str(fst)+'\n') 

allFstOut.close()
