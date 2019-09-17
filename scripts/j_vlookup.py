#!/usr/bin/env python
# encoding: utf-8
"""
To get reverse complementary sequences and cyclic sequences
for instance:
ACG -> ACG CGA GAC CGT TCG GTC
"""

import os,sys,string,gzip,argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i','--inp',help="Input vcf file name",required=False)
parser.add_argument('-o','--out',help='Output file name', required=True)

args = parser.parse_args() 

f = str(args.inp)

allFstOut = open(args.inp+'.'+args.out,'w')

def demo(seq):
    for k in xrange(len(seq)):
        seqk = seq[k:]+seq[:k]
        allFstOut.write(seqk+'\n')


def rev_c(seq):
    rev     = seq[::-1]
    mapping = {"A":"T", "T":"A", "C":"G", "G":"C", "N":"N"}
    res     = ""
    for i in xrange(len(rev)):
        res = res+mapping[rev[i]]
    demo(res)
    demo(seq)

rev_c(f)

