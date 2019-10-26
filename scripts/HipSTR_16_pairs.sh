#!/bin/bash
wd=/home/mywd/HipSTR/
cd $wd
HipSTR --bams bam1,bam2,...,bam32 \
--fasta /home/mywd/HipSTR/pip_reference_bed/all_chroms.fa \
--regions /home/mywd/HipSTR/lobSTR_hipSTR.bed \
--str-vcf str_f6_16.vcf.gz \
--log out.log
