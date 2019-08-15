###
### STR genotype in population-scale
### lobSTR pipeline
###

### 00. constract lobSTR index files

python2.7 /home/biosoft/lobSTR-bin-Linux-x86_64-3.0.3/scripts/lobstr_index.py \
  --str all.str26.bed \
  --ref Sus_scrofa.Sscrofa11.1.dna.chromosome.all.fa \
  --out /home/pig_genome11/index/ 
 
### 01. constract target_locus_info.tab

python2.7 /home/biosoft/lobSTR-bin-Linux-x86_64-3.0.3/scripts/GetSTRInfo.py all.str26.bed Sus_scrofa.Sscrofa11.1.dna.chromosome.all.fa > info.tab

### 02. bwa mem mapping NGS reads

bwa index Sus_scrofa.Sscrofa11.1.dna.chromosome.all.fa
for i in $(cat sample.id.list)
do
Sample=$i
FLOWCELL=1
LANE=1
RG="@RG\\tID:${Sample}_${FLOWCELL}_${LANE}\\tLB:${Sample}\\tPL:ILLUMINA\\tPU:${FLOWCELL}\\tSM:${Sample}"
bwa mem -R ${RG} -t 30 Sus_scrofa.Sscrofa11.1.dna.chromosome.all.fa $i_R1.fa $i_R2.fa | samtools view -bSu - > $i.bam
samtools sort $i.bam -o $i.sorted.bam
samtools index $i.sorted.bam
done

### 03. To genotype STR using lobSTR

#!/bin/bash
/home/biosoft/lobSTR-bin-Linux-x86_64-3.0.3/bin/allelotype classify \
  --bam S1.sorted.bam,...,Sn.sorted.bam \
  --out STRoutcome \
  --noise_model /home/biosoft/lobSTR-bin-Linux-x86_64-3.0.3/share/lobSTR/models/illumina_v3.pcrfree \
  --strinfo /home/pig_genome11.1/delete_all_index/strinfo.tab  \
  --index-prefix /home/pig_genome11.1/delete_all_index/lobSTR_
  --noweb \
  --filter-mapq0 \
  --filter-clipped \
  --max-repeats-in-ends 3 \
  --min-read-end-match 8
  
### 04. To filter STR outcomes (vcf.files)

#!/bin/bash
for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 X Y
do
python2.7 /home/lobSTR_filter_vcf.py --vcf /home/myproject/STRoutcome.chr${i}.raw.vcf \
  --loc-cov 5 \
  --loc-log-score 0.8 \
  --loc-call-rate 0.6 \
  --loc-max-ref-length 80 \
  --call-cov 3 \
  --call-log-score 0.6 > /home/myproject/filtered_vcf/STRoutcome.${i}.filter.vcf
done

### 05. filtered vcf files for next analysis

#!/bin/bash
cd /home/myproject/filtered_vcf/
for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 X Y
do
awk '{if($7=="PASS") print $0}' /home/myproject/filtered_vcf/STRoutcome.${i}.filter.vcf > /home/myproject/filtered_vcf/STRoutcome.${i}.tmp
done
sed -n '1,30p' /home/myproject/filtered_vcf/STRoutcome.1.tmp > header
cat STRoutcome.*.tmp > chr.all.tmp
cat header chr.all.tmp > STRoutcome.filtered.vcf
rm *.tmp
