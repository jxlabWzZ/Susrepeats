#
# pipeline for detect STR from sus scrofa11.1
# TRF + bedtools
#

#download genome
#!/bin/bash
wd=/home/pig_genome11.1/str/
N=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 X Y)
cd $wd
for chr in ${N[@]}; do
echo $chr
wget ftp://ftp.ensembl.org/pub/release-95/fasta/sus_scrofa/dna/Sus_scrofa.Sscrofa11.1.dna.chromosome.${chr}.fa.gz
done

unzip *.fa.gz

#run TRF
#!/bin/bash
wd=/home/pig_genome11.1/str/
trf=/home/biosoft/trf_genome/trf409.legacylinux64
N=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 X Y)
cd $wd
for chr in ${N[@]}; do
echo $chr
${trf} /home/pig_genome11.1/lobGenome/Sus_scrofa.Sscrofa11.1.dna.chromosome.${chr}.fa 2 7 7 80 10 20 100 -d -h
done

# filter step_1 # at least 10bp long; at least 3 repeats
#!/bin/bash
wd=/home/pig_genome11.1/str
cd $wd
N=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 X Y)
for chr in ${N[@]}; do
sed '1,15d' /home/pig_genome11.1/str/Sus_scrofa.Sscrofa11.1.dna.chromosome.${chr}.fa.2.7.7.80.10.20.100.dat | sed 's/ /\t/g' \
| sed 's/^/chr'${chr}'\t&/g' | awk '{if(($4==$6) && (($4==1 && $5>=10) || ($4==2 && $5>=5) || ($4==3 && $5>=4) || ($4>=4 && $5>=3))) print $0}' > Sus_scrofa.Sscrofa11.${chr}.tmp.bed
done

# bedtools merge to get information of overlap or unoverlap,
# to be exact, where is non-overlap region and where is overlap region
#!/bin/bash
wd=/home/pig_genome11.1/str
cd $wd
N=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 X Y)
for chr in ${N[@]}; do
sort -k1,1 -k2,2n Sus_scrofa.Sscrofa11.${chr}.tmp.bed > Sus_scrofa.Sscrofa11.${chr}.tmp2.bed
bedtools merge -i Sus_scrofa.Sscrofa11.${chr}.tmp2.bed -c 1 -o count > Sus_scrofa.Sscrofa11.${chr}.tmp3.bed
done

#overlap info $4>=2; non-overlap info $4==1
cat Sus_scrofa.Sscrofa11.*.tmp3.bed > all.ovl_info.bed
awk '{if($4==1) print $0}' all.ovl_info.bed > all.ovl_non_tmp.bed
awk '{if($4>=2) print $0}' all.ovl_info.bed > all.ovl_is_tmp.bed

#info of repeats, divide into overlap sets, and non-overlap sets
cat Sus_scrofa.Sscrofa11.*.tmp2.bed > all.rep_info.bed
bedtools intersect -a all.rep_info.bed -b all.ovl_non_tmp.bed -wa > all.rep_nonovl_info.bed
bedtools intersect -a all.rep_info.bed -b all.ovl_is_tmp.bed -wa > all.rep_isovl_info.bed

# non_overlap classify ; vntr unit>=7 ,str unit <=6;
awk '{if($4==1) print $0}' all.rep_nonovl_info.bed > all.rep_nonovl_info_str1.bed
awk '{if($4>=2 && $4<=6) print $0}' all.rep_nonovl_info.bed > all.rep_nonovl_info_str26.bed
awk '{if($4>=7) print $0}' all.rep_nonovl_info.bed > all.rep_nonovl_info_vntr7.bed

# for overlaped region, select max trf scores
File=all.ovl_is_tmp.bed
if [ -f $File ]; then
    cat $File | while read line
    do
    echo $line | sed 's/ /\t/g' > m.bed; bedtools intersect -a all.rep_isovl_info.bed -b m.bed -wa | sort -r -nk9 | awk 'NR==1{print $0}' >> my_ov_rep_out
    done
else
    echo "File $File not exist."
fi

# overlap region classify
awk '{if($4==1) print $0}' my_ov_rep_out > all.rep_nonovl_info_str1.ovl.bed
awk '{if($4>=2 && $4<=6) print $0}' my_ov_rep_out > all.rep_nonovl_info_str26.ovl.bed
awk '{if($4>=7) print $0}' my_ov_rep_out > all.rep_nonovl_info_vntr7.ovl.bed

# merge non_overlap and overlap types
mkdir final
cat my_ov_rep_out.chr*.bed > final/my_ov_rep_out.all.bed
awk '{if($4>1 && $4<=6) print $0}' final/my_ov_rep_out.all.bed > final/overlap.26.bed
awk '{if($4==1) print $0}' final/my_ov_rep_out.all.bed > final/overlap.1.bed
awk '{if($4>=7) print $0}' final/my_ov_rep_out.all.bed > final/overlap.vntr.bed
cd final
cat all.rep_nonovl_info_str26.bed overlap.26.bed > total.str.26.bed
cat all.rep_nonovl_info_str1.bed overlap.1.bed > total.str.1.bed
cat all.rep_nonovl_info_vntr7.bed overlap.vntr.bed > total.vntr.bed

# #!/bin/bash
# #PBS -N combine_01
# #PBS -l nodes=1:ppn=1,mem=10gb
# #PBS -e /home/qsublog/
# #PBS -o /home/qsublog/
# #PBS -q cu01
# #PBS -t 0-20
# wd=/home/pig_genome11.1/str
# cd $wd
# # for overlaped region, select max trf scores
# for i in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chrX chrY
# do
# awk -v i=$i '{if($1==i) print $0}' all.ovl_is_tmp.bed > all.ovl_is_tmp.bed.${i}.bed
# done
# N=(`seq 1 18` X Y)
# File=all.ovl_is_tmp.bed.chr${N[$PBS_ARRAYID]}.bed
# if [ -f $File ]; then
# cat $File | while read line
# do
# echo $line | sed 's/ /\t/g' > m.chr${N[$PBS_ARRAYID]}.bed; bedtools intersect -a all.rep_isovl_info.bed -b m.chr${N[$PBS_ARRAYID]}.bed -wa | sort -r -nk9 | awk 'NR==1{print $0}' >> my_ov_rep_out.chr${N[$PBS_ARRAYID]}.bed
# done
# else
# echo "File $File not exist."
# fi
