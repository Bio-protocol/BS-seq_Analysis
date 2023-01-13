#!/bin/sh
##usage: bash file1 file2

SampleX=$1
SampleY=$2

cd /path/to/aligntoREF/DMRs

bedtools makewindows -g /path/to/REFfile/chrom.sizes -w 100 > REFfile_w100.bed

#The next 3 steps are same between samples.
awk 'NR>1{print $1"\t"$2-1"\t"$2"\t"$3"\t"$4"\t"$7"\t"$8}' /path/to/aligntoREF/BSMAPratio/${SampleX} > ${SampleX}.bed
awk 'NR>1{print $1"\t"$2-1"\t"$2"\t"$3"\t"$4"\t"$7"\t"$8}' /path/to/aligntoREF/BSMAPratio/${SampleY} > ${SampleY}.bed

bedtools intersect -a REFfile_w100.bed -b ${SampleX}.bed -wa -wb > ${SampleX}_w100.bed
bedtools intersect -a REFfile_w100.bed -b ${SampleY}.bed -wa -wb > ${SampleY}_w100.bed

awk '{n[$1"\t"$2"\t"$3"\t"$8]+=1; sum1[$1"\t"$2"\t"$3"\t"$8]+=$9; sum2[$1"\t"$2"\t"$3"\t"$8]+=$10}END{for(i in n) print i"\t"n[i]"\t"sum1[i]"\t"sum2[i]"\t"sum1[i]/sum2[i]}' ${SampleX}_w100.bed | sort -n -k1 -k2 > ${SampleX}_w100_smry.bed
awk '{n[$1"\t"$2"\t"$3"\t"$8]+=1; sum1[$1"\t"$2"\t"$3"\t"$8]+=$9; sum2[$1"\t"$2"\t"$3"\t"$8]+=$10}END{for(i in n) print i"\t"n[i]"\t"sum1[i]"\t"sum2[i]"\t"sum1[i]/sum2[i]}' ${SampleY}_w100.bed | sort -n -k1 -k2 > ${SampleY}_w100_smry.bed

#Detect the DMRs between SampleX and SampleY.
bedtools intersect -a ${SampleX}_w100_smry.bed -b ${SampleY}_w100_smry.bed -wa -wb > ${SampleX}_${SampleY}_cmpr_w100_smry.bed

awk '{if($4=="CG" && $12=="CG" && $5>=3 && $13>=3 && $7>=$5*2 && $15>=$13*2 && (($8-$16)>=0.6 || ($16-$8)>=0.6)) print}' ${SampleX}_${SampleY}_cmpr_w100_smry.bed > ${SampleX}_${SampleY}_cmpr_w100_smry_CG_DMRs.bed

awk '{if($4=="CHG" && $12=="CHG" && $5>=3 && $13>=3 && $7>=$5*2 && $15>=$13*2 && (($8-$16)>=0.6 || ($16-$8)>=0.6)) print}' ${SampleX}_${SampleY}_cmpr_w100_smry.bed > ${SampleX}_${SampleY}_cmpr_w100_smry_CHG_DMRs.bed

awk '{if($4=="CHH" && $12=="CHH" && $5>=6 && $13>=6 && $7>=$5*2 && $15>=$13*2 && (($8<0.05 && $16>0.25) || ($8>0.25 && $16<0.05))) print}' ${SampleX}_${SampleY}_cmpr_w100_smry.bed > ${SampleX}_${SampleY}_cmpr_w100_smry_CHH_DMRs.bed

#Recount the methylation levels in merged DMRs. The final 3 steps for CHG and CHH are same as CG.
bedtools merge -i ${SampleX}_${SampleY}_cmpr_w100_smry_CG_DMRs.bed > ${SampleX}_${SampleY}_cmpr_w100_smry_CG_DMRs_merged.bed
bedtools merge -i ${SampleX}_${SampleY}_cmpr_w100_smry_CHG_DMRs.bed > ${SampleX}_${SampleY}_cmpr_w100_smry_CHG_DMRs_merged.bed
bedtools merge -i ${SampleX}_${SampleY}_cmpr_w100_smry_CHH_DMRs.bed > ${SampleX}_${SampleY}_cmpr_w100_smry_CHH_DMRs_merged.bed

#The final 2 steps are same between samples.
bedtools intersect -a ${SampleX}_${SampleY}_cmpr_w100_smry_CG_DMRs_merged.bed -b ${SampleX}.bed -wa -wb > ${SampleX}_${SampleY}_cmpr_w100_smry_CG_DMRs_merged_${SampleX}.bed
bedtools intersect -a ${SampleX}_${SampleY}_cmpr_w100_smry_CG_DMRs_merged.bed -b ${SampleY}.bed -wa -wb > ${SampleX}_${SampleY}_cmpr_w100_smry_CG_DMRs_merged_${SampleY}.bed

bedtools intersect -a ${SampleX}_${SampleY}_cmpr_w100_smry_CHG_DMRs_merged.bed -b ${SampleX}.bed -wa -wb > ${SampleX}_${SampleY}_cmpr_w100_smry_CHG_DMRs_merged_${SampleX}.bed
bedtools intersect -a ${SampleX}_${SampleY}_cmpr_w100_smry_CHG_DMRs_merged.bed -b ${SampleY}.bed -wa -wb > ${SampleX}_${SampleY}_cmpr_w100_smry_CHG_DMRs_merged_${SampleY}.bed

bedtools intersect -a ${SampleX}_${SampleY}_cmpr_w100_smry_CHH_DMRs_merged.bed -b ${SampleX}.bed -wa -wb > ${SampleX}_${SampleY}_cmpr_w100_smry_CHH_DMRs_merged_${SampleX}.bed
bedtools intersect -a ${SampleX}_${SampleY}_cmpr_w100_smry_CHH_DMRs_merged.bed -b ${SampleY}.bed -wa -wb > ${SampleX}_${SampleY}_cmpr_w100_smry_CHH_DMRs_merged_${SampleY}.bed

awk '{n[$1"\t"$2"\t"$3"\t"$8]+=1; sum1[$1"\t"$2"\t"$3"\t"$8]+=$9; sum2[$1"\t"$2"\t"$3"\t"$8]+=$10}END{for(i in n) print i"\t"n[i]"\t"sum1[i]"\t"sum2[i]"\t"sum1[i]/sum2[i]}' ${SampleX}_${SampleY}_cmpr_w100_smry_CG_DMRs_merged_${SampleX}.bed | sort -n -k1 -k2 > ${SampleX}_${SampleY}_cmpr_w100_smry_CG_DMRs_merged_${SampleX}_smry.bed
awk '{n[$1"\t"$2"\t"$3"\t"$8]+=1; sum1[$1"\t"$2"\t"$3"\t"$8]+=$9; sum2[$1"\t"$2"\t"$3"\t"$8]+=$10}END{for(i in n) print i"\t"n[i]"\t"sum1[i]"\t"sum2[i]"\t"sum1[i]/sum2[i]}' ${SampleX}_${SampleY}_cmpr_w100_smry_CG_DMRs_merged_${SampleY}.bed | sort -n -k1 -k2 > ${SampleX}_${SampleY}_cmpr_w100_smry_CG_DMRs_merged_${SampleY}_smry.bed

awk '{n[$1"\t"$2"\t"$3"\t"$8]+=1; sum1[$1"\t"$2"\t"$3"\t"$8]+=$9; sum2[$1"\t"$2"\t"$3"\t"$8]+=$10}END{for(i in n) print i"\t"n[i]"\t"sum1[i]"\t"sum2[i]"\t"sum1[i]/sum2[i]}' ${SampleX}_${SampleY}_cmpr_w100_smry_CHG_DMRs_merged_${SampleX}.bed | sort -n -k1 -k2 > ${SampleX}_${SampleY}_cmpr_w100_smry_CHG_DMRs_merged_${SampleX}_smry.bed
awk '{n[$1"\t"$2"\t"$3"\t"$8]+=1; sum1[$1"\t"$2"\t"$3"\t"$8]+=$9; sum2[$1"\t"$2"\t"$3"\t"$8]+=$10}END{for(i in n) print i"\t"n[i]"\t"sum1[i]"\t"sum2[i]"\t"sum1[i]/sum2[i]}' ${SampleX}_${SampleY}_cmpr_w100_smry_CHG_DMRs_merged_${SampleY}.bed | sort -n -k1 -k2 > ${SampleX}_${SampleY}_cmpr_w100_smry_CHG_DMRs_merged_${SampleY}_smry.bed

awk '{n[$1"\t"$2"\t"$3"\t"$8]+=1; sum1[$1"\t"$2"\t"$3"\t"$8]+=$9; sum2[$1"\t"$2"\t"$3"\t"$8]+=$10}END{for(i in n) print i"\t"n[i]"\t"sum1[i]"\t"sum2[i]"\t"sum1[i]/sum2[i]}' ${SampleX}_${SampleY}_cmpr_w100_smry_CHH_DMRs_merged_${SampleX}.bed | sort -n -k1 -k2 > ${SampleX}_${SampleY}_cmpr_w100_smry_CHH_DMRs_merged_${SampleX}_smry.bed
awk '{n[$1"\t"$2"\t"$3"\t"$8]+=1; sum1[$1"\t"$2"\t"$3"\t"$8]+=$9; sum2[$1"\t"$2"\t"$3"\t"$8]+=$10}END{for(i in n) print i"\t"n[i]"\t"sum1[i]"\t"sum2[i]"\t"sum1[i]/sum2[i]}' ${SampleX}_${SampleY}_cmpr_w100_smry_CHH_DMRs_merged_${SampleY}.bed | sort -n -k1 -k2 > ${SampleX}_${SampleY}_cmpr_w100_smry_CHH_DMRs_merged_${SampleY}_smry.bed

awk '{if($4=="CG") print}' ${SampleX}_${SampleY}_cmpr_w100_smry_CG_DMRs_merged_${SampleX}_smry.bed > ${SampleX}_${SampleY}_cmpr_w100_smry_CG_DMRs_merged_${SampleX}_smry_CG.bed
awk '{if($4=="CG") print}' ${SampleX}_${SampleY}_cmpr_w100_smry_CG_DMRs_merged_${SampleY}_smry.bed > ${SampleX}_${SampleY}_cmpr_w100_smry_CG_DMRs_merged_${SampleY}_smry_CG.bed

awk '{if($4=="CHG") print}' ${SampleX}_${SampleY}_cmpr_w100_smry_CHG_DMRs_merged_${SampleX}_smry.bed > ${SampleX}_${SampleY}_cmpr_w100_smry_CHG_DMRs_merged_${SampleX}_smry_CHG.bed
awk '{if($4=="CHG") print}' ${SampleX}_${SampleY}_cmpr_w100_smry_CHG_DMRs_merged_${SampleY}_smry.bed > ${SampleX}_${SampleY}_cmpr_w100_smry_CHG_DMRs_merged_${SampleY}_smry_CHG.bed

awk '{if($4=="CHH") print}' ${SampleX}_${SampleY}_cmpr_w100_smry_CHH_DMRs_merged_${SampleX}_smry.bed > ${SampleX}_${SampleY}_cmpr_w100_smry_CHH_DMRs_merged_${SampleX}_smry_CHH.bed
awk '{if($4=="CHH") print}' ${SampleX}_${SampleY}_cmpr_w100_smry_CHH_DMRs_merged_${SampleY}_smry.bed > ${SampleX}_${SampleY}_cmpr_w100_smry_CHH_DMRs_merged_${SampleY}_smry_CHH.bed
