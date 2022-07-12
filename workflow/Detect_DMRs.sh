#!/bin/sh

cd /path/to/aligntoREF/DMRs

bedtools makewindows -g /path/to/REFfile/chrom.sizes -w 100 > REFfile_w100.bed

#The next 3 steps are same between samples.
awk 'NR>1{print $1"\t"$2-1"\t"$2"\t"$3"\t"$4"\t"$7"\t"$8}' /path/to/aligntoREF/BSMAPratio/SampleX > SampleX.bed
awk 'NR>1{print $1"\t"$2-1"\t"$2"\t"$3"\t"$4"\t"$7"\t"$8}' /path/to/aligntoREF/BSMAPratio/SampleY > SampleY.bed

bedtools intersect -a REFfile_w100.bed -b SampleX.bed -wa -wb > SampleX_w100.bed
bedtools intersect -a REFfile_w100.bed -b SampleY.bed -wa -wb > SampleY_w100.bed

awk '{n[$1"\t"$2"\t"$3"\t"$8]+=1; sum1[$1"\t"$2"\t"$3"\t"$8]+=$9; sum2[$1"\t"$2"\t"$3"\t"$8]+=$10}END{for(i in n) print i"\t"n[i]"\t"sum1[i]"\t"sum2[i]"\t"sum1[i]/sum2[i]}' SampleX_w100.bed | sort -n -k1 -k2 > SampleX_w100_smry.bed
awk '{n[$1"\t"$2"\t"$3"\t"$8]+=1; sum1[$1"\t"$2"\t"$3"\t"$8]+=$9; sum2[$1"\t"$2"\t"$3"\t"$8]+=$10}END{for(i in n) print i"\t"n[i]"\t"sum1[i]"\t"sum2[i]"\t"sum1[i]/sum2[i]}' SampleY_w100.bed | sort -n -k1 -k2 > SampleY_w100_smry.bed

#Detect the DMRs between SampleX and SampleY.
bedtools intersect -a SampleX_w100_smry.bed -b SampleY_w100_smry.bed -wa -wb > SampleXYcmpr_w100_smry.bed

awk '{if($4=="CG" && $12=="CG" && $5>=3 && $13>=3 && $7>=$5*2 && $15>=$13*2 && (($8-$16)>=0.6 || ($16-$8)>=0.6)) print}' SampleXYcmpr_w100_smry.bed > SampleXYcmpr_w100_smry_CG_DMRs.bed

awk '{if($4=="CHG" && $12=="CHG" && $5>=3 && $13>=3 && $7>=$5*2 && $15>=$13*2 && (($8-$16)>=0.6 || ($16-$8)>=0.6)) print}' SampleXYcmpr_w100_smry.bed > SampleXYcmpr_w100_smry_CHG_DMRs.bed

awk '{if($4=="CHH" && $12=="CHH" && $5>=6 && $13>=6 && $7>=$5*2 && $15>=$13*2 && (($8<0.05 && $16>0.25) || ($8>0.25 && $16<0.05))) print}' SampleXYcmpr_w100_smry.bed > SampleXYcmpr_w100_smry_CHH_DMRs.bed

#Recount the methylation levels in merged DMRs. The final 3 steps for CHG and CHH are same as CG.
bedtools merge -i SampleXYcmpr_w100_smry_CG_DMRs.bed > SampleXYcmpr_w100_smry_CG_DMRs_merged.bed
bedtools merge -i SampleXYcmpr_w100_smry_CHG_DMRs.bed > SampleXYcmpr_w100_smry_CHG_DMRs_merged.bed
bedtools merge -i SampleXYcmpr_w100_smry_CHH_DMRs.bed > SampleXYcmpr_w100_smry_CHH_DMRs_merged.bed

#The final 2 steps are same between samples.
bedtools intersect -a SampleXYcmpr_w100_smry_CG_DMRs_merged.bed -b SampleX.bed -wa -wb > SampleXYcmpr_w100_smry_CG_DMRs_merged_sampleX.bed
bedtools intersect -a SampleXYcmpr_w100_smry_CG_DMRs_merged.bed -b SampleY.bed -wa -wb > SampleXYcmpr_w100_smry_CG_DMRs_merged_sampleY.bed

bedtools intersect -a SampleXYcmpr_w100_smry_CHG_DMRs_merged.bed -b SampleX.bed -wa -wb > SampleXYcmpr_w100_smry_CHG_DMRs_merged_sampleX.bed
bedtools intersect -a SampleXYcmpr_w100_smry_CHG_DMRs_merged.bed -b SampleY.bed -wa -wb > SampleXYcmpr_w100_smry_CHG_DMRs_merged_sampleY.bed

bedtools intersect -a SampleXYcmpr_w100_smry_CHH_DMRs_merged.bed -b SampleX.bed -wa -wb > SampleXYcmpr_w100_smry_CHH_DMRs_merged_sampleX.bed
bedtools intersect -a SampleXYcmpr_w100_smry_CHH_DMRs_merged.bed -b SampleY.bed -wa -wb > SampleXYcmpr_w100_smry_CHH_DMRs_merged_sampleY.bed

awk '{n[$1"\t"$2"\t"$3"\t"$8]+=1; sum1[$1"\t"$2"\t"$3"\t"$8]+=$9; sum2[$1"\t"$2"\t"$3"\t"$8]+=$10}END{for(i in n) print i"\t"n[i]"\t"sum1[i]"\t"sum2[i]"\t"sum1[i]/sum2[i]}' SampleXYcmpr_w100_smry_CG_DMRs_merged_sampleX.bed | sort -n -k1 -k2 > SampleXYcmpr_w100_smry_CG_DMRs_merged_sampleX_smry.bed
awk '{n[$1"\t"$2"\t"$3"\t"$8]+=1; sum1[$1"\t"$2"\t"$3"\t"$8]+=$9; sum2[$1"\t"$2"\t"$3"\t"$8]+=$10}END{for(i in n) print i"\t"n[i]"\t"sum1[i]"\t"sum2[i]"\t"sum1[i]/sum2[i]}' SampleXYcmpr_w100_smry_CG_DMRs_merged_sampleY.bed | sort -n -k1 -k2 > SampleXYcmpr_w100_smry_CG_DMRs_merged_sampleY_smry.bed

awk '{n[$1"\t"$2"\t"$3"\t"$8]+=1; sum1[$1"\t"$2"\t"$3"\t"$8]+=$9; sum2[$1"\t"$2"\t"$3"\t"$8]+=$10}END{for(i in n) print i"\t"n[i]"\t"sum1[i]"\t"sum2[i]"\t"sum1[i]/sum2[i]}' SampleXYcmpr_w100_smry_CHG_DMRs_merged_sampleX.bed | sort -n -k1 -k2 > SampleXYcmpr_w100_smry_CHG_DMRs_merged_sampleX_smry.bed
awk '{n[$1"\t"$2"\t"$3"\t"$8]+=1; sum1[$1"\t"$2"\t"$3"\t"$8]+=$9; sum2[$1"\t"$2"\t"$3"\t"$8]+=$10}END{for(i in n) print i"\t"n[i]"\t"sum1[i]"\t"sum2[i]"\t"sum1[i]/sum2[i]}' SampleXYcmpr_w100_smry_CHG_DMRs_merged_sampleY.bed | sort -n -k1 -k2 > SampleXYcmpr_w100_smry_CHG_DMRs_merged_sampleY_smry.bed

awk '{n[$1"\t"$2"\t"$3"\t"$8]+=1; sum1[$1"\t"$2"\t"$3"\t"$8]+=$9; sum2[$1"\t"$2"\t"$3"\t"$8]+=$10}END{for(i in n) print i"\t"n[i]"\t"sum1[i]"\t"sum2[i]"\t"sum1[i]/sum2[i]}' SampleXYcmpr_w100_smry_CHH_DMRs_merged_sampleX.bed | sort -n -k1 -k2 > SampleXYcmpr_w100_smry_CHH_DMRs_merged_sampleX_smry.bed
awk '{n[$1"\t"$2"\t"$3"\t"$8]+=1; sum1[$1"\t"$2"\t"$3"\t"$8]+=$9; sum2[$1"\t"$2"\t"$3"\t"$8]+=$10}END{for(i in n) print i"\t"n[i]"\t"sum1[i]"\t"sum2[i]"\t"sum1[i]/sum2[i]}' SampleXYcmpr_w100_smry_CHH_DMRs_merged_sampleY.bed | sort -n -k1 -k2 > SampleXYcmpr_w100_smry_CHH_DMRs_merged_sampleY_smry.bed

awk '{if($4=="CG") print}' SampleXYcmpr_w100_smry_CG_DMRs_merged_sampleX_smry.bed > SampleXYcmpr_w100_smry_CG_DMRs_merged_sampleX_smry_CG.bed
awk '{if($4=="CG") print}' SampleXYcmpr_w100_smry_CG_DMRs_merged_sampleY_smry.bed > SampleXYcmpr_w100_smry_CG_DMRs_merged_sampleY_smry_CG.bed

awk '{if($4=="CHG") print}' SampleXYcmpr_w100_smry_CHG_DMRs_merged_sampleX_smry.bed > SampleXYcmpr_w100_smry_CHG_DMRs_merged_sampleX_smry_CHG.bed
awk '{if($4=="CHG") print}' SampleXYcmpr_w100_smry_CHG_DMRs_merged_sampleY_smry.bed > SampleXYcmpr_w100_smry_CHG_DMRs_merged_sampleY_smry_CHG.bed

awk '{if($4=="CHH") print}' SampleXYcmpr_w100_smry_CHH_DMRs_merged_sampleX_smry.bed > SampleXYcmpr_w100_smry_CHH_DMRs_merged_sampleX_smry_CHH.bed
awk '{if($4=="CHH") print}' SampleXYcmpr_w100_smry_CHH_DMRs_merged_sampleY_smry.bed > SampleXYcmpr_w100_smry_CHH_DMRs_merged_sampleY_smry_CHH.bed
