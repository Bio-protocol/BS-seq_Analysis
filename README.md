# Bioinformatic Analysis for Bisulfite Sequencing Data

## Introduction
BSseq-scripts contains scripts for bisulfite sequencing data analysis including data quality control, alignment, quantification and visualization.

## Workflow
![image](https://user-images.githubusercontent.com/108569109/178256142-6a255767-24f3-4341-b97e-bced6e03264a.png)

## Installation
Users should first install the following software.

1.	Trim Galore (https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) 
2.	BSMAP (Xi and Li, 2009; https://code.google.com/archive/p/bsmap/)
3.	Picard (https://broadinstitute.github.io/picard/)
4.	Bamtools (Barnett et al., 2011; https://github.com/pezmaster31/bamtools)
5.	bamUtil (https://github.com/statgen/bamUtil)
6.	BS-SNPer (Gao et al., 2015; https://github.com/hellbelly/BS-Snper)
7.	bedtools (https://bedtools.readthedocs.io/en/latest/)
8.	ViewBS (Huang et al., 2018; https://github.com/xie186/ViewBS)

## Input data
```
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR850/SRR850328/SRR850328_1.fastq.gz  
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR850/SRR850328/SRR850328_2.fastq.gz  
wget http://ftp.ensemblgenomes.org/pub/plants/release-45/fasta/zea_mays/dna/Zea_mays.B73_RefGen_v4.dna.toplevel.fa.gz
```

## Major Steps
### Perform quality control analysis for sequencing data
`trim_galore --phred33 --fastqc --fastqc_args "--noextract --outdir /path/to/rawdata/fastqc/" -o /path/to/rawdata/trimgalore/ --paired /path/to/rawdata/fastq/SampleX_R1.fq.gz /path/to/rawdata/fastq/SampleX_R2.fq.gz`

### Align reads to reference genome
`bsmap -a /path/to/rawdata/trimgalore/SampleX_R1_val_1.fq.gz -b /path/to/rawdata/trimgalore/SampleX_R2_val_2.fq.gz -d /path/to/REFfile/REFfile.fa -o /path/to/aligntoREF/BAM/SampleX.bam -v 5 -r 0 -p 8 -q 20 -A AGATCGGAAGAGCGGTTCAGCAGGAATGCCG`(adapter sequence according to each sample)

### Remove PCR duplicates
`java -jar picard.jar SortSam -I /path/to/aligntoREF/BAM/SampleX.bam -O /path/to/aligntoREF/BAM/SampleX_sorted.bam -SORT_ORDER coordinate`

`java -jar picard.jar MarkDuplicates -I /path/to/aligntoREF/BAM/SampleX_sorted.bam -O /path/to/aligntoREF/BAM/SampleX_sorted_MarkDup.bam -M /path/to/aligntoREF/DuplicateRate/SampleX_MarkDup.txt --REMOVE_DUPLICATES true`

### Keep properly mapped pairs
`bamtools filter -isMapped true -isPaired true -isProperPair true -in /path/to/aligntoREF/BAM/SampleX_sorted_MarkDup.bam -out /path/to/aligntoREF/BAM/SampleX_sorted_MarkDup_paired.bam`

### Clip overlapping read pairs
`bamUtil/bin/bam clipOverlap --in /path/to/aligntoREF/BAM/SampleX_sorted_MarkDup_paired.bam --out /path/to/aligntoREF/BAM/SampleX_sorted_MarkDup_paired_clipOverlap.bam --stats`

### Calculate DNA methylation levels in each Cytosines
`/bsmap-2.90/methratio.py -o /path/to/aligntoREF/BSMAPratio/SampleX -d /path/to/REFfile/REFfile.fa -u -z -r /path/to/aligntoREF/BAM/SampleX_sorted_MarkDup_paired_clipOverlap.bam`

### Calculate the conversion efficiency of unmethylated cytosine
`awk -F "\t" '{if($1=="Pt") print}' /path/to/aligntoREF/BSMAPratio/SampleX | awk '{sum1 += $7; sum2 +=$8}END{print sum1"\t"sum2"\t"100-sum1/sum2\*100}' > /path/to/aligntoREF/ConversionRate/SampleX_conversion_rate.txt`

### Detect single nucleotide polymorphisms (SNPs) 
`perl BS-Snper.pl --fa /path/to/REFfile/REFfile.fa --input /path/to/aligntoREF/BAM/SampleX_sorted_MarkDup_paired_clipOverlap.bam --output /path/to/aligntoREF/BSnper/SampleX_candidate.out --methcg /path/to/aligntoREF/BSnper/SampleX.cg --methchg /path/to/aligntoREF/BSnper/SampleX.chg --methchh /path/to/aligntoREF/BSnper/SampleX.chh --minhetfreq 0.01 --minhomfreq 0.99 --minquali 30 --mincover 4 --maxcover 50 --minread2 2 --errorate 0.02 --mapvalue 30 >/path/to/aligntoREF/BSnper/SampleX_SNP.out 2> /path/to/aligntoREF/BSnper/SampleX_ERR.log`

### Detect DMRs between different samples
`cd /path/to/aligntoREF/DMRs`

`bedtools makewindows -g /path/to/REFfile/chrom.sizes -w 100 > REFfile_w100.bed`

#The next 3 steps are same between samples.  
`awk 'NR>1{print $1"\t"$2-1"\t"$2"\t"$3"\t"$4"\t"$7"\t"$8}' /path/to/aligntoREF/BSMAPratio/SampleX > SampleX.bed`

`bedtools intersect -a REFfile_w100.bed -b SampleX.bed -wa -wb > SampleX_w100.bed`

`awk '{n[$1"\t"$2"\t"$3"\t"$8]+=1; sum1[$1"\t"$2"\t"$3"\t"$8]+=$9; sum2[$1"\t"$2"\t"$3"\t"$8]+=$10}END{for(i in n) print i"\t"n[i]"\t"sum1[i]"\t"sum2[i]"\t"sum1[i]/sum2[i]}' SampleX_w100.bed | sort -n -k1 -k2 > SampleX_w100_smry.bed`

#Detect the DMRs between SampleX and SampleY.  
`bedtools intersect -a SampleX_w100_smry.bed -b SampleY_w100_smry.bed -wa -wb > SampleXYcmpr_w100_smry.bed`

`awk '{if($4=="CG" && $12=="CG" && $5>=3 && $13>=3 && $7>=$5*2 && $15>=$13*2 && (($8-$16)>=0.6 || ($16-$8)>=0.6)) print}' SampleXYcmpr_w100_smry.bed > SampleXYcmpr_w100_smry_CG_DMRs.bed`

`awk '{if($4=="CHG" && $12=="CHG" && $5>=3 && $13>=3 && $7>=$5*2 && $15>=$13*2 && (($8-$16)>=0.6 || ($16-$8)>=0.6)) print}' SampleXYcmpr_w100_smry.bed > SampleXYcmpr_w100_smry_CHG_DMRs.bed`

`awk '{if($4=="CHH" && $12=="CHH" && $5>=6 && $13>=6 && $7>=$5*2 && $15>=$13*2 && (($8<0.05 && $16>0.25) || ($8>0.25 && $16<0.05))) print}' SampleXYcmpr_w100_smry.bed > SampleXYcmpr_w100_smry_CHH_DMRs.bed`

#Recount the methylation levels in merged DMRs. The final 3 steps for CHG and CHH are same as CG.  
`bedtools merge -i SampleXYcmpr_w100_smry_CG_DMRs.bed > SampleXYcmpr_w100_smry_CG_DMRs_merged.bed`

#The final 2 steps are same between samples.  
`bedtools intersect -a SampleXYcmpr_w100_smry_CG_DMRs_merged.bed -b SampleX.bed -wa -wb > SampleXYcmpr_w100_smry_CG_DMRs_merged_sampleX.bed`

`awk '{n[$1"\t"$2"\t"$3"\t"$8]+=1; sum1[$1"\t"$2"\t"$3"\t"$8]+=$9; sum2[$1"\t"$2"\t"$3"\t"$8]+=$10}END{for(i in n) print i"\t"n[i]"\t"sum1[i]"\t"sum2[i]"\t"sum1[i]/sum2[i]}' SampleXYcmpr_w100_smry_CG_DMRs_merged_sampleX.bed | sort -n -k1 -k2 > SampleXYcmpr_w100_smry_CG_DMRs_merged_sampleX_smry.bed`

`awk '{if($4=="CG") print}' SampleXYcmpr_w100_smry_CG_DMRs_merged_sampleX_smry.bed > SampleXYcmpr_w100_smry_CG_DMRs_merged_sampleX_smry_CG.bed`

### The visualization of DNA methylation level
`awk '{if($1!="chr") print $1"\t"$2"\t"$3"\t"$7"\t"$8-$7"\t"$4"\t"$4}' /path/to/aligntoREF/BSMAPratio/SampleX > /path/to/aligntoREF/ViewBS/SampleX_forViewBS.tab`

`bgzip /path/to/aligntoREF/ViewBS/SampleX_forViewBS.tab`

`tabix -C -p vcf /path/to/aligntoREF/ViewBS/SampleX_forViewBS.tab.gz`

`ViewBS MethOverRegion --sample /path/to/aligntoREF/ViewBS/SampleX_forViewBS.tab.gz,"SampleX_CG" --region /path/to/REFfile/target_gene.bed --prefix SampleX_target_gene --context CG --outdir /path/to/aligntoREF/ViewBS/MethOverRegion --binNumber 20 --binLength 100`

`ViewBS MethOverRegion --sample /path/to/aligntoREF/ViewBS/SampleX_forViewBS.tab.gz,"SampleX_CHG" --region /path/to/REFfile/target_gene.bed --prefix SampleX_target_gene --context CHG --outdir /path/to/aligntoREF/ViewBS/MethOverRegion --binNumber 20 --binLength 100`

`ViewBS MethOverRegion --sample /path/to/aligntoREF/ViewBS/SampleX_forViewBS.tab.gz,"SampleX_CHH" --region /path/to/REFfile/target_gene.bed --prefix SampleX_target_gene --context CHH --outdir /path/to/aligntoREF/ViewBS/MethOverRegion --binNumber 20 --binLength 100`

## Expected results
![Rplot04](https://user-images.githubusercontent.com/108569109/178265151-8ca83c79-3a01-4f68-bc26-e9bc4266ed4d.png)
