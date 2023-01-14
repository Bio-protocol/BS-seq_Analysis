# Bioinformatic Analysis for Bisulfite Sequencing Data

## Introduction
BS-seq_Analysis contains scripts for bisulfite sequencing data analysis including data quality control, alignment, quantification and visualization.

## Workflow
![Abstract1](https://user-images.githubusercontent.com/108569109/212268254-fbf19d05-ee8a-44e2-93fd-67d609ef3e05.png)

## Installation
Users should first install the following software.

1.	Trim Galore (https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) 
2.	BSMAP (Xi and Li, 2009; https://anaconda.org/bioconda/bsmap/files)
3.	Picard (https://broadinstitute.github.io/picard/)
4.	Bamtools (Barnett et al., 2011; https://github.com/pezmaster31/bamtools)
5.	bamUtil (https://github.com/statgen/bamUtil)
6.	BS-SNPer (Gao et al., 2015; https://github.com/hellbelly/BS-Snper)
7.	bedtools (https://bedtools.readthedocs.io/en/latest/)
8.	ViewBS (Huang et al., 2018; https://github.com/xie186/ViewBS)

## Reference data
### Download raw data of B73 inbred line
```
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR850/SRR850328/SRR850328_1.fastq.gz  
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR850/SRR850328/SRR850328_2.fastq.gz
```
### Download raw data of Mo17 inbred line
```
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR850/SRR850332/SRR850332_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR850/SRR850332/SRR850332_2.fastq.gz
```
### Download reference genome sequence
```
wget http://ftp.ensemblgenomes.org/pub/plants/release-45/fasta/zea_mays/dna/Zea_mays.B73_RefGen_v4.dna.toplevel.fa.gz
```

## Major Steps
### Perform quality control analysis for sequencing data
```
trim_galore --phred33 --fastqc --fastqc_args "--noextract --outdir /path/to/rawdata/fastqc/" -o /path/to/rawdata/trimgalore/ --paired /path/to/rawdata/fastq/SampleX_R1.fq.gz /path/to/rawdata/fastq/SampleX_R2.fq.gz
```

### Align reads to reference genome
```
bsmap -a /path/to/rawdata/trimgalore/SampleX_R1_val_1.fq.gz -b /path/to/rawdata/trimgalore/SampleX_R2_val_2.fq.gz -d /path/to/REFfile/REFfile.fa -o /path/to/aligntoREF/BAM/SampleX.bam -v 5 -r 0 -p 8 -q 20 -A AGATCGGAAGAGCGGTTCAGCAGGAATGCCG
```
Adapter sequence according to each sample.

### Remove PCR duplicates
```
java -jar picard.jar SortSam -I /path/to/aligntoREF/BAM/SampleX.bam -O /path/to/aligntoREF/BAM/SampleX_sorted.bam -SORT_ORDER coordinate

java -jar picard.jar MarkDuplicates -I /path/to/aligntoREF/BAM/SampleX_sorted.bam -O /path/to/aligntoREF/BAM/SampleX_sorted_MarkDup.bam -M /path/to/aligntoREF/DuplicateRate/SampleX_MarkDup.txt --REMOVE_DUPLICATES true
```

### Keep properly mapped pairs
```
bamtools filter -isMapped true -isPaired true -isProperPair true -in /path/to/aligntoREF/BAM/SampleX_sorted_MarkDup.bam -out /path/to/aligntoREF/BAM/SampleX_sorted_MarkDup_paired.bam
```

### Clip overlapping read pairs
```
bamUtil/bin/bam clipOverlap --in /path/to/aligntoREF/BAM/SampleX_sorted_MarkDup_paired.bam --out /path/to/aligntoREF/BAM/SampleX_sorted_MarkDup_paired_clipOverlap.bam --stats
```

### Calculate DNA methylation levels in each Cytosines
```
/bsmap-2.90/methratio.py -o /path/to/aligntoREF/BSMAPratio/SampleX -d /path/to/REFfile/REFfile.fa -u -z -r /path/to/aligntoREF/BAM/SampleX_sorted_MarkDup_paired_clipOverlap.bam
```
Python2 and samtools in version with the -X parameter of samtools view are needed.

### Calculate the conversion efficiency of unmethylated cytosine
```
awk -F "\t" '{if($1=="Pt") print}' /path/to/aligntoREF/BSMAPratio/SampleX | awk '{sum1 += $7; sum2 +=$8}END{print sum1"\t"sum2"\t"100-sum1/sum2*100}' > /path/to/aligntoREF/ConversionRate/SampleX_conversion_rate.txt
```

### Detect single nucleotide polymorphisms (SNPs) 
```
perl BS-Snper.pl --fa /path/to/REFfile/REFfile.fa --input /path/to/aligntoREF/BAM/SampleX_sorted_MarkDup_paired_clipOverlap.bam --output /path/to/aligntoREF/BSnper/SampleX_candidate.out --methcg /path/to/aligntoREF/BSnper/SampleX.cg --methchg /path/to/aligntoREF/BSnper/SampleX.chg --methchh /path/to/aligntoREF/BSnper/SampleX.chh --minhetfreq 0.01 --minhomfreq 0.99 --minquali 30 --mincover 4 --maxcover 50 --minread2 2 --errorate 0.02 --mapvalue 30 >/path/to/aligntoREF/BSnper/SampleX_SNP.out 2> /path/to/aligntoREF/BSnper/SampleX_ERR.log
```

### Detect DMRs between different samples
```
bash 1_Detect_DMRs.sh file1 file2
```
file1: The output of methratio.py for sampleX  
file2: The output of methratio.py for another sample  
Users need to pay attention to the file name and path when running the code.

### The visualization of DNA methylation level
```
bash 2_MethOverRegion.sh file1
```
file1: The output of methratio.py for sampleX

The metagene plot drawn in R language.
```
Rscript 3_Metagene_plot.R
```
## Expected results
![Rplot04](https://user-images.githubusercontent.com/108569109/178265151-8ca83c79-3a01-4f68-bc26-e9bc4266ed4d.png)
