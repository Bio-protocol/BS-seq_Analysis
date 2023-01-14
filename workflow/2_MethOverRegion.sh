#!/bin/bash
##usage: bash file1

SampleX=$1

##Transfer the data format for ViewBS
awk '{if($1!="chr") print $1"\t"$2"\t"$3"\t"$7"\t"$8-$7"\t"$4"\t"$4}' /path/to/aligntoREF/BSMAPratio/${SampleX} > /path/to/aligntoREF/ViewBS/${SampleX}_forViewBS.tab

##ViewBS requires input of compressed file format
bgzip /path/to/aligntoREF/ViewBS/${SampleX}_forViewBS.tab

tabix -C -p vcf /path/to/aligntoREF/ViewBS/${SampleX}_forViewBS.tab.gz

##This step can calculate DNA methylation average level in target region. The value is output in the txt file, and the picture is output in the pdf file.
##The target genes here including total genes in maize.
##Users can modify the figure using value output by 3_Metagene_plot.R
ViewBS MethOverRegion --sample /path/to/aligntoREF/ViewBS/${SampleX}_forViewBS.tab.gz,"CG" --region /path/to/REFfile/target_gene.bed --prefix ${SampleX}_target_gene --context CG --outdir /path/to/aligntoREF/ViewBS/MethOverRegion --binNumber 20 --binLength 100

ViewBS MethOverRegion --sample /path/to/aligntoREF/ViewBS/${SampleX}_forViewBS.tab.gz,"CHG" --region /path/to/REFfile/target_gene.bed --prefix ${SampleX}_target_gene --context CHG --outdir /path/to/aligntoREF/ViewBS/MethOverRegion --binNumber 20 --binLength 100

ViewBS MethOverRegion --sample /path/to/aligntoREF/ViewBS/${SampleX}_forViewBS.tab.gz,"CHH" --region /path/to/REFfile/target_gene.bed --prefix ${SampleX}_target_gene --context CHH --outdir /path/to/aligntoREF/ViewBS/MethOverRegion --binNumber 20 --binLength 100
