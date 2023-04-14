#!/bin/bash
##usage: bash file1 file2

SampleX=$1
Region=$2

##ViewBS requires input of compressed file format
bgzip ${SampleX}_forViewBS.tab

tabix -C -p vcf ${SampleX}_forViewBS.tab.gz

##This step can calculate DNA methylation average level in target region. The value is output in the txt file, and the picture is output in the pdf file.
##The target genes here including total genes in maize.
##Users can modify the figure using value output by 3_Metagene_plot.R
ViewBS MethOverRegion --sample ${SampleX}_forViewBS.tab.gz,"CG" --region ${Region} --prefix ${SampleX}_target_gene --context CG --outdir ./MethOverRegion --binNumber 20 --binLength 100

ViewBS MethOverRegion --sample ${SampleX}_forViewBS.tab.gz,"CHG" --region ${Region} --prefix ${SampleX}_target_gene --context CHG --outdir ./MethOverRegion --binNumber 20 --binLength 100

ViewBS MethOverRegion --sample ${SampleX}_forViewBS.tab.gz,"CHH" --region ${Region} --prefix ${SampleX}_target_gene --context CHH --outdir ./MethOverRegion --binNumber 20 --binLength 100
