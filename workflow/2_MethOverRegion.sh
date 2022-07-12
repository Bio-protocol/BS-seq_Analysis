#!/bin/bash

awk '{if($1!="chr") print $1"\t"$2"\t"$3"\t"$7"\t"$8-$7"\t"$4"\t"$4}' /path/to/aligntoREF/BSMAPratio/SampleX > /path/to/aligntoREF/ViewBS/SampleX_forViewBS.tab

bgzip /path/to/aligntoREF/ViewBS/SampleX_forViewBS.tab

tabix -C -p vcf /path/to/aligntoREF/ViewBS/SampleX_forViewBS.tab.gz

ViewBS MethOverRegion --sample /path/to/aligntoREF/ViewBS/SampleX_forViewBS.tab.gz,"SampleX_CG" --region /path/to/REFfile/target_gene.bed --prefix SampleX_target_gene --context CG --outdir /path/to/aligntoREF/ViewBS/MethOverRegion --binNumber 20 --binLength 100

ViewBS MethOverRegion --sample /path/to/aligntoREF/ViewBS/SampleX_forViewBS.tab.gz,"SampleX_CHG" --region /path/to/REFfile/target_gene.bed --prefix SampleX_target_gene --context CHG --outdir /path/to/aligntoREF/ViewBS/MethOverRegion --binNumber 20 --binLength 100

ViewBS MethOverRegion --sample /path/to/aligntoREF/ViewBS/SampleX_forViewBS.tab.gz,"SampleX_CHH" --region /path/to/REFfile/target_gene.bed --prefix SampleX_target_gene --context CHH --outdir /path/to/aligntoREF/ViewBS/MethOverRegion --binNumber 20 --binLength 100
