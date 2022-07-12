library(ggplot2)

cg<-read.table("6_B73_BSMAP_MethOverRegion_CG.txt",header=T)
chg<-read.table("6_B73_BSMAP_MethOverRegion_CHG.txt",header=T)
chh<-read.table("6_B73_BSMAP_MethOverRegion_CHH.txt",header=T)

data<-rbind(cg,chg,chh)
data$bin_num<-factor(data$bin_num,levels=c('-19','-18','-17','-16','-15','-14','-13','-12','-11','-10','-9','-8','-7','-6','-5','-4','-3','-2','-1','0','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29','30','31','32','33','34','35','36','37','38','39','40'))

png("Metagene_plot.png")
ggplot(data,aes(x=bin_num,y=Methylation_level,fill=sample_name,color=sample_name)) +geom_line(aes(group=sample_name),size=2) + theme(panel.background=element_rect(fill='transparent')) + theme(axis.line = element_line()) + theme(axis.title.x = element_blank(),axis.title.y = element_text(size=18,color='black'),axis.text.y = element_text(size=12,color='black'),axis.text.x = element_text(size=6,angle=45,hjust=1,color='black')) + scale_color_manual(values=c('#6FAC9C','#C78EB9','#73698C')) + theme(strip.text = element_text( size=15),strip.background = element_blank())+ylab("Methylation level")+coord_cartesian(ylim=c(0,1))+geom_vline(xintercept=c(21,40),linetype="dashed",color='grey',size=1)
dev.off()
