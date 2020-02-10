#! /usr/bin/Rscript --no-save --no-restore

.libPaths("/home/jingwang/R/x86_64-redhat-linux-gnu-library/3.2")
library(utils)
library(lattice)
library(RColorBrewer)
library(ggplot2)

setwd("/proj/b2011141/nobackup/eQTL_paper/gene_list/Potra_Potri/gkaks/core_eGenes")

dnds=read.table("gene.dnds.txt",header=F)
dnds$V4=gsub("\\.[0-9]","",dnds$V2)


###eGenes_core
eGene_core=read.table("eGene_core.txt")
eGene_core_dnds=dnds[which(dnds$V4 %in% eGene_core$V1),c("V4","V3")]
names(eGene_core_dnds)=c("gene","kaks")
eGene_core_dnds$eGene="eGene"
eGene_core_dnds$core="core"

###eGenes_noncore
eGene_noncore=read.table("eGene_noncore.txt")
eGene_noncore_dnds=dnds[which(dnds$V4 %in% eGene_noncore$V1),c("V4","V3")]
names(eGene_noncore_dnds)=c("gene","kaks")
eGene_noncore_dnds$eGene="eGene"
eGene_noncore_dnds$core="Non-core"

###Non_eGene_core
Non_eGene_core=read.table("non_eGene_core.txt")
Non_eGene_core_dnds=dnds[which(dnds$V4 %in% Non_eGene_core$V1),c("V4","V3")]
names(Non_eGene_core_dnds)=c("gene","kaks")
Non_eGene_core_dnds$eGene="Non_eGene"
Non_eGene_core_dnds$core="core"

###Non_eGene_noncore
Non_eGene_noncore=read.table("non_eGene_noncore.txt")
Non_eGene_noncore_dnds=dnds[which(dnds$V4 %in% Non_eGene_noncore$V1),c("V4","V3")]
names(Non_eGene_noncore_dnds)=c("gene","kaks")
Non_eGene_noncore_dnds$eGene="Non_eGene"
Non_eGene_noncore_dnds$core="Non-core"

total_dnds=rbind(eGene_core_dnds,eGene_noncore_dnds,Non_eGene_core_dnds,Non_eGene_noncore_dnds)


###make the plot
png(filename="eGene_core.dnds.png",width=5,height=5,units='in',res=300)
a=ggplot(total_dnds)+geom_boxplot(aes(x=core,y=kaks,fill=eGene),outlier.size=NA)+labs(y=expression(d[N]/d[S]))+scale_y_continuous(limits=c(0,1.75))

###notify the significance for different paris of groups
###No-outliers
###Mark the significance
##0.01<P<0.05 *
##0.001<P<0.01 **
##P<0.001  ***

#> nrow(eGene_core_dnds)
#[1] 133
#> nrow(eGene_noncore_dnds)
#[1] 5758
#> nrow(Non_eGene_noncore_dnds)
#[1] 12828
#> nrow(Non_eGene_core_dnds)
#[1] 1599

wilcox.test(eGene_core_dnds$kaks,Non_eGene_core_dnds$kaks)  ###p-value = 0.8271
wilcox.test(eGene_noncore_dnds$kaks,Non_eGene_noncore_dnds$kaks)   ###p-value = 6.425e-05
wilcox.test(eGene_core_dnds$kaks,eGene_noncore_dnds$kaks)  ###p-value = 0.003441
wilcox.test(Non_eGene_core_dnds$kaks,Non_eGene_noncore_dnds$kaks)   ###p-value = 6.986e-13

df1 <- data.frame(a = c(0.8, 0.8,1.2,1.2), b = c(1.3, 1.35, 1.35, 1.3))
df2 <- data.frame(a = c(1.8, 1.8,2.2,2.2), b = c(1.3, 1.35, 1.35, 1.3))
df3 <- data.frame(a = c(0.8,0.8,1.8,1.8), b = c(1.45, 1.5, 1.5, 1.45))
df4 <- data.frame(a = c(1.2,1.2,2.2,2.2), b = c(1.6, 1.65, 1.65, 1.6))

a+geom_line(data = df1, aes(x = a, y = b)) + annotate("text", x = 1, y = 1.4, label = "n.s.", size = 6)+geom_line(data = df2, aes(x = a, y = b)) + annotate("text", x = 2, y = 1.4, label = "***", size = 6)+geom_line(data = df3, aes(x = a, y = b)) + annotate("text", x = 1.3, y = 1.55, label = "**", size = 6)+geom_line(data = df4, aes(x = a, y = b)) + annotate("text", x = 1.7, y = 1.7, label = "***", size = 6)

dev.off()








