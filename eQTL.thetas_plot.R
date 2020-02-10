#! /usr/bin/Rscript --no-save --no-restore

# File: plotThetas.R, the following is one example of how to use this script
#Rscript angsd_plot_10kbthetas.R /proj/b2011141/nobackup/tremula_vs_tremuloides_paper/bwa_mem_alignment/ANGSD/unfolded_SFS/tremula/tremula_19/tremula_19.thetas10kbwindow.gz.pestPG /proj/b2011141/nobackup/tremula_vs_tremuloides_paper/bwa_mem_alignment/ANGSD/unfolded_SFS/tremuloides/tremuloides_19/tremuloides_19.thetas10kbwindow.gz.pestPG /proj/b2011141/nobackup/tremula_vs_tremuloides_paper/bwa_mem_alignment/ANGSD/summary/Fst/Chr19/tremula_tremuloides.chr19.win.fst.statwindow10000_step10000.fst.txt /proj/b2011141/nobackup/tremula_vs_tremuloides_paper/bwa_mem_alignment/ANGSD/summary/Fst/Chr19/Chr19.tremula_fixed.fst.txt /proj/b2011141/nobackup/tremula_vs_tremuloides_paper/bwa_mem_alignment/ANGSD/summary/Fst/Chr19/Chr19.tremuloides_fixed.fst.txt 10000 10000

library(utils)
library(lattice)
.libPaths("/home/jingwang/R/x86_64-redhat-linux-gnu-library/3.2")
library(RColorBrewer)
library(ggplot2)

setwd("/proj/b2011141/nobackup/eQTL_paper/gene_list/thetas")
eGene_core=read.table("eGene_core.thetas.txt",header=T)
eGene_noncore=read.table("eGene_noncore.thetas.txt",header=T)
non_eGene_core=read.table("non_eGene_core.thetas.txt",header=T)
non_eGene_noncore=read.table("non_eGene_noncore.thetas.txt",header=T)

###merge tables
eGene_core=eGene_core[,c("tP.norm","tajD")]
eGene_core$eGene="eGene"
eGene_core$core="core"

eGene_noncore=eGene_noncore[,c("tP.norm","tajD")]
eGene_noncore$eGene="eGene"
eGene_noncore$core="Non-core"

non_eGene_core=non_eGene_core[,c("tP.norm","tajD")]
non_eGene_core$eGene="Non-eGene"
non_eGene_core$core="core"

non_eGene_noncore=non_eGene_noncore[,c("tP.norm","tajD")]
non_eGene_noncore$eGene="Non-eGene"
non_eGene_noncore$core="Non-core"

total=rbind(eGene_core,eGene_noncore,non_eGene_core,non_eGene_noncore)

png(filename="eGene_core.pi.png",width=5,height=5,units='in',res=300)
ggplot(total)+geom_boxplot(aes(x=eGene,y=tP.norm,fill=core))+labs(y=expression(theta[pi]))+scale_fill_manual(values=c("#999999", "#E69F00"))
dev.off()
###No Outliers
#png(filename="eGene_core.pi.png",width=5,height=5,units='in',res=300)
#ggplot(total)+geom_boxplot(aes(x=eGene,y=tP.norm,fill=core),outlier.size=NA)+labs(y=expression(theta[pi]))+scale_fill_manual(values=c("#999999", "#E69F00"))+scale_y_continuous(limits=c(0,0.04))


png(filename="core_eGene.pi.png",width=5,height=5,units='in',res=300)
ggplot(total)+geom_boxplot(aes(x=core,y=tP.norm,fill=eGene))+labs(y=expression(theta[pi]))
dev.off()

###No-outliers
###Mark the significance 
##0.01<P<0.05 *
##0.001<P<0.01 **
##P<0.001  ***

#> nrow(eGene_core)
#[1] 132
#> nrow(non_eGene_core)
#[1] 1569
#> nrow(non_eGene_noncore)
#[1] 12952
#> nrow(eGene_noncore)
#[1] 5795


png(filename="core_eGene.pi.no_Outliers.png",width=5,height=5,units='in',res=300)
a=ggplot(total)+geom_boxplot(aes(x=core,y=tP.norm,fill=eGene),outlier.size=NA)+labs(y=expression(theta[pi]))+scale_y_continuous(limits=c(0,0.035))
###notify the significance for different paris of groups
wilcox.test(eGene_core$tP.norm,non_eGene_core$tP.norm)  ###p-value = 0.001804
wilcox.test(eGene_noncore$tP.norm,non_eGene_noncore$tP.norm)   ###p-value < 2.2e-16
wilcox.test(eGene_core$tP.norm,eGene_noncore$tP.norm)  ###p-value = 0.01028
wilcox.test(non_eGene_core$tP.norm,non_eGene_noncore$tP.norm)   ###p-value = 2.546e-16

df1 <- data.frame(a = c(0.8, 0.8,1.2,1.2), b = c(0.025, 0.026, 0.026, 0.025))
df2 <- data.frame(a = c(1.8, 1.8,2.2,2.2), b = c(0.025, 0.026, 0.026, 0.025))
df3 <- data.frame(a = c(0.8,0.8,1.8,1.8), b = c(0.029, 0.03, 0.03, 0.029))
df4 <- data.frame(a = c(1.2,1.2,2.2,2.2), b = c(0.032, 0.033, 0.033, 0.032))

a+geom_line(data = df1, aes(x = a, y = b)) + annotate("text", x = 1, y = 0.027, label = "**", size = 6)+geom_line(data = df2, aes(x = a, y = b)) + annotate("text", x = 2, y = 0.027, label = "***", size = 6)+geom_line(data = df3, aes(x = a, y = b)) + annotate("text", x = 1.3, y = 0.031, label = "*", size = 6)+geom_line(data = df4, aes(x = a, y = b)) + annotate("text", x = 1.7, y = 0.034, label = "***", size = 6)

dev.off()

png(filename="eGene_core.tajD.png",width=5,height=5,units='in',res=300)
ggplot(total)+geom_boxplot(aes(x=eGene,y=tajD,fill=core))+labs(y="Tajima's D")+scale_fill_manual(values=c("#999999", "#E69F00"))
dev.off()

png(filename="core_eGene.tajD.png",width=5,height=5,units='in',res=300)
ggplot(total)+geom_boxplot(aes(x=core,y=tajD,fill=eGene))+labs(y="Tajima's D")
dev.off()


###plot without outliers
png(filename="core_eGene.tajD.no_Outliers.png",width=5,height=5,units='in',res=300)
a=ggplot(total)+geom_boxplot(aes(x=core,y=tajD,fill=eGene),outlier.size=NA)+labs(y="Tajima's D")+scale_y_continuous(limits=c(-3,2))
###notify the significance for different paris of groups
wilcox.test(eGene_core$tajD,non_eGene_core$tajD)  ###p-value = 0.001664
wilcox.test(eGene_noncore$tajD,non_eGene_noncore$tajD)   ###p-value < 2.2e-16
wilcox.test(eGene_core$tajD,eGene_noncore$tajD)  ###p-value = 0.003984
wilcox.test(non_eGene_core$tajD,non_eGene_noncore$tajD)   ###p-value < 2.2e-16

df1 <- data.frame(a = c(0.8, 0.8,1.2,1.2), b = c(0.9,1,1,0.9))
df2 <- data.frame(a = c(1.8, 1.8,2.2,2.2), b = c(0.9,1,1,0.9))
df3 <- data.frame(a = c(0.8,0.8,1.8,1.8), b = c(1.2, 1.3, 1.3, 1.2))
df4 <- data.frame(a = c(1.2,1.2,2.2,2.2), b = c(1.5, 1.6, 1.6, 1.5))

a+geom_line(data = df1, aes(x = a, y = b)) + annotate("text", x = 1, y = 1.1, label = "**", size = 6)+geom_line(data = df2, aes(x = a, y = b)) + annotate("text", x = 2, y = 1.1, label = "***", size = 6)+geom_line(data = df3, aes(x = a, y = b)) + annotate("text", x = 1.3, y = 1.4, label = "**", size = 6)+geom_line(data = df4, aes(x = a, y = b)) + annotate("text", x = 1.7, y = 1.7, label = "***", size = 6)

dev.off()




