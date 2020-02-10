#! /usr/bin/Rscript --no-save --no-restore

# File: plotThetas.R, the following is one example of how to use this script
#Rscript angsd_plot_10kbthetas.R /proj/b2011141/nobackup/tremula_vs_tremuloides_paper/bwa_mem_alignment/ANGSD/unfolded_SFS/tremula/tremula_19/tremula_19.thetas10kbwindow.gz.pestPG /proj/b2011141/nobackup/tremula_vs_tremuloides_paper/bwa_mem_alignment/ANGSD/unfolded_SFS/tremuloides/tremuloides_19/tremuloides_19.thetas10kbwindow.gz.pestPG /proj/b2011141/nobackup/tremula_vs_tremuloides_paper/bwa_mem_alignment/ANGSD/summary/Fst/Chr19/tremula_tremuloides.chr19.win.fst.statwindow10000_step10000.fst.txt /proj/b2011141/nobackup/tremula_vs_tremuloides_paper/bwa_mem_alignment/ANGSD/summary/Fst/Chr19/Chr19.tremula_fixed.fst.txt /proj/b2011141/nobackup/tremula_vs_tremuloides_paper/bwa_mem_alignment/ANGSD/summary/Fst/Chr19/Chr19.tremuloides_fixed.fst.txt 10000 10000

library(utils)
library(lattice)
library(RColorBrewer)
library(ggplot2)

#setwd("/proj/b2011141/nobackup/eQTL_paper/gene_list/thetas/0_4_fold")
setwd("~/Dropbox/eQTL_paper/data/0_4_fold")
eGene_core_zero=read.table("eGene_core.zero_fold.thetas.txt",header=T)
eGene_core_zero$core="core"
eGene_core_zero$eGene="eGene"
eGene_core_zero$fold="0-fold"
eGene_core_four=read.table("eGene_core.four_fold.thetas.txt",header=T)
eGene_core_four$core="core"
eGene_core_four$eGene="eGene"
eGene_core_four$fold="4-fold"

eGene_noncore_zero=read.table("eGene_noncore.zero_fold.thetas.txt",header=T)
eGene_noncore_zero$core="Non-core"
eGene_noncore_zero$eGene="eGene"
eGene_noncore_zero$fold="0-fold"
eGene_noncore_four=read.table("eGene_noncore.four_fold.thetas.txt",header=T)
eGene_noncore_four$core="Non-core"
eGene_noncore_four$eGene="eGene"
eGene_noncore_four$fold="4-fold"

non_eGene_core_zero=read.table("non_eGene_core.zero_fold.thetas.txt",header=T)
non_eGene_core_zero$core="core"
non_eGene_core_zero$eGene="Non-eGene"
non_eGene_core_zero$fold="0-fold"
non_eGene_core_four=read.table("non_eGene_core.four_fold.thetas.txt",header=T)
non_eGene_core_four$core="core"
non_eGene_core_four$eGene="Non-eGene"
non_eGene_core_four$fold="4-fold"

non_eGene_noncore_zero=read.table("non_eGene_noncore.zero_fold.thetas.txt",header=T)
non_eGene_noncore_zero$core="Non-core"
non_eGene_noncore_zero$eGene="Non-eGene"
non_eGene_noncore_zero$fold="0-fold"
non_eGene_noncore_four=read.table("non_eGene_noncore.four_fold.thetas.txt",header=T)
non_eGene_noncore_four$core="Non-core"
non_eGene_noncore_four$eGene="Non-eGene"
non_eGene_noncore_four$fold="4-fold"


#####total
total_eGene_core=rbind(eGene_core_zero,eGene_core_four,eGene_noncore_zero,eGene_noncore_four,non_eGene_core_zero,non_eGene_core_four,non_eGene_noncore_zero,non_eGene_noncore_four)

###merge tables
eGene_core=cbind(eGene_core_zero[,c("tP.norm","tajD")],eGene_core_four[,c("tP.norm","tajD")])
names(eGene_core)=c("tP.norm.zero_fold","tajD.zero_fold","tP.norm.four_fold","tajD.four_fold")
eGene_core$pi0_pi4=eGene_core$tP.norm.zero_fold/eGene_core$tP.norm.four_fold
eGene_core$eGene="eGene"
eGene_core$core="core"

eGene_noncore=cbind(eGene_noncore_zero[,c("tP.norm","tajD")],eGene_noncore_four[,c("tP.norm","tajD")])
names(eGene_noncore)=c("tP.norm.zero_fold","tajD.zero_fold","tP.norm.four_fold","tajD.four_fold")
eGene_noncore$pi0_pi4=eGene_noncore$tP.norm.zero_fold/eGene_noncore$tP.norm.four_fold
eGene_noncore$eGene="eGene"
eGene_noncore$core="Non-core"

non_eGene_core=cbind(non_eGene_core_zero[,c("tP.norm","tajD")],non_eGene_core_four[,c("tP.norm","tajD")])
names(non_eGene_core)=c("tP.norm.zero_fold","tajD.zero_fold","tP.norm.four_fold","tajD.four_fold")
non_eGene_core$pi0_pi4=non_eGene_core$tP.norm.zero_fold/non_eGene_core$tP.norm.four_fold
non_eGene_core$eGene="Non-eGene"
non_eGene_core$core="core"

non_eGene_noncore=cbind(non_eGene_noncore_zero[,c("tP.norm","tajD")],non_eGene_noncore_four[,c("tP.norm","tajD")])
names(non_eGene_noncore)=c("tP.norm.zero_fold","tajD.zero_fold","tP.norm.four_fold","tajD.four_fold")
non_eGene_noncore$pi0_pi4=non_eGene_noncore$tP.norm.zero_fold/non_eGene_noncore$tP.norm.four_fold
non_eGene_noncore$eGene="Non-eGene"
non_eGene_noncore$core="Non-core"

total=rbind(eGene_core,eGene_noncore,non_eGene_core,non_eGene_noncore)


colors=brewer.pal(12,"Paired")


png(filename="eGene_core.pi0_pi4.png",width=5,height=5,units='in',res=300)
ggplot(total)+geom_boxplot(aes(x=eGene,y=pi0_pi4,fill=core))+labs(y=expression(theta[0-fold]/theta[4-fold]))+scale_y_continuous(limits=c(0,2), breaks=seq(0,2,0.5), expand = c(0, 0))+scale_fill_manual(values=c("#999999", "#E69F00"))
dev.off()

png(filename="core_eGene.pi0_pi4.png",width=5,height=5,units='in',res=300)
ggplot(total)+geom_boxplot(aes(x=core,y=pi0_pi4,fill=eGene))+labs(y=expression(theta[0-fold]/theta[4-fold]))+scale_y_continuous(limits=c(0,2), breaks=seq(0,2,0.5), expand = c(0, 0))
dev.off()


png(filename="eGene_core_fold.tajD.png",width=5,height=5,units='in',res=300)
ggplot(total_eGene_core)+geom_boxplot(aes(x=eGene,y=tajD,fill=core))+labs(y="Tajima's D")+facet_wrap(~fold)+scale_fill_manual(values=c("#999999", "#E69F00"))
dev.off()

png(filename="core_eGene_fold.tajD.png",width=5,height=5,units='in',res=300)
ggplot(total_eGene_core)+geom_boxplot(aes(x=core,y=tajD,fill=eGene))+labs(y="Tajima's D")+facet_wrap(~fold)
dev.off()


png(filename="eGene_core_fold.tP.png",width=5,height=5,units='in',res=300)
ggplot(total_eGene_core)+geom_boxplot(aes(x=eGene,y=tP.norm,fill=core))+labs(y=expression(theta[pi]))+facet_wrap(~fold)+scale_fill_manual(values=c("#999999", "#E69F00"))
dev.off()

png(filename="core_eGene_fold.tP.png",width=5,height=5,units='in',res=300)
ggplot(total_eGene_core)+geom_boxplot(aes(x=core,y=tP.norm,fill=eGene))+labs(y=expression(theta[pi]))+facet_wrap(~fold)
dev.off()





