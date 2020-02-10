#! /usr/bin/Rscript --no-save --no-restore

.libPaths("/home/jingwang/R/x86_64-redhat-linux-gnu-library/3.2")
library(utils)
library(lattice)
library(RColorBrewer)
library(ggplot2)
library(data.table)

setwd("/proj/b2011141/nobackup/eQTL_paper/gene_list")

eqtls=read.table("S2_file_network_eqtl_gene_statistics.tsv",header=T)
eqtls$tP=NA
eqtls$tajD=NA
eqtls$tP_zero=NA
eqtls$tP_four=NA
eqtls$tajD_zero=NA
eqtls$tajD_four=NA
eqtls$tP_0_4=NA
eqtls$dnds=NA
eqtls$dn=NA
eqtls$ds=NA


eGene_core=fread("./thetas/eGene_core.thetas.gene.txt",header=T)
eGene_core_zero_fold=fread("./thetas/0_4_fold/new/eGene_core.zero_fold.thetas.gene.txt",header=T)
eGene_core_four_fold=fread("./thetas/0_4_fold/new/eGene_core.four_fold.thetas.gene.txt",header=T)

eGene_noncore=fread("./thetas/eGene_noncore.thetas.gene.txt",header=T)
eGene_noncore_zero_fold=fread("./thetas/0_4_fold/new/eGene_noncore.zero_fold.thetas.gene.txt",header=T)
eGene_noncore_four_fold=fread("./thetas/0_4_fold/new/eGene_noncore.four_fold.thetas.gene.txt",header=T)

non_eGene_core=fread("./thetas/non_eGene_core.thetas.gene.txt",header=T)
non_eGene_core_zero_fold=fread("./thetas/0_4_fold/new/non_eGene_core.zero_fold.thetas.gene.txt",header=T)
non_eGene_core_four_fold=fread("./thetas/0_4_fold/new/non_eGene_core.four_fold.thetas.gene.txt",header=T)

non_eGene_noncore=fread("./thetas/non_eGene_noncore.thetas.gene.txt",header=T)
non_eGene_noncore_zero_fold=fread("./thetas/0_4_fold/new/non_eGene_noncore.zero_fold.thetas.gene.txt",header=T)
non_eGene_noncore_four_fold=fread("./thetas/0_4_fold/new/non_eGene_noncore.four_fold.thetas.gene.txt",header=T)

##eGene_core
eqtls[which(eqtls$gene %in% eGene_core$Gene),"tP"]=eGene_core$tP.norm
eqtls[which(eqtls$gene %in% eGene_core$Gene),"tajD"]=eGene_core$tajD

eqtls[which(eqtls$gene %in% eGene_core_zero_fold$Gene),"tP_zero"]=eGene_core_zero_fold$tP.norm
eqtls[which(eqtls$gene %in% eGene_core_zero_fold$Gene),"tajD_zero"]=eGene_core_zero_fold$tajD

eqtls[which(eqtls$gene %in% eGene_core_four_fold$Gene),"tP_four"]=eGene_core_four_fold$tP.norm
eqtls[which(eqtls$gene %in% eGene_core_four_fold$Gene),"tajD_four"]=eGene_core_four_fold$tajD

#eGene_noncore
eqtls[which(eqtls$gene %in% eGene_noncore$Gene),"tP"]=eGene_noncore[-5105,]$tP.norm
eqtls[which(eqtls$gene %in% eGene_noncore$Gene),"tajD"]=eGene_noncore[-5105,]$tajD

eqtls[which(eqtls$gene %in% eGene_noncore_zero_fold$Gene),"tP_zero"]=eGene_noncore_zero_fold[-which(duplicated(eGene_noncore_zero_fold$Gene)==TRUE),]$tP.norm
eqtls[which(eqtls$gene %in% eGene_noncore_zero_fold$Gene),"tajD_zero"]=eGene_noncore_zero_fold[-which(duplicated(eGene_noncore_zero_fold$Gene)==TRUE),]$tajD


eqtls[which(eqtls$gene %in% eGene_noncore_four_fold$Gene),"tP_four"]=eGene_noncore_four_fold[-which(duplicated(eGene_noncore_four_fold$Gene)==TRUE),]$tP.norm
eqtls[which(eqtls$gene %in% eGene_noncore_four_fold$Gene),"tajD_four"]=eGene_noncore_four_fold[-which(duplicated(eGene_noncore_four_fold$Gene)==TRUE),]$tajD

#non_eGene_core
eqtls[which(eqtls$gene %in% non_eGene_core$Gene),"tP"]=non_eGene_core$tP.norm
eqtls[which(eqtls$gene %in% non_eGene_core$Gene),"tajD"]=non_eGene_core$tajD

eqtls[which(eqtls$gene %in% non_eGene_core_zero_fold$Gene),"tP_zero"]=non_eGene_core_zero_fold[-which(duplicated(non_eGene_core_zero_fold$Gene)==TRUE),]$tP.norm
eqtls[which(eqtls$gene %in% non_eGene_core_zero_fold$Gene),"tajD_zero"]=non_eGene_core_zero_fold[-which(duplicated(non_eGene_core_zero_fold$Gene)==TRUE),]$tajD

eqtls[which(eqtls$gene %in% non_eGene_core_four_fold$Gene),"tP_four"]=non_eGene_core_four_fold[-which(duplicated(non_eGene_core_four_fold$Gene)==TRUE),]$tP.norm
eqtls[which(eqtls$gene %in% non_eGene_core_four_fold$Gene),"tajD_four"]=non_eGene_core_four_fold[-which(duplicated(non_eGene_core_four_fold$Gene)==TRUE),]$tajD

#non_eGene_noncore
eqtls[which(eqtls$gene %in% non_eGene_noncore$Gene),"tP"]=non_eGene_noncore[-which(duplicated(non_eGene_noncore$Gene)==TRUE),]$tP.norm
eqtls[which(eqtls$gene %in% non_eGene_noncore$Gene),"tajD"]=non_eGene_noncore[-which(duplicated(non_eGene_noncore$Gene)==TRUE),]$tajD

eqtls[which(eqtls$gene %in% non_eGene_noncore_zero_fold$Gene),"tP_zero"]=non_eGene_noncore_zero_fold[-which(duplicated(non_eGene_noncore_zero_fold$Gene)==TRUE),]$tP.norm
eqtls[which(eqtls$gene %in% non_eGene_noncore_zero_fold$Gene),"tajD_zero"]=non_eGene_noncore_zero_fold[-which(duplicated(non_eGene_noncore_zero_fold$Gene)==TRUE),]$tajD

eqtls[which(eqtls$gene %in% non_eGene_noncore_four_fold$Gene),"tP_four"]=non_eGene_noncore_four_fold[-which(duplicated(non_eGene_noncore_four_fold$Gene)==TRUE),]$tP.norm
eqtls[which(eqtls$gene %in% non_eGene_noncore_four_fold$Gene),"tajD_four"]=non_eGene_noncore_four_fold[-which(duplicated(non_eGene_noncore_four_fold$Gene)==TRUE),]$tajD

######read.kaks
#dnds=fread("/proj/b2011141/nobackup/eQTL_paper/gene_list/Potra_Potri/gkaks/core_eGenes/gene.dnds.txt",header=F)
#dnds$V4=gsub("\\.[0-9]","",dnds$V2)
dnds=fread("/proj/b2011141/nobackup/eQTL_paper/gene_list/Potra_Potri/gkaks/summary/gene.summary.ka.ks.kaks.txt",header=F)

eqtls[which(eqtls$gene %in% dnds$V2),]$dnds=dnds[which(dnds$V2 %in% eqtls$gene),]$V5
eqtls[which(eqtls$gene %in% dnds$V2),]$dn=dnds[which(dnds$V2 %in% eqtls$gene),]$V3
eqtls[which(eqtls$gene %in% dnds$V2),]$ds=dnds[which(dnds$V2 %in% eqtls$gene),]$V4
####tP 0/4
eqtls$tP_0_4=eqtls$tP_zero/eqtls$tP_four

#####write.table
write.table(eqtls,file="eQTLs.gene.summary.txt",sep="\t", quote=F, row.names=F, col.names=T)


###Step1: estimate diversity and tajima's D separately for non-eGenes, eGenes with only distant loci, and eGenes with local
non_egenes=eqtls[which(eqtls$is_egene==FALSE),c("tP","tajD","tP_0_4","kaks")]
egenes_distant=eqtls[which(eqtls$is_egene==TRUE & eqtls$is_local==FALSE & eqtls$is_distant==TRUE),c("tP","tajD","tP_0_4","kaks")]
egenes_local=eqtls[which(eqtls$is_egene==TRUE & eqtls$is_local==TRUE),c("tP","tajD","tP_0_4","kaks")]


#####plot
##* p<0.05
##** p<0.01
##*** p<0.001


wilcox.test(non_egenes$tP,egenes_distant$tP)  ##p-value = 2.57e-11
wilcox.test(non_egenes$tP,egenes_local$tP)  ##p-value < 2.2e-16
wilcox.test(egenes_distant$tP,egenes_local$tP)  ##p-value = 0.006154


png(filename="eQTLs.egene_core.summary.png",width=8,height=6,units='in',res=300)

par(mfrow=c(2,4))
par(mar=c(4,4.5,1.5,1.5))
par(cex.axis=0.85)

boxplot(non_egenes$tP,egenes_distant$tP,egenes_local$tP,notch=TRUE,outline=FALSE,col=c("brown1","dodgerblue2","lightblue"),names=c("non-\neGene","eGene\ndistant","eGene\nlocal"),ylim=c(0,0.035),ylab=expression(theta[pi]))
segments(1,0.03,1,0.031)
segments(1,0.031,1.95,0.031)
segments(1.95,0.03,1.95,0.031)
text(1.5,0.032,"***")
segments(2.05,0.03,2.05,0.031)
segments(2.05,0.031,3,0.031)
segments(3,0.03,3,0.031)
text(2.5,0.032,"**")
segments(1,0.033,1,0.034)
segments(1,0.034,3,0.034)
segments(3,0.033,3,0.034)
text(2,0.035,"***")
mtext("a",side=3,line=0.05,adj=-0.4,font=2,cex=1)

boxplot(non_egenes$tajD,egenes_distant$tajD,egenes_local$tajD,notch=TRUE,outline=FALSE,col=c("brown1","dodgerblue2","lightblue"),names=c("non-\neGene","eGene\ndistant","eGene\nlocal"),ylim=c(-3,2),ylab="Tajima's D")

wilcox.test(non_egenes$tajD,egenes_distant$tajD)  ##p-value = 0.000231
wilcox.test(non_egenes$tajD,egenes_local$tajD)   ##p-value < 2.2e-16
wilcox.test(egenes_distant$tajD,egenes_local$tajD)   ##p-value < 2.2e-16

segments(1,1,1,1.15)
segments(1,1.15,1.95,1.15)
segments(1.95,1,1.95,1.15)
text(1.5,1.3,"***")
segments(2.05,1,2.05,1.15)
segments(2.05,1.15,3,1.15)
segments(3,1,3,1.15)
text(2.5,1.3,"***")
segments(1,1.45,1,1.6)
segments(1,1.6,3,1.6)
segments(3,1.45,3,1.6)
text(2,1.75,"***")
mtext("b",side=3,line=0.05,adj=-0.4,font=2,cex=1)

boxplot(non_egenes$tP_0_4,egenes_distant$tP_0_4,egenes_local$tP_0_4,notch=TRUE,outline=FALSE,col=c("brown1","dodgerblue2","lightblue"),names=c("non-\neGene","eGene\ndistant","eGene\nlocal"),ylim=c(0,2.2),ylab=expression(theta[0-fold]/theta[4-fold]))

wilcox.test(non_egenes$tP_0_4,egenes_distant$tP_0_4)  ##p-value = 1.717e-07
wilcox.test(non_egenes$tP_0_4,egenes_local$tP_0_4)   ## p-value = 3.222e-11
wilcox.test(egenes_distant$tP_0_4,egenes_local$tP_0_4)   ##p-value = 0.02538

segments(1,1.8,1,1.85)
segments(1,1.85,1.95,1.85)
segments(1.95,1.8,1.95,1.85)
text(1.5,1.9,"***")
segments(2.05,1.8,2.05,1.85)
segments(2.05,1.85,3,1.85)
segments(3,1.8,3,1.85)
text(2.5,1.9,"*")
segments(1,1.95,1,2)
segments(1,2,3,2)
segments(3,1.95,3,2)
text(2,2.05,"***")
mtext("c",side=3,line=0.05,adj=-0.4,font=2,cex=1)

boxplot(non_egenes$kaks,egenes_distant$kaks,egenes_local$kaks,notch=TRUE,outline=FALSE,col=c("brown1","dodgerblue2","lightblue"),names=c("non-\neGene","eGene\ndistant","eGene\nlocal"),ylim=c(0,1.7),ylab=expression(d[N]/d[S]))

wilcox.test(non_egenes$kaks,egenes_distant$kaks)  ##p-value = 0.0002779
wilcox.test(non_egenes$kaks,egenes_local$kaks)   ## p-value = 0.001139
wilcox.test(egenes_distant$kaks,egenes_local$kaks)   ##p-value = 0.05447

segments(1,1.3,1,1.35)
segments(1,1.35,1.95,1.35)
segments(1.95,1.3,1.95,1.35)
text(1.5,1.4,"***")
segments(2.05,1.3,2.05,1.35)
segments(2.05,1.35,3,1.35)
segments(3,1.3,3,1.35)
text(2.5,1.4,"n.s.")
segments(1,1.45,1,1.5)
segments(1,1.5,3,1.5)
segments(3,1.45,3,1.5)
text(2,1.55,"**")
mtext("d",side=3,line=0.05,adj=-0.4,font=2,cex=1)

######core vs. non-core
core=eqtls[which(eqtls$is_egene==TRUE & eqtls$is_core_gene==TRUE),c("tP","tajD","tP_0_4","kaks")]
non_core=eqtls[which(eqtls$is_egene==TRUE & eqtls$is_core_gene==FALSE),c("tP","tajD","tP_0_4","kaks")]

###tP
boxplot(core$tP,non_core$tP,notch=TRUE,outline=FALSE,col=c("gold","grey"),names=c("core","non-core"),ylim=c(0,0.035),ylab=expression(theta[pi]))
wilcox.test(core$tP,non_core$tP)  ##p-value = 0.01036

segments(1,0.03,1,0.032)
segments(1,0.032,1.95,0.032)
segments(1.95,0.03,1.95,0.032)
text(1.5,0.033,"*")
mtext("e",side=3,line=0.05,adj=-0.4,font=2,cex=1)

###tajD
boxplot(core$tajD,non_core$tajD,notch=TRUE,outline=FALSE,col=c("gold","grey"),names=c("core","non-core"),ylim=c(-3,2),ylab="Tajima's D")
wilcox.test(core$tajD,non_core$tajD)  ##p-value = 0.004016

segments(1,1,1,1.2)
segments(1,1.2,1.95,1.2)
segments(1.95,1,1.95,1.2)
text(1.5,1.35,"**")
mtext("f",side=3,line=0.05,adj=-0.4,font=2,cex=1)

###tP_0_4
boxplot(core$tP_0_4,non_core$tP_0_4,notch=TRUE,outline=FALSE,col=c("gold","grey"),names=c("core","non-core"),ylim=c(0,2),ylab=expression(theta[0-fold]/theta[4-fold]))
wilcox.test(core$tP_0_4,non_core$tP_0_4)  ## p-value = 0.002808

segments(1,1.6,1,1.7)
segments(1,1.7,1.95,1.7)
segments(1.95,1.6,1.95,1.7)
text(1.5,1.8,"**")
mtext("g",side=3,line=0.05,adj=-0.4,font=2,cex=1)

##dNdS
boxplot(core$kaks,non_core$kaks,notch=TRUE,outline=FALSE,col=c("gold","grey"),names=c("core","non-core"),ylim=c(0,1.8),ylab=expression(d[N]/d[S]))
wilcox.test(core$kaks,non_core$kaks)  ## p-value = 0.003408

segments(1,1.4,1,1.5)
segments(1,1.5,1.95,1.5)
segments(1.95,1.4,1.95,1.5)
text(1.5,1.6,"**")
mtext("h",side=3,line=0.05,adj=-0.4,font=2,cex=1)

dev.off()











