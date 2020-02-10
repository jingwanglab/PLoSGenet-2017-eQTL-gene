#! /usr/bin/Rscript --no-save --no-restore
library(dplyr)
library(utils)
library(lattice)
library(RColorBrewer)
library(ggplot2)
library(data.table)

setwd("/Users/jingwang/Dropbox/eQTL_paper/data")

eqtls=fread("eQTLs.gene.summary.txt",header=T)

###Step1: estimate diversity and tajima's D separately for non-eGenes, eGenes with only distant loci, and eGenes with local
non_egenes=eqtls[which(eqtls$is_egene==FALSE),c("tP","tP_zero","tP_four","tajD","tajD_zero","tajD_four","tP_0_4","dn","ds","dnds"),with=F]
egenes_distant=eqtls[which(eqtls$is_egene==TRUE & eqtls$is_local==FALSE & eqtls$is_distant==TRUE),c("tP","tP_zero","tP_four","tajD","tajD_zero","tajD_four","tP_0_4","dn","ds","dnds"),with=F]
egenes_local=eqtls[which(eqtls$is_egene==TRUE & eqtls$is_local==TRUE),c("tP","tP_zero","tP_four","tajD","tajD_zero","tajD_four","tP_0_4","dn","ds","dnds"),with=F]


#####plot
##* p<0.05
##** p<0.01
##*** p<0.001


###basic summary for the statistics

wilcox.test(non_egenes$tP,egenes_distant$tP,conf.int = TRUE)  ##p-value = 2.57e-11
wilcox.test(non_egenes$tP,egenes_local$tP,conf.int = TRUE)  ##p-value < 2.2e-16
wilcox.test(egenes_distant$tP,egenes_local$tP,conf.int = TRUE)  ##p-value = 0.006154

median(non_egenes$tP,na.rm=T)
quantile(non_egenes$tP,na.rm=T,c(0.025,0.975))

median(egenes_distant$tP,na.rm=T)
quantile(egenes_distant$tP,na.rm=T,c(0.025,0.975))

median(egenes_local$tP,na.rm=T)
quantile(egenes_local$tP,na.rm=T,c(0.025,0.975))


wilcox.test(non_egenes$tP_four,egenes_distant$tP_four)  ##p-value = 2.57e-11
wilcox.test(non_egenes$tP_four,egenes_local$tP_four)  ##p-value < 2.2e-16
wilcox.test(egenes_distant$tP_four,egenes_local$tP_four)  ##p-value = 0.006154

wilcox.test(non_egenes$tP_zero,egenes_distant$tP_zero)  ##p-value = 2.57e-11
wilcox.test(non_egenes$tP_zero,egenes_local$tP_zero)  ##p-value < 2.2e-16
wilcox.test(egenes_distant$tP_zero,egenes_local$tP_zero)  ##p-value = 0.006154


wilcox.test(non_egenes$ds,egenes_distant$ds)  ##p-value = 2.57e-11
wilcox.test(non_egenes$ds,egenes_local$ds)  ##p-value < 2.2e-16
wilcox.test(egenes_distant$ds,egenes_local$ds)  ##p-value = 0.006154

wilcox.test(non_egenes$dn,egenes_distant$dn)  ##p-value = 2.57e-11
wilcox.test(non_egenes$dn,egenes_local$dn)  ##p-value < 2.2e-16
wilcox.test(egenes_distant$dn,egenes_local$dn)  ##p-value = 0.006154



boxplot(non_egenes$tP_four,egenes_distant$tP_four,egenes_local$tP_four,notch=TRUE,outline=FALSE,col=c("brown1","dodgerblue2","lightblue"),names=c("non-\neGene","eGene\ndistant","eGene\nlocal"),ylim=c(0,0.035),ylab=expression(theta[4-fold]))
boxplot(non_egenes$tP_zero,egenes_distant$tP_zero,egenes_local$tP_zero,notch=TRUE,outline=FALSE,col=c("brown1","dodgerblue2","lightblue"),names=c("non-\neGene","eGene\ndistant","eGene\nlocal"),ylim=c(0,0.035),ylab=expression(theta[0-fold]))

boxplot(non_egenes$ds,egenes_distant$ds,egenes_local$ds,notch=TRUE,outline=FALSE,col=c("brown1","dodgerblue2","lightblue"),names=c("non-\neGene","eGene\ndistant","eGene\nlocal"),ylim=c(0,0.2),ylab=expression(d[S]))
boxplot(non_egenes$dn,egenes_distant$dn,egenes_local$dn,notch=TRUE,outline=FALSE,col=c("brown1","dodgerblue2","lightblue"),names=c("non-\neGene","eGene\ndistant","eGene\nlocal"),ylim=c(0,0.1),ylab=expression(d[N]))


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
mtext("A",side=3,line=0.05,adj=-0.4,font=2,cex=1)







#png(filename="eQTLs.egene_core.summary.png",width=9,height=6,units='in',res=300)
pdf("eQTLs.egene_core.summary.pdf",width=9.5,height=6.5)

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
mtext("A",side=3,line=0.05,adj=-0.4,font=2,cex=1)

boxplot(non_egenes$tajD,egenes_distant$tajD,egenes_local$tajD,notch=TRUE,outline=FALSE,col=c("brown1","dodgerblue2","lightblue"),names=c("non-\neGene","eGene\ndistant","eGene\nlocal"),ylim=c(-3,2),ylab="Tajima's D")

wilcox.test(non_egenes$tajD,egenes_distant$tajD,conf.int = TRUE)  ##p-value = 0.000231
wilcox.test(non_egenes$tajD,egenes_local$tajD,conf.int = TRUE)   ##p-value < 2.2e-16
wilcox.test(egenes_distant$tajD,egenes_local$tajD,conf.int = TRUE)   ##p-value < 2.2e-16

median(non_egenes$tajD,na.rm=T)
quantile(non_egenes$tajD,na.rm=T,c(0.025,0.975))

median(egenes_distant$tajD,na.rm=T)
quantile(egenes_distant$tajD,na.rm=T,c(0.025,0.975))

median(egenes_local$tajD,na.rm=T)
quantile(egenes_local$tajD,na.rm=T,c(0.025,0.975))


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
mtext("C",side=3,line=0.05,adj=-0.4,font=2,cex=1)

boxplot(non_egenes$tP_0_4,egenes_distant$tP_0_4,egenes_local$tP_0_4,notch=TRUE,outline=FALSE,col=c("brown1","dodgerblue2","lightblue"),names=c("non-\neGene","eGene\ndistant","eGene\nlocal"),ylim=c(0,2.2),ylab=expression(theta[0-fold]/theta[4-fold]))

wilcox.test(non_egenes$tP_0_4,egenes_distant$tP_0_4,conf.int = TRUE)  ##p-value = 1.717e-07
wilcox.test(non_egenes$tP_0_4,egenes_local$tP_0_4,conf.int = TRUE)   ## p-value = 3.222e-11
wilcox.test(egenes_distant$tP_0_4,egenes_local$tP_0_4,conf.int = TRUE)   ##p-value = 0.02538

median(non_egenes$tP_0_4,na.rm=T)
quantile(non_egenes$tP_0_4,na.rm=T,c(0.025,0.975))

median(egenes_distant$tP_0_4,na.rm=T)
quantile(egenes_distant$tP_0_4,na.rm=T,c(0.025,0.975))

median(egenes_local$tP_0_4,na.rm=T)
quantile(egenes_local$tP_0_4,na.rm=T,c(0.025,0.975))



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
mtext("E",side=3,line=0.05,adj=-0.4,font=2,cex=1)

boxplot(non_egenes$dnds,egenes_distant$dnds,egenes_local$dnds,notch=TRUE,outline=FALSE,col=c("brown1","dodgerblue2","lightblue"),names=c("non-\neGene","eGene\ndistant","eGene\nlocal"),ylim=c(0,1.7),ylab=expression(d[N]/d[S]))

wilcox.test(non_egenes$dnds,egenes_distant$dnds,conf.int = TRUE)  ##p-value = 0.0002779
wilcox.test(non_egenes$dnds,egenes_local$dnds,conf.int = TRUE)   ## p-value = 0.001139
wilcox.test(egenes_distant$dnds,egenes_local$dnds,conf.int = TRUE)   ##p-value = 0.05447

median(non_egenes$dnds,na.rm=T)
quantile(non_egenes$dnds,na.rm=T,c(0.025,0.975))

median(egenes_distant$dnds,na.rm=T)
quantile(egenes_distant$dnds,na.rm=T,c(0.025,0.975))

median(egenes_local$dnds,na.rm=T)
quantile(egenes_local$dnds,na.rm=T,c(0.025,0.975))




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
mtext("G",side=3,line=0.05,adj=-0.4,font=2,cex=1)

######core vs. non-core
egene_core=eqtls[which(eqtls$is_egene==TRUE & eqtls$is_core_gene==TRUE),c("tP","tajD","tP_0_4","dnds"),with=F]
egene_non_core=eqtls[which(eqtls$is_egene==TRUE & eqtls$is_core_gene==FALSE),c("tP","tajD","tP_0_4","dnds"),with=F]

nonegene_core=eqtls[which(eqtls$is_egene==FALSE & eqtls$is_core_gene==TRUE),c("tP","tajD","tP_0_4","dnds"),with=F]
nonegene_non_core=eqtls[which(eqtls$is_egene==FALSE & eqtls$is_core_gene==FALSE),c("tP","tajD","tP_0_4","dnds"),with=F]

colors <- brewer.pal(9,"Paired")[c(1,2,5,6)]
###tP
boxplot(nonegene_core$tP,nonegene_non_core$tP,egene_core$tP,egene_non_core$tP,notch=TRUE,outline=FALSE,col=colors[c(3,4,1,2)],names=c("core","non-\ncore","core","non-\ncore"),ylim=c(0,0.035),ylab=expression(theta[pi]))
wilcox.test(core$tP,non_core$tP)  ##p-value = 0.01036
wilcox.test(egene_core$tP,egene_non_core$tP,conf.int = TRUE) #p-value = 0.01036
wilcox.test(nonegene_core$tP,nonegene_non_core$tP,conf.int = TRUE) #p-value = 2.471e-16

median(egene_core$tP,na.rm=T)
quantile(egene_core$tP,na.rm=T,c(0.025,0.975))

median(egene_non_core$tP,na.rm=T)
quantile(egene_non_core$tP,na.rm=T,c(0.025,0.975))

median(nonegene_core$tP,na.rm=T)
quantile(nonegene_core$tP,na.rm=T,c(0.025,0.975))

median(nonegene_non_core$tP,na.rm=T)
quantile(nonegene_non_core$tP,na.rm=T,c(0.025,0.975))



segments(1,0.028,1,0.03)
segments(1,0.03,1.95,0.03)
segments(1.95,0.028,1.95,0.03)
text(1.5,0.031,"***")

segments(3,0.028,3,0.03)
segments(3,0.03,3.95,0.03)
segments(3.95,0.028,3.95,0.03)
text(3.5,0.031,"*")

segments(1,-0.0063,2,-0.0063,xpd=NA,lwd=1,col="black")
segments(3,-0.0063,4,-0.0063,xpd=NA,lwd=1,col="black")
mtext(c("non-eGene","eGene"),1,line=2.2,at=c(1.5,3.5),cex=0.6)

mtext("B",side=3,line=0.05,adj=-0.4,font=2,cex=1)

###tajD
boxplot(nonegene_core$tajD,nonegene_non_core$tajD,egene_core$tajD,egene_non_core$tajD,notch=TRUE,outline=FALSE,col=colors[c(3,4,1,2)],names=c("core","non-\ncore","core","non-\ncore"),ylim=c(-3,2),ylab="Tajima's D")
#boxplot(core$tajD,non_core$tajD,notch=TRUE,outline=FALSE,col=c("gold","grey"),names=c("core","non-core"),ylim=c(-3,2),ylab="Tajima's D")
wilcox.test(egene_core$tajD,egene_non_core$tajD,conf.int = TRUE)  ##p-value = 0.004016
wilcox.test(nonegene_core$tajD,nonegene_non_core$tajD,conf.int = TRUE)  ##p-value < 2.2e-16

median(egene_core$tajD,na.rm=T)
quantile(egene_core$tajD,na.rm=T,c(0.025,0.975))

median(egene_non_core$tajD,na.rm=T)
quantile(egene_non_core$tajD,na.rm=T,c(0.025,0.975))

median(nonegene_core$tajD,na.rm=T)
quantile(nonegene_core$tajD,na.rm=T,c(0.025,0.975))

median(nonegene_non_core$tajD,na.rm=T)
quantile(nonegene_non_core$tajD,na.rm=T,c(0.025,0.975))


segments(1,1,1,1.2)
segments(1,1.2,1.95,1.2)
segments(1.95,1,1.95,1.2)
text(1.5,1.35,"***")

segments(3,1,3,1.2)
segments(3,1.2,3.95,1.2)
segments(3.95,1,3.95,1.2)
text(3.5,1.35,"**")
mtext(c("non-eGene","eGene"),1,line=2.2,at=c(1.5,3.5),cex=0.6)
segments(1,-3.9,2,-3.9,xpd=NA,lwd=1,col="black")
segments(3,-3.9,4,-3.9,xpd=NA,lwd=1,col="black")

mtext("D",side=3,line=0.05,adj=-0.4,font=2,cex=1)

###tP_0_4
boxplot(nonegene_core$tP_0_4,nonegene_non_core$tP_0_4,egene_core$tP_0_4,egene_non_core$tP_0_4,notch=TRUE,outline=FALSE,col=colors[c(3,4,1,2)],names=c("core","non-\ncore","core","non-\ncore"),ylim=c(0,2),ylab=expression(theta[0-fold]/theta[4-fold]))
wilcox.test(egene_core$tP_0_4,egene_non_core$tP_0_4,conf.int = TRUE)  ## p-value = 0.002808
wilcox.test(nonegene_core$tP_0_4,nonegene_non_core$tP_0_4,conf.int = TRUE)  ## p-value = 6.544e-14

median(egene_core$tP_0_4,na.rm=T)
quantile(egene_core$tP_0_4,na.rm=T,c(0.025,0.975))

median(egene_non_core$tP_0_4,na.rm=T)
quantile(egene_non_core$tP_0_4,na.rm=T,c(0.025,0.975))

median(nonegene_core$tP_0_4,na.rm=T)
quantile(nonegene_core$tP_0_4,na.rm=T,c(0.025,0.975))

median(nonegene_non_core$tP_0_4,na.rm=T)
quantile(nonegene_non_core$tP_0_4,na.rm=T,c(0.025,0.975))


segments(1,1.6,1,1.7)
segments(1,1.7,1.95,1.7)
segments(1.95,1.6,1.95,1.7)
text(1.5,1.8,"***")

segments(3,1.6,3,1.7)
segments(3,1.7,3.95,1.7)
segments(3.95,1.6,3.95,1.7)
text(3.5,1.8,"**")
mtext(c("non-eGene","eGene"),1,line=2.2,at=c(1.5,3.5),cex=0.6)

segments(1,-0.35,2,-0.35,xpd=NA,lwd=1,col="black")
segments(3,-0.35,4,-0.35,xpd=NA,lwd=1,col="black")

mtext("F",side=3,line=0.05,adj=-0.4,font=2,cex=1)

##dNdS
boxplot(nonegene_core$dnds,nonegene_non_core$dnds,egene_core$dnds,egene_non_core$dnds,notch=TRUE,outline=FALSE,col=colors[c(3,4,1,2)],names=c("core","non-\ncore","core","non-\ncore"),ylim=c(0,1.8),ylab=expression(d[N]/d[S]))

wilcox.test(egene_core$dnds,egene_non_core$dnds,conf.int = TRUE)  ## p-value = 0.003408
wilcox.test(nonegene_core$dnds,nonegene_non_core$dnds,conf.int = TRUE)  ## p-value = 1.453e-14

median(egene_core$dnds,na.rm=T)
quantile(egene_core$dnds,na.rm=T,c(0.025,0.975))

median(egene_non_core$dnds,na.rm=T)
quantile(egene_non_core$dnds,na.rm=T,c(0.025,0.975))

median(nonegene_core$dnds,na.rm=T)
quantile(nonegene_core$dnds,na.rm=T,c(0.025,0.975))

median(nonegene_non_core$dnds,na.rm=T)
quantile(nonegene_non_core$dnds,na.rm=T,c(0.025,0.975))


segments(1,1.4,1,1.5)
segments(1,1.5,1.95,1.5)
segments(1.95,1.4,1.95,1.5)
text(1.5,1.6,"***")

segments(3,1.4,3,1.5)
segments(3,1.5,3.95,1.5)
segments(3.95,1.4,3.95,1.5)
text(3.5,1.6,"**")

mtext(c("non-eGene","eGene"),1,line=2.2,at=c(1.5,3.5),cex=0.6)

segments(1,-0.32,2,-0.32,xpd=NA,lwd=1,col="black")
segments(3,-0.32,4,-0.32,xpd=NA,lwd=1,col="black")

mtext("H",side=3,line=0.05,adj=-0.4,font=2,cex=1)

dev.off()




