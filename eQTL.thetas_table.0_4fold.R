#! /usr/bin/Rscript --no-save --no-restore

# File: plotThetas.R, the following is one example of how to use this script
#Rscript angsd_plot_10kbthetas.R /proj/b2011141/nobackup/tremula_vs_tremuloides_paper/bwa_mem_alignment/ANGSD/unfolded_SFS/tremula/tremula_19/tremula_19.thetas10kbwindow.gz.pestPG /proj/b2011141/nobackup/tremula_vs_tremuloides_paper/bwa_mem_alignment/ANGSD/unfolded_SFS/tremuloides/tremuloides_19/tremuloides_19.thetas10kbwindow.gz.pestPG /proj/b2011141/nobackup/tremula_vs_tremuloides_paper/bwa_mem_alignment/ANGSD/summary/Fst/Chr19/tremula_tremuloides.chr19.win.fst.statwindow10000_step10000.fst.txt /proj/b2011141/nobackup/tremula_vs_tremuloides_paper/bwa_mem_alignment/ANGSD/summary/Fst/Chr19/Chr19.tremula_fixed.fst.txt /proj/b2011141/nobackup/tremula_vs_tremuloides_paper/bwa_mem_alignment/ANGSD/summary/Fst/Chr19/Chr19.tremuloides_fixed.fst.txt 10000 10000

library(utils)
library(lattice)
library(RColorBrewer)

args=commandArgs(TRUE)
folder <- args[1] ##folder
group <- args[2] ##group :eGene_core; eGene_noncore,non_eGene_core, non_eGene_noncore
zero_fold_table=paste(group,".zero_fold.thetas.gz.pestPG",sep="")
four_fold_table=paste(group,".four_fold.thetas.gz.pestPG",sep="")
#window <- 10000
#step <- 10000


#set working directory
setwd(folder)

#set window size and step size
#window <- gsub("^[a-z]+_[0-9]+.[a-z]+([0-9]+)[a-z]+.*","\\1",tremula_input)
#step <- gsub("^[a-z]+_[0-9]+.[a-z]+([0-9]+)[a-z]+.*","\\1",tremula_input)

#read table
zero_fold_thetas <- read.table(zero_fold_table)
four_fold_thetas <- read.table(four_fold_table)

#set names of columns
names(zero_fold_thetas) <- c("range", "chr", "pos", "tW", "tP", "tF", "tH", "tL", "tajD", "fulif", "fuliD", "fayH", "zengsE", "numSites")
names(four_fold_thetas) <- c("range", "chr", "pos", "tW", "tP", "tF", "tH", "tL", "tajD", "fulif", "fuliD", "fayH", "zengsE", "numSites")

#normarize theta by numSites
zero_fold_thetas <- cbind(zero_fold_thetas, zero_fold_thetas[,c("tW", "tP", "tF", "tH", "tL")] /zero_fold_thetas$numSites)
names(zero_fold_thetas) <- c("range", "chr", "pos", "tW", "tP", "tF", "tH", "tL", "tajD", "fulif", "fuliD", "fayH", "zengsE", "numSites","tW.norm", "tP.norm", "tF.norm", "tH.norm", "tL.norm")

four_fold_thetas <- cbind(four_fold_thetas, four_fold_thetas[,c("tW", "tP", "tF", "tH", "tL")] /four_fold_thetas$numSites)
names(four_fold_thetas) <- c("range", "chr", "pos", "tW", "tP", "tF", "tH", "tL", "tajD", "fulif", "fuliD", "fayH", "zengsE", "numSites","tW.norm", "tP.norm", "tF.norm", "tH.norm", "tL.norm")


#only selected the sites that have to been 20% left for analysis
minimum_sites=20
zero_fold_thetas[which(zero_fold_thetas$numSites < minimum_sites),c(9,10,11,12,13,16,17,18,19)]=rep(NA,9)
four_fold_thetas[which(four_fold_thetas$numSites < minimum_sites),c(9,10,11,12,13,16,17,18,19)]=rep(NA,9)

#four_fold_thetas=four_fold_thetas[four_fold_thetas$gene %in% zero_fold_thetas$gene,]
#zero_fold_thetas=zero_fold_thetas[zero_fold_thetas$gene %in% four_fold_thetas$gene,]

zero_fold_theta=data.frame(cbind(Chr=as.character(zero_fold_thetas$chr), Pos=zero_fold_thetas$pos, numSites=zero_fold_thetas$numSites,tW.norm=zero_fold_thetas$tW.norm,tP.norm=zero_fold_thetas$tP.norm,tajD=zero_fold_thetas$tajD,fulif=zero_fold_thetas$fulif,fuliD=zero_fold_thetas$fuliD,fayH=zero_fold_thetas$fayH,zengsE=zero_fold_thetas$zengsE))
four_fold_theta=data.frame(cbind(Chr=as.character(four_fold_thetas$chr), Pos=four_fold_thetas$pos, numSites=four_fold_thetas$numSites,tW.norm=four_fold_thetas$tW.norm,tP.norm=four_fold_thetas$tP.norm,tajD=four_fold_thetas$tajD,fulif=four_fold_thetas$fulif,fuliD=four_fold_thetas$fuliD,fayH=four_fold_thetas$fayH,zengsE=four_fold_thetas$zengsE))

write.table(zero_fold_theta, file=paste(group,".zero_fold.thetas.txt",sep="",collapse=""), sep="\t", quote=F, row.names=F, col.names=T)
write.table(four_fold_theta, file=paste(group,".four_fold.thetas.txt",sep="",collapse=""), sep="\t", quote=F, row.names=F, col.names=T)



