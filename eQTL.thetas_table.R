#! /usr/bin/Rscript --no-save --no-restore

# File: plotThetas.R, the following is one example of how to use this script
#Rscript angsd_plot_10kbthetas.R /proj/b2011141/nobackup/tremula_vs_tremuloides_paper/bwa_mem_alignment/ANGSD/unfolded_SFS/tremula/tremula_19/tremula_19.thetas10kbwindow.gz.pestPG /proj/b2011141/nobackup/tremula_vs_tremuloides_paper/bwa_mem_alignment/ANGSD/unfolded_SFS/tremuloides/tremuloides_19/tremuloides_19.thetas10kbwindow.gz.pestPG /proj/b2011141/nobackup/tremula_vs_tremuloides_paper/bwa_mem_alignment/ANGSD/summary/Fst/Chr19/tremula_tremuloides.chr19.win.fst.statwindow10000_step10000.fst.txt /proj/b2011141/nobackup/tremula_vs_tremuloides_paper/bwa_mem_alignment/ANGSD/summary/Fst/Chr19/Chr19.tremula_fixed.fst.txt /proj/b2011141/nobackup/tremula_vs_tremuloides_paper/bwa_mem_alignment/ANGSD/summary/Fst/Chr19/Chr19.tremuloides_fixed.fst.txt 10000 10000

library(utils)
library(lattice)
library(RColorBrewer)

args=commandArgs(TRUE)
folder <- args[1] ##folder
group <- args[2] ##group :eGene_core; eGene_noncore,non_eGene_core, non_eGene_noncore
table=paste(group,".thetas.gz.pestPG",sep="")
#window <- 10000
#step <- 10000


#set working directory
setwd(folder)

#set window size and step size
#window <- gsub("^[a-z]+_[0-9]+.[a-z]+([0-9]+)[a-z]+.*","\\1",tremula_input)
#step <- gsub("^[a-z]+_[0-9]+.[a-z]+([0-9]+)[a-z]+.*","\\1",tremula_input)

#read table
thetas <- read.table(table)

#set names of columns
names(thetas) <- c("range", "chr", "pos", "tW", "tP", "tF", "tH", "tL", "tajD", "fulif", "fuliD", "fayH", "zengsE", "numSites")

#normarize theta by numSites
thetas <- cbind(thetas, thetas[,c("tW", "tP", "tF", "tH", "tL")] / thetas$numSites)
names(thetas) <- c("range", "chr", "pos", "tW", "tP", "tF", "tH", "tL", "tajD", "fulif", "fuliD", "fayH", "zengsE", "numSites", "tW.norm", "tP.norm", "tF.norm", "tH.norm", "tL.norm")

#only selected the sites that have to been 20% left for analysis
minimum_sites=50
thetas[which(thetas$numSites < minimum_sites),c(9,10,11,12,13,15,16,17,18,19)]=rep(NA,10)

theta=data.frame(cbind(Chr=as.character(thetas$chr), Pos=thetas$pos, numSites=thetas$numSites,tW.norm=thetas$tW.norm,tP.norm=thetas$tP.norm,tajD=thetas$tajD,fulif=thetas$fulif,fuliD=thetas$fuliD,fayH=thetas$fayH,zengsE=thetas$zengsE))

write.table(theta, file=paste(group,".thetas.txt",sep="",collapse=""), sep="\t", quote=F, row.names=F, col.names=T)



