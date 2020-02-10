#! /usr/bin/Rscript --no-save --no-restore

.libPaths("/home/jingwang/R/x86_64-redhat-linux-gnu-library/3.2")
library(utils)
library(lattice)
library(RColorBrewer)
library(ggplot2)
library(data.table)

setwd("/proj/b2011141/nobackup/eQTL_paper/gene_list")

eqtls=read.table("swasp_gene_table.txt",header=T)
#eqtls=read.table("S2_file_network_eqtl_gene_statistics.tsv",header=T)

eqtls$utr3_tP=NA
eqtls$utr3_tajD=NA
eqtls$utr5_tP=NA
eqtls$utr5_tajD=NA
eqtls$intron_tP=NA
eqtls$intron_tajD=NA
eqtls$upstream_tP=NA
eqtls$upstream_tajD=NA
eqtls$downstream_tP=NA
eqtls$downstream_tajD=NA
eqtls$zero_fold_tP=NA
eqtls$zero_fold_tajD=NA
eqtls$four_fold_tP=NA
eqtls$four_fold_tajD=NA

utr3=fread("/proj/b2011141/nobackup/eQTL_paper/angsd_related/out/thetas/summary/utr3.thetas.gene.txt",header=T)
utr3=utr3[!duplicated(utr3[,"Gene",with=F],fromLast=F),]

utr5=fread("/proj/b2011141/nobackup/eQTL_paper/angsd_related/out/thetas/summary/utr5.thetas.gene.txt",header=T)
utr5=utr5[!duplicated(utr5[,"Gene",with=F],fromLast=F),]

intron=fread("/proj/b2011141/nobackup/eQTL_paper/angsd_related/out/thetas/summary/intron.thetas.gene.txt",header=T)
intron=intron[!duplicated(intron[,"Gene",with=F],fromLast=F),]

upstream=fread("/proj/b2011141/nobackup/eQTL_paper/angsd_related/out/thetas/summary/upstream.thetas.gene.txt",header=T)
upstream=upstream[!duplicated(upstream[,"Gene",with=F],fromLast=F),]

downstream=fread("/proj/b2011141/nobackup/eQTL_paper/angsd_related/out/thetas/summary/downstream.thetas.gene.txt",header=T)
downstream=downstream[!duplicated(downstream[,"Gene",with=F],fromLast=F),]

zero_fold=fread("/proj/b2011141/nobackup/eQTL_paper/angsd_related/out/thetas/summary/zero_fold.thetas.gene.txt",header=T)
zero_fold=zero_fold[!duplicated(zero_fold[,"Gene",with=F],fromLast=F),]

four_fold=fread("/proj/b2011141/nobackup/eQTL_paper/angsd_related/out/thetas/summary/four_fold.thetas.gene.txt",header=T)
four_fold=four_fold[!duplicated(four_fold[,"Gene",with=F],fromLast=F),]


eqtls[which(eqtls$gene %in% utr3$Gene),"utr3_tP"]=utr3$tP.norm
eqtls[which(eqtls$gene %in% utr3$Gene),"utr3_tajD"]=utr3$tajD

eqtls[which(eqtls$gene %in% utr5$Gene),"utr5_tP"]=utr5$tP.norm
eqtls[which(eqtls$gene %in% utr5$Gene),"utr5_tajD"]=utr5$tajD

eqtls[which(eqtls$gene %in% intron$Gene),"intron_tP"]=intron$tP.norm
eqtls[which(eqtls$gene %in% intron$Gene),"intron_tajD"]=intron$tajD

eqtls[which(eqtls$gene %in% upstream$Gene),"upstream_tP"]=upstream[which(upstream$Gene %in% eqtls$gene),]$tP.norm
eqtls[which(eqtls$gene %in% upstream$Gene),"upstream_tajD"]=upstream[which(upstream$Gene %in% eqtls$gene),]$tajD

eqtls[which(eqtls$gene %in% downstream$Gene),"downstream_tP"]=downstream[which(downstream$Gene %in% eqtls$gene),]$tP.norm
eqtls[which(eqtls$gene %in% downstream$Gene),"downstream_tajD"]=downstream[which(downstream$Gene %in% eqtls$gene),]$tajD

eqtls[which(eqtls$gene %in% zero_fold$Gene),"zero_fold_tP"]=zero_fold$tP.norm
eqtls[which(eqtls$gene %in% zero_fold$Gene),"zero_fold_tajD"]=zero_fold$tajD

eqtls[which(eqtls$gene %in% four_fold$Gene),"four_fold_tP"]=four_fold$tP.norm
eqtls[which(eqtls$gene %in% four_fold$Gene),"four_fold_tajD"]=four_fold$tajD


#####write.table
write.table(eqtls,file="eQTLs.gene.sites_category.summary.txt",sep="\t", quote=F, row.names=F, col.names=T)






