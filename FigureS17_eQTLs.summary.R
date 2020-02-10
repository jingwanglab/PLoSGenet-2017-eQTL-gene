#! /usr/bin/Rscript --no-save --no-restore

.libPaths("/home/jingwang/R/x86_64-redhat-linux-gnu-library/3.2")

library(data.table)
setwd("/proj/b2011141/nobackup/eQTL_paper/eQTLs/maf")

all_snp=fread("/proj/b2011141/nobackup/eQTL_paper/vcf/freq2/maf/SwAsp_94_snps.freq_0_05.maf",header=T)
all_snp_new=all_snp[which(all_snp$MAF_class!="0.00-0.05" & all_snp$MAF_class!="0.05-0.10"),]

eqtls=fread("SwAsp_94_snps.freq_0_05.eqtls.no_blank.maf",header=T)
eqtls_new=eqtls[which(eqtls$MAF_class!="0.00-0.05" & eqtls$MAF_class!="0.05-0.10"),]

eqtls_new$abs_beta=abs(eqtls_new$beta)

local_eqtls=eqtls_new[which(eqtls_new$is_local==TRUE),]
distant_eqtls=eqtls_new[which(eqtls_new$is_local==FALSE),]

#table
##all_snp
all_snp_table=as.data.frame(table(all_snp_new$MAF_class))
all_snp_table$type="All snps"

names(all_snp_table)=c("MAF_class","Freq","Type")
all_snp_table$Freq=all_snp_table$Freq/sum(all_snp_table$Freq)
##local_eqtls
local_eqtls_table=as.data.frame(table(local_eqtls$MAF_class))
local_eqtls_table$type="Local eSNPs"
names(local_eqtls_table)=c("MAF_class","Freq","Type")
local_eqtls_table$Freq=local_eqtls_table$Freq/sum(local_eqtls_table$Freq)
##distant_eqtls
distant_eqtls_table=as.data.frame(table(distant_eqtls$MAF_class))
distant_eqtls_table$type="Distant eSNPs"
names(distant_eqtls_table)=c("MAF_class","Freq","Type")
distant_eqtls_table$Freq=distant_eqtls_table$Freq/sum(distant_eqtls_table$Freq)

##total
total_table=rbind(all_snp_table,local_eqtls_table,distant_eqtls_table)

##make the plot
library(gplots)
library(RColorBrewer)
library(ggplot2)

colors <- brewer.pal(9,"Set1")
png(filename="eSNP_maf.png",width=7,height=4,units='in',res=300)
ggplot(total_table,aes(x=MAF_class,y=Freq,fill=Type))+geom_bar(position="dodge",stat="identity")
dev.off()


###plot maf vs. effect size 
pal <- colorRampPalette(c("light blue", "yellow", "red"))

colors_maf=densCols(eqtls_new$MAF,eqtls_new$abs_beta,colramp=pal)
png(filename="eQTL_maf_effect.png",width=5,height=5,units='in',res=300)
par(mar=c(4,4,1,1))
###fit a loess line or the fit of a non-linear regression
loess_fit=loess(eqtls_new$abs_beta~eqtls_new$MAF)
plot(eqtls_new$MAF,eqtls_new$abs_beta,cex=.3,col=colors_maf,pch=19,xlab="eQTL minor allele frequency",ylab="Absolute effect size")
j=order(eqtls_new$MAF)
lines(eqtls_new$MAF[j],loess_fit$fitted[j],col="black",lwd=3)
text(0.4,4,expression(paste(rho,"=-0.285")^"***"))
dev.off()



