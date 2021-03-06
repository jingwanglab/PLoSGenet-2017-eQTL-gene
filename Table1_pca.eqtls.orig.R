library(utils)
library(lattice)
library(RColorBrewer)
library(ggplot2)
library(data.table)
library(dplyr)
library(reshape2)

setwd("/Users/jingwang/Dropbox/eQTL_paper/data")

#eqtls=fread("eQTLs.gene.summary.txt",header=T)
eqtls=fread("S2_file_network_eqtl_gene_statistics.tsv",header=T)

eqtls_all=fread("S2_file_all_gene_statistics.tsv",header=T)
###choose the new one
eqtls_new=eqtls_all[which(eqtls_all$gene %in% eqtls$gene),]
eqtls$orig_expr_mean=eqtls_new$expr_mean
eqtls$orig_expr_var=eqtls_new$expr_var

setwd("orig/")

###PCA analysis for the eQTLs
#pca=prcomp(eqtls[,c(2,9,10,13,16),with=F],scale=T,center=T)
pca=prcomp(eqtls[,c(2,9,46,47,15),with=F],scale=T,center=T)
pca_variance=pca$rotation^2
var=c(0.3703,0.2211,0.2010,0.1454,0.0622)
###adjust the bar length to be consistent with the total variance
pca_variance_new=matrix(0,ncol=5,nrow=5)
for (i in c(1:5))
{
  pca_variance_new[,6-i]=pca_variance[,i]*var[i]
}
colnames(pca_variance_new)=c("PC5","PC4","PC3","PC2","PC1")
###making the plot
library(RColorBrewer)
colors <- brewer.pal(12,"Paired")

#png("pca.eqtls.png",width=7,height=4.5,units='in',res=400)
pdf("pca.eqtls.pdf",width=7.5,height=5)
par(mar=c(2,4,4,4))
barplot(pca_variance_new,horiz=T,col=colors[1:5],axes=F,xlim=c(0,0.4))
axis(3)
legend("bottomright",inset=c(-0.0,-0.0),fill=colors[1:5],legend=c("Connectivity","Core vs. non-core","Expression level","Expression variance","eGene vs. Non-eGene"),bty="n")
dev.off()

#loadings=eigen(cov(eqtls[,c(2,9,11,12,15)]))$vectors
#explvar=loadings^2

##scores
scores=pca$x
loadings=pca$rotation

####################################################
##principal component regression 
pcr=lm(eqtls$tP~scores[,1]+scores[,2]+scores[,3]+scores[,4]+scores[,5])

summary(lm(log(eqtls$tP)~scores[,1]))$r.squared
summary(lm(log(eqtls$tP)~scores[,2]))$r.squared
summary(lm(log(eqtls$tP)~scores[,3]))$r.squared
summary(lm(log(eqtls$tP)~scores[,4]))$r.squared
summary(lm(log(eqtls$tP)~scores[,5]))$r.squared

summary(lm(log(eqtls$kaks)~scores[,1]))$r.squared
  summary(lm(log(eqtls$kaks)~scores[,2]))$r.squared
  summary(lm(log(eqtls$kaks)~scores[,3]))$r.squared
  summary(lm(log(eqtls$kaks)~scores[,4]))$r.squared
  summary(lm(log(eqtls$kaks)~scores[,5]))$r.squared

  summary(lm(eqtls$tP_0_4~scores[,1]))$r.squared
  summary(lm(eqtls$tP_0_4~scores[,2]))$r.squared
  summary(lm(eqtls$tP_0_4~scores[,3]))$r.squared
  summary(lm(eqtls$tP_0_4~scores[,4]))$r.squared
  summary(lm(eqtls$tP_0_4~scores[,5]))$r.squared
 
  summary(lm(log(eqtls$tP_0_4)~scores[,1]))$r.squared
  summary(lm(log(eqtls$tP_0_4)~scores[,2]))$r.squared
  summary(lm(log(eqtls$tP_0_4)~scores[,3]))$r.squared
  summary(lm(log(eqtls$tP_0_4)~scores[,4]))$r.squared
  summary(lm(log(eqtls$tP_0_4)~scores[,5]))$r.squared

  summary(lm(eqtls$tajD~scores[,1]))$r.squared
  summary(lm(eqtls$tajD~scores[,2]))$r.squared  
  summary(lm(eqtls$tajD~scores[,3]))$r.squared
  summary(lm(eqtls$tajD~scores[,4]))$r.squared
  summary(lm(eqtls$tajD~scores[,5]))$r.squared
  
####################################################
######PCA-with loadings  
  
#biplot(scores[,1:2],loadings[,1:2],cex=NULL)
png("pca.eqtls.loadings.pc1_2.png",width=5,height=4,units='in',res=300)
par(mar=c(4,4,1,1))
plot(loadings[,"PC1"],loadings[,"PC2"],col=colors[1:5],pch=NA,xlim=c(-0.3,0.7),ylim=c(-0.1,0.8),xlab="PC1 (37.0%)",ylab="PC2 (22.1%)",cex=1.5)
abline(v=0,lty=2,col="grey")
abline(h=0,lty=2,col="grey")
for (i in c(1:5))
{
  segments(0,0,loadings[i,"PC1"],loadings[i,"PC2"])
}
points(loadings[,"PC1"],loadings[,"PC2"],pch=19,col=colors[1:5],cex=1.5)
legend("topright",inset=c(-0.0,-0.0),col=colors[1:5],pch=19,legend=c("Connectivity","Core vs. non-core","Expression level","Expression variance","eGene vs. Non-eGene"),bty="n")
dev.off()

png("pca.eqtls.loadings.pc1_3.png",width=4,height=4,units='in',res=300)
par(mar=c(4,4,1,1))
plot(loadings[,"PC1"],loadings[,"PC3"],col=colors[1:5],pch=NA,xlim=c(-0.3,0.7),ylim=c(-0.6,0.8),xlab="PC1 (37.0%)",ylab="PC3 (20.1%)",cex=1.5)
abline(v=0,lty=2,col="grey")
abline(h=0,lty=2,col="grey")
for (i in c(1:5))
{
  segments(0,0,loadings[i,"PC1"],loadings[i,"PC3"])
}
points(loadings[,"PC1"],loadings[,"PC3"],pch=19,col=colors[1:5],cex=1.5)
#legend("bottomleft",inset=c(-0.0,-0.0),col=colors[1:5],pch=19,legend=c("Connectivity","Core vs. non-core","Expression level","Expression variance","eGene vs. Non-eGene"),bty="n")
dev.off()





png("pca.eqtls.loadings.pc3_4.png",width=5,height=4,units='in',res=300)
par(mar=c(4,4,1,1))
plot(loadings[,"PC3"],loadings[,"PC4"],col=colors[1:5],pch=NA,xlim=c(-0.6,0.8),ylim=c(-0.6,0.6),xlab="PC3 (20.1%)",ylab="PC4 (14.5%)",cex=1.5)
abline(v=0,lty=2,col="grey")
abline(h=0,lty=2,col="grey")
for (i in c(1:5))
{
  segments(0,0,loadings[i,"PC3"],loadings[i,"PC4"])
}
points(loadings[,"PC3"],loadings[,"PC4"],pch=19,col=colors[1:5],cex=1.5)
#legend(-0.8,0.2,cex=0.8,inset=c(-0.0,-0.0),col=colors[1:5],pch=19,legend=c("Connectivity","Core vs. non-core","Expression level","Expression variance","eGene vs. Non-eGene"),bty="n")
dev.off()

png("pca.eqtls.loadings.pc1_5.png",width=5,height=4,units='in',res=300)
par(mar=c(4,4,1,1))
plot(loadings[,"PC1"],loadings[,"PC5"],col=colors[1:5],pch=NA,xlim=c(-0.3,0.7),ylim=c(-0.8,0.7),xlab="PC1 (37.0%)",ylab="PC5 (6.3%)",cex=1.5)
abline(v=0,lty=2,col="grey")
abline(h=0,lty=2,col="grey")
for (i in c(1:5))
{
  segments(0,0,loadings[i,"PC1"],loadings[i,"PC5"])
}
points(loadings[,"PC1"],loadings[,"PC5"],pch=19,col=colors[1:5],cex=1.5)
#legend("bottomright",cex=0.8,inset=c(-0.0,-0.0),col=colors[1:5],pch=19,legend=c("Connectivity","Core vs. non-core","Expression level","Expression variance","eGene vs. Non-eGene"),bty="n")
dev.off()



####diversity
###the first is the correlation between diversity and pc scores
cor.test(scores[,1],eqtls$tP,method="spearman")  #p-value < 2.2e-16
cor.test(scores[,2],eqtls$tP,method="spearman")  #p-value < 2.2e-16
cor.test(scores[,3],eqtls$tP,method="spearman")  #p-value < 2.2e-16
cor.test(scores[,4],eqtls$tP,method="spearman")  #p-value < 2.2e-16
cor.test(scores[,5],eqtls$tP,method="spearman")  #p-value = 4.986e-15

diversity_cor=c(-0.2380,0.1350,-0.1053,-0.0695,0.0567)
pca_variance_diversity=matrix(0,ncol=5,nrow=5)
for (i in c(1:5))
{
  pca_variance_diversity[,i]=pca_variance[,i]*diversity_cor[i]
}
colnames(pca_variance_diversity)=c("PC1","PC2","PC3","PC4","PC5")
##P<0.0001***
#png("pca.eqtls.tP.png",width=4,height=4,units='in',res=400)
pdf("pca.eqtls.tP.pdf",width=4,height=5)

par(mar=c(4,5,2,1))
barplot(pca_variance_diversity,horiz=F,col=colors[1:5],axes=T,ylim=c(-0.3,0.2),ylab=expression(paste("Rank correlation with ",theta[pi],sep="")))
text(0.7,0.19,"***")
text(1.9,0.19,"***")
text(3.1,0.19,"***")
text(4.3,0.19,"***")
text(5.5,0.19,"***")
###add correlations of pcs and the corresponding explanatory variables
#pc1
text(0.7,-0.05,"+",cex=0.8)
text(0.7,-0.15,"+",cex=0.8)
text(0.7,-0.21,"+",cex=0.8)
#text(0.7,-0.225,"+",cex=0.8)
text(0.7,-0.234,"-",cex=0.8)
#pc2
#text(1.9,0.003,"+",cex=0.8)
text(1.9,0.04,"+",cex=0.8) #variance
text(1.9,0.1,"+",cex=0.8)
#pc3
#text(3.1,-0.001,"-",cex=0.8)
text(3.1,-0.035,"+",cex=0.8)
text(3.1,-0.08,"-",cex=0.8)
text(3.1,-0.1,"+",cex=0.8)
#pc4
text(4.3,-0.007,"-",cex=0.8)
text(4.3,-0.02,"+",cex=0.8)
text(4.3,-0.04,"+",cex=0.8)
text(4.3,-0.06,"-",cex=0.8)
#pc5
text(5.5,0.015,"-",cex=0.8)
text(5.5,0.04,"+",cex=0.8)
#axis(3)
#legend("bottomright",inset=c(-0.0,-0.05),fill=colors[1:5],legend=c("Connectivity","Core vs. non-core","Expression level","Expression variance","eGene vs. Non-eGene"),bty="n")
dev.off()

###tajima's D
cor.test(scores[,1],eqtls$tajD,method="spearman")   #p-value < 2.2e-16
cor.test(scores[,2],eqtls$tajD,method="spearman")   #p-value < 2.2e-16
cor.test(scores[,3],eqtls$tajD,method="spearman")   #p-value = 0.0008127
cor.test(scores[,4],eqtls$tajD,method="spearman")   #p-value < 2.2e-16
cor.test(scores[,5],eqtls$tajD,method="spearman")   #p-value = 5.588e-13

tajD_cor=c(-0.2107,0.1261,-0.0235,-0.0816,0.0505)
pca_variance_tajD=matrix(0,ncol=5,nrow=5)
for (i in c(1:5))
{
  pca_variance_tajD[,i]=pca_variance[,i]*tajD_cor[i]
}
colnames(pca_variance_tajD)=c("PC1","PC2","PC3","PC4","PC5")
##P<0.0001***
#png("pca.eqtls.tajD.png",width=4,height=4,units='in',res=400)
pdf("pca.eqtls.tajD.pdf",width=4,height=5)
par(mar=c(4,5,2,1))
barplot(pca_variance_tajD,horiz=F,col=colors[1:5],axes=T,ylim=c(-0.3,0.2),ylab="Rank correlation with Tajima's D")
text(0.7,0.19,"***")
text(1.9,0.19,"***")
text(3.1,0.19,"**")
text(4.3,0.19,"***")
text(5.5,0.19,"***")

#pc1
text(0.7,-0.05,"+",cex=0.8)
text(0.7,-0.13,"+",cex=0.8)
text(0.7,-0.185,"+",cex=0.8)
#text(0.7,-0.225,"+",cex=0.8)
text(0.7,-0.205,"-",cex=0.8)
#pc2
#text(1.9,0.003,"+",cex=0.8)
text(1.9,0.035,"+",cex=0.8) #variance
text(1.9,0.095,"+",cex=0.8)
#pc3
#text(3.1,-0.001,"-",cex=0.8)
text(3.1,-0.008,"+",cex=0.8)
text(3.1,-0.018,"-",cex=0.8)
#text(3.1,-0.1,"+",cex=0.8)
#pc4
text(4.3,-0.01,"-",cex=0.8)
text(4.3,-0.028,"+",cex=0.8)
text(4.3,-0.045,"+",cex=0.8)
text(4.3,-0.07,"-",cex=0.8)
#pc5
text(5.5,0.015,"-",cex=0.8)
text(5.5,0.04,"+",cex=0.8)


#axis(3)
#legend("bottomright",inset=c(-0.0,-0.05),fill=colors[1:5],legend=c("Connectivity","Core vs. non-core","Expression level","Expression variance","eGene vs. Non-eGene"),bty="n")
dev.off()

####tP_0_4
cor.test(scores[,1],eqtls$tP_0_4,method="spearman")   #p-value < 2.2e-16
cor.test(scores[,2],eqtls$tP_0_4,method="spearman")   #p-value = 0.0001959
cor.test(scores[,3],eqtls$tP_0_4,method="spearman")   # p-value < 2.2e-16
cor.test(scores[,4],eqtls$tP_0_4,method="spearman")   #p-value < 2.2e-16
cor.test(scores[,5],eqtls$tP_0_4,method="spearman")   #p-value = 4.7e-06

tP_0_4_cor=c(-0.1908,-0.0268,-0.1746,-0.1030,-0.0329)  
pca_variance_tP_0_4=matrix(0,ncol=5,nrow=5)
for (i in c(1:5))
{
  pca_variance_tP_0_4[,i]=pca_variance[,i]*tP_0_4_cor[i]
}
colnames(pca_variance_tP_0_4)=c("PC1","PC2","PC3","PC4","PC5")
##P<0.0001***
#png("pca.eqtls.tP_0_4.png",width=4,height=4,units='in',res=400)
pdf("pca.eqtls.tP_0_4.pdf",width=4,height=5)
par(mar=c(4,5,2,1))
barplot(pca_variance_tP_0_4,horiz=F,col=colors[1:5],axes=T,ylim=c(-0.3,0.2),ylab=expression(paste("Rank correlation with ",theta[0-fold]/theta[4-fold],sep="")))
text(0.7,0.19,"***")
text(1.9,0.19,"**")
text(3.1,0.19,"***")
text(4.3,0.19,"***")
text(5.5,0.19,"***")

#pc1
text(0.7,-0.045,"+",cex=0.8)
text(0.7,-0.12,"+",cex=0.8)
text(0.7,-0.17,"+",cex=0.8)
#text(0.7,-0.225,"+",cex=0.8)
text(0.7,-0.185,"-",cex=0.8)
#pc2
#text(1.9,0.003,"+",cex=0.8)
text(1.9,-0.008,"+",cex=0.8) #variance
text(1.9,-0.02,"+",cex=0.8)
#pc3
#text(3.1,-0.001,"-",cex=0.8)
text(3.1,-0.06,"+",cex=0.8)
text(3.1,-0.13,"-",cex=0.8)
text(3.1,-0.165,"+",cex=0.8)
#pc4
text(4.3,-0.01,"-",cex=0.8)
text(4.3,-0.032,"+",cex=0.8)
text(4.3,-0.06,"+",cex=0.8)
text(4.3,-0.09,"-",cex=0.8)
#pc5
text(5.5,-0.01,"-",cex=0.8)
text(5.5,-0.025,"+",cex=0.8)

#axis(3)
#legend("bottomright",inset=c(-0.0,-0.05),fill=colors[1:5],legend=c("Connectivity","Core vs. non-core","Expression level","Expression variance","eGene vs. Non-eGene"),bty="n")
dev.off()


###dnds
cor.test(scores[,1],eqtls$kaks,method="spearman")  #p-value < 2.2e-16
cor.test(scores[,2],eqtls$kaks,method="spearman")  #p-value = 6.672e-08
cor.test(scores[,3],eqtls$kaks,method="spearman")  #p-value < 2.2e-16
cor.test(scores[,4],eqtls$kaks,method="spearman")  #p-value < 2.2e-16
cor.test(scores[,5],eqtls$kaks,method="spearman")  #p-value = 0.001008

dnds_cor=c(-0.2277,-0.0670,-0.2446,-0.1273,-0.0226)  
pca_variance_dnds=matrix(0,ncol=5,nrow=5)
for (i in c(1:5))
{
  pca_variance_dnds[,i]=pca_variance[,i]*dnds_cor[i]
}
colnames(pca_variance_dnds)=c("PC1","PC2","PC3","PC4","PC5")
##P<0.0001***
#png("pca.eqtls.dnds.png",width=4,height=4,units='in',res=400)
pdf("pca.eqtls.dnds.pdf",width=4,height=5)
par(mar=c(4,5,2,1))
barplot(pca_variance_dnds,horiz=F,col=colors[1:5],axes=T,ylim=c(-0.3,0.2),ylab=expression(paste("Rank correlation with ",d[N]/d[S],sep="")))
text(0.7,0.19,"***")
text(1.9,0.19,"***")
text(3.1,0.19,"***")
text(4.3,0.19,"***")
text(5.5,0.19,"*")

#pc1
text(0.7,-0.05,"+",cex=0.8)
text(0.7,-0.15,"+",cex=0.8)
text(0.7,-0.2,"+",cex=0.8)
#text(0.7,-0.225,"+",cex=0.8)
text(0.7,-0.22,"-",cex=0.8)
#pc2
#text(1.9,0.003,"+",cex=0.8)
text(1.9,-0.02,"+",cex=0.8) #variance
text(1.9,-0.05,"+",cex=0.8)
#pc3
#text(3.1,-0.001,"-",cex=0.8)
text(3.1,-0.08,"+",cex=0.8)
text(3.1,-0.18,"-",cex=0.8)
text(3.1,-0.23,"+",cex=0.8)
#pc4
text(4.3,-0.01,"-",cex=0.8)
text(4.3,-0.04,"+",cex=0.8)
text(4.3,-0.07,"+",cex=0.8)
text(4.3,-0.11,"-",cex=0.8)
#pc5
text(5.5,-0.007,"-",cex=0.8)
text(5.5,-0.015,"+",cex=0.8)
#axis(3)
#legend("bottomright",inset=c(-0.0,-0.05),fill=colors[1:5],legend=c("Connectivity","Core vs. non-core","Expression level","Expression variance","eGene vs. Non-eGene"),bty="n")
dev.off()




########################################################################
##partial correlations
library(ppcor)
cor.test(eqtls$tP,eqtls$kTotal)
eqtls_cor=eqtls[,c(2,9,46,47,15,38,39,44,45)]
eqtls_cor_new=eqtls_cor[complete.cases(eqtls_cor)]

##tP vs. five measures of expression
cor.test(eqtls_cor_new$tP,eqtls_cor_new$kTotal,method="spearman")
pcor.test(eqtls_cor_new$tP,eqtls_cor_new$kTotal,eqtls_cor_new[,c("is_core_gene","orig_expr_mean","orig_expr_var","is_egene")],method="spearman")

cor.test(eqtls_cor_new$tP,as.numeric(eqtls_cor_new$is_core_gene),method="spearman")
pcor.test(eqtls_cor_new$tP,as.numeric(eqtls_cor_new$is_core_gene),eqtls_cor_new[,c("kTotal","orig_expr_mean","orig_expr_var","is_egene")],method="spearman")

cor.test(eqtls_cor_new$tP,eqtls_cor_new$orig_expr_mean,method="spearman")
pcor.test(eqtls_cor_new$tP,eqtls_cor_new$orig_expr_mean,eqtls_cor_new[,c("kTotal","is_core_gene","orig_expr_var","is_egene")],method="spearman")

cor.test(eqtls_cor_new$tP,eqtls_cor_new$orig_expr_var,method="spearman")
pcor.test(eqtls_cor_new$tP,eqtls_cor_new$orig_expr_var,eqtls_cor_new[,c("kTotal","is_core_gene","orig_expr_mean","is_egene")],method="spearman")

cor.test(eqtls_cor_new$tP,as.numeric(eqtls_cor_new$is_egene),method="spearman")
pcor.test(eqtls_cor_new$tP,as.numeric(eqtls_cor_new$is_egene),eqtls_cor_new[,c("kTotal","orig_expr_mean","orig_expr_var","is_core_gene")],method="spearman")

##tajD vs. five measures of expression
cor.test(eqtls_cor_new$tajD,eqtls_cor_new$kTotal,method="spearman")
pcor.test(eqtls_cor_new$tajD,eqtls_cor_new$kTotal,eqtls_cor_new[,c("is_core_gene","orig_expr_mean","orig_expr_var","is_egene")],method="spearman")

cor.test(eqtls_cor_new$tajD,as.numeric(eqtls_cor_new$is_core_gene),method="spearman")
pcor.test(eqtls_cor_new$tajD,as.numeric(eqtls_cor_new$is_core_gene),eqtls_cor_new[,c("kTotal","orig_expr_mean","orig_expr_var","is_egene")],method="spearman")

cor.test(eqtls_cor_new$tajD,eqtls_cor_new$orig_expr_mean,method="spearman")
pcor.test(eqtls_cor_new$tajD,eqtls_cor_new$orig_expr_mean,eqtls_cor_new[,c("kTotal","is_core_gene","orig_expr_var","is_egene")],method="spearman")

cor.test(eqtls_cor_new$tajD,eqtls_cor_new$orig_expr_var,method="spearman")
pcor.test(eqtls_cor_new$tajD,eqtls_cor_new$orig_expr_var,eqtls_cor_new[,c("kTotal","is_core_gene","orig_expr_mean","is_egene")],method="spearman")

cor.test(eqtls_cor_new$tajD,as.numeric(eqtls_cor_new$is_egene),method="spearman")
pcor.test(eqtls_cor_new$tajD,as.numeric(eqtls_cor_new$is_egene),eqtls_cor_new[,c("kTotal","orig_expr_mean","orig_expr_var","is_core_gene")],method="spearman")

##tP_0_4 vs. five measures of expression
cor.test(eqtls_cor_new$tP_0_4,eqtls_cor_new$kTotal,method="spearman")
pcor.test(eqtls_cor_new$tP_0_4,eqtls_cor_new$kTotal,eqtls_cor_new[,c("is_core_gene","orig_expr_mean","orig_expr_var","is_egene")],method="spearman")
pcor.test(eqtls_cor_new$tP_0_4,eqtls_cor_new$kTotal,eqtls_cor_new[,c("orig_expr_mean")],method="spearman")

cor.test(eqtls_cor_new$tP_0_4,as.numeric(eqtls_cor_new$is_core_gene),method="spearman")
pcor.test(eqtls_cor_new$tP_0_4,as.numeric(eqtls_cor_new$is_core_gene),eqtls_cor_new[,c("kTotal","orig_expr_mean","orig_expr_var","is_egene")],method="spearman")

cor.test(eqtls_cor_new$tP_0_4,eqtls_cor_new$orig_expr_mean,method="spearman")
pcor.test(eqtls_cor_new$tP_0_4,eqtls_cor_new$orig_expr_mean,eqtls_cor_new[,c("kTotal","is_core_gene","orig_expr_var","is_egene")],method="spearman")

cor.test(eqtls_cor_new$tP_0_4,eqtls_cor_new$orig_expr_var,method="spearman")
pcor.test(eqtls_cor_new$tP_0_4,eqtls_cor_new$orig_expr_var,eqtls_cor_new[,c("kTotal","is_core_gene","orig_expr_mean","is_egene")],method="spearman")

cor.test(eqtls_cor_new$tP_0_4,as.numeric(eqtls_cor_new$is_egene),method="spearman")
pcor.test(eqtls_cor_new$tP_0_4,as.numeric(eqtls_cor_new$is_egene),eqtls_cor_new[,c("kTotal","orig_expr_mean","orig_expr_var","is_core_gene")],method="spearman")

##dnds vs. five measures of expression
cor.test(eqtls_cor_new$kaks,eqtls_cor_new$kTotal,method="spearman")
pcor.test(eqtls_cor_new$kaks,eqtls_cor_new$kTotal,eqtls_cor_new[,c("is_core_gene","orig_expr_mean","orig_expr_var","is_egene")],method="spearman")

cor.test(eqtls_cor_new$kaks,as.numeric(eqtls_cor_new$is_core_gene),method="spearman")
pcor.test(eqtls_cor_new$kaks,as.numeric(eqtls_cor_new$is_core_gene),eqtls_cor_new[,c("kTotal","orig_expr_mean","orig_expr_var","is_egene")],method="spearman")

cor.test(eqtls_cor_new$kaks,eqtls_cor_new$orig_expr_mean,method="spearman")
pcor.test(eqtls_cor_new$kaks,eqtls_cor_new$orig_expr_mean,eqtls_cor_new[,c("kTotal","is_core_gene","orig_expr_var","is_egene")],method="spearman")

cor.test(eqtls_cor_new$kaks,eqtls_cor_new$orig_expr_var,method="spearman")
pcor.test(eqtls_cor_new$kaks,eqtls_cor_new$orig_expr_var,eqtls_cor_new[,c("kTotal","is_core_gene","orig_expr_mean","is_egene")],method="spearman")

cor.test(eqtls_cor_new$kaks,as.numeric(eqtls_cor_new$is_egene),method="spearman")
pcor.test(eqtls_cor_new$kaks,as.numeric(eqtls_cor_new$is_egene),eqtls_cor_new[,c("kTotal","orig_expr_mean","orig_expr_var","is_core_gene")],method="spearman")





