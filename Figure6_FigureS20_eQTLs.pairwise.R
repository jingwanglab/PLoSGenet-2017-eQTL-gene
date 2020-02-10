library(utils)
library(lattice)
library(RColorBrewer)
library(ggplot2)
library(data.table)
library(dplyr)

setwd("/Users/Jing/Dropbox/eQTL_paper/data")


#eqtls=fread("eQTLs.gene.summary.txt",header=T)
eqtls=fread("S2_file_network_eqtl_gene_statistics.tsv",header=T)
eqtls_all=fread("S2_file_all_gene_statistics.tsv",header=T)
###choose the new one
eqtls_new=eqtls_all[which(eqtls_all$gene %in% eqtls$gene),]
eqtls$orig_expr_mean=eqtls_new$expr_mean
eqtls$orig_expr_var=eqtls_new$expr_var

###step1:-----set pairwise plot function---------
panel.cor <- function(x, y, digits=2, cex.cor)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- as.numeric(as.character(cor.test(x, y,method="spearman",na.rm=T)$estimate))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  test <- cor.test(x,y,method="spearman")
  Signif <- ifelse(round(test$p.value,3)<0.001,"P<0.001",paste("P=",round(test$p.value,3)))
  text(0.5, 0.5,cex=1,paste("r=",txt))
  text(.5, .75,cex=1,Signif)
}

pal <- colorRampPalette(c("light blue", "yellow", "red"))

panel.smooth<-function (x, y, col = "gray", bg = NA, pch = 18,
                        cex = 0.6, col.smooth = "red", span = 2/3, iter = 3, ...)
{
  colors=densCols(x,y,colramp=pal)
  points(x, y, pch = pch, col = colors, bg = bg, cex = cex)
  ok <- is.finite(x) & is.finite(y)
  if (any(ok))
    lines(stats::lowess(x[ok], y[ok], f = span, iter = iter),
          col = col.smooth, ...)
}
panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE, breaks=15)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col="blue", ...)
}

eqtls_select=filter(eqtls,tP_0_4<10,kaks<10)
eqtls_all_select=filter(eqtls_all,tP_0_4<10,kaks<10)


png("gene_expression.orig.pairwise.png",width=8,height=8,units='in',res=400)
#pairs(select(eqtls_select,kTotal,expr_mean,expr_var,H2,tP,tajD,tP_0_4,dnds),lower.panel=panel.smooth, upper.panel=panel.cor,diag.panel=panel.hist)
pairs(select(eqtls,kTotal,is_core_gene,orig_expr_mean,orig_expr_var,is_egene,tP,tajD,tP_0_4,kaks),lower.panel=panel.smooth, upper.panel=panel.cor,diag.panel=panel.hist)
dev.off()

png("gene_expression.pairwise.select.png",width=8,height=8,units='in',res=400)
#pairs(select(eqtls_select,kTotal,expr_mean,expr_var,H2,tP,tajD,tP_0_4,dnds),lower.panel=panel.smooth, upper.panel=panel.cor,diag.panel=panel.hist)
pairs(select(eqtls_select,kTotal,is_core_gene,orig_expr_mean,expr_var,is_egene,tP,tajD,tP_0_4,kaks),lower.panel=panel.smooth, upper.panel=panel.cor,diag.panel=panel.hist)
dev.off()



###specificly to eGenes
egenes=filter(eqtls,is_egene==TRUE)
egenes_select=filter(egenes,tP_0_4<10,dnds<10)

png("egenes.pairwise.png",width=8,height=8,units='in',res=400)
pairs(select(egenes_select,kTotal,expr_mean,expr_var,tP,tajD,tP_0_4,dnds),lower.panel=panel.smooth, upper.panel=panel.cor,diag.panel=panel.hist)
dev.off()


###Specifically to non-eGenes
nonegenes=filter(eqtls,is_egene==FALSE)

nonegenes_select=filter(nonegenes,tP_0_4<10,dnds<10)

png("nonegenes.pairwise.png",width=8,height=8,units='in',res=400)
pairs(select(nonegenes_select,kTotal,expr_mean,expr_var,tP,tajD,tP_0_4,dnds),lower.panel=panel.smooth, upper.panel=panel.cor,diag.panel=panel.hist)
dev.off()

###PCA analysis for the eQTLs
pca=prcomp(eqtls[,c(2,9,10,13,16),with=F],scale=T,center=T)
pca_variance=pca$rotation^2
var=c(0.3718,0.2226,0.2048,0.1387,0.0622)
###adjust the bar length to be consistent with the total variance
pca_variance_new=matrix(0,ncol=5,nrow=5)
for (i in c(1:5))
{
pca_variance_new[,6-i]=pca_variance[,i]*var[i]
}
colnames(pca_variance_new)=c("PC5","PC4","PC3","PC2","PC1")
###making the plot
library(RColorBrewer)
colors <- brewer.pal(8,"Set2")

png("pca.eqtls.png",width=7,height=4,units='in',res=300)
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

#biplot(scores[,1:2],loadings[,1:2],cex=NULL)
png("pca.eqtls.loadings.pc1_2.png",width=4,height=4,units='in',res=300)
par(mar=c(4,4,1,1))
plot(loadings[,"PC1"],loadings[,"PC2"],col=colors[1:5],pch=NA,xlim=c(-0.7,0.3),ylim=c(-0.1,0.8),xlab="PC1",ylab="PC2",cex=1.5)
abline(v=0,lty=2,col="grey")
abline(h=0,lty=2,col="grey")
for (i in c(1:5))
{
segments(0,0,loadings[i,"PC1"],loadings[i,"PC2"])
}
points(loadings[,"PC1"],loadings[,"PC2"],pch=19,col=colors[1:5],cex=1.5)
legend("topleft",inset=c(-0.0,-0.0),col=colors[1:5],pch=19,legend=c("Connectivity","Core vs. non-core","Expression level","Expression variance","eGene vs. Non-eGene"),bty="n")
dev.off()

png("pca.eqtls.loadings.pc3_4.png",width=4,height=4,units='in',res=300)
par(mar=c(4,4,1,1))
plot(loadings[,"PC3"],loadings[,"PC4"],col=colors[1:5],pch=NA,xlim=c(-0.8,0.45),ylim=c(-0.6,0.6),xlab="PC3 (20.5%)",ylab="PC4 (13.9%)",cex=1.5)
abline(v=0,lty=2,col="grey")
abline(h=0,lty=2,col="grey")
for (i in c(1:5))
{
  segments(0,0,loadings[i,"PC3"],loadings[i,"PC4"])
}
points(loadings[,"PC3"],loadings[,"PC4"],pch=19,col=colors[1:5],cex=1.5)
#legend(-0.8,0.2,cex=0.8,inset=c(-0.0,-0.0),col=colors[1:5],pch=19,legend=c("Connectivity","Core vs. non-core","Expression level","Expression variance","eGene vs. Non-eGene"),bty="n")
dev.off()

png("pca.eqtls.loadings.pc1_5.png",width=4,height=4,units='in',res=300)
par(mar=c(4,4,1,1))
plot(loadings[,"PC1"],loadings[,"PC5"],col=colors[1:5],pch=NA,xlim=c(-0.7,0.3),ylim=c(-0.8,0.7),xlab="PC1 (37.2%)",ylab="PC5 (6.2%)",cex=1.5)
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

diversity_cor=c(0.2785,0.1643,0.0678,-0.0621,0.0549)
pca_variance_diversity=matrix(0,ncol=5,nrow=5)
for (i in c(1:5))
{
  pca_variance_diversity[,i]=pca_variance[,i]*diversity_cor[i]
}
colnames(pca_variance_diversity)=c("PC1","PC2","PC3","PC4","PC5")
##
#***P<0.0001
#**P<0.01
#*P<0.05

png("pca.eqtls.tP.png",width=4,height=4,units='in',res=400)
par(mar=c(4,5,2,1))
barplot(pca_variance_diversity,horiz=F,col=colors[1:5],axes=T,ylim=c(-0.1,0.3),ylab=expression(paste("Rank correlation with ",theta[pi],sep="")))
text(0.7,0.29,"***")
text(1.9,0.29,"***")
text(3.1,0.29,"***")
text(4.3,0.29,"***")
text(5.5,0.29,"***")

#axis(3)
#legend("bottomright",inset=c(-0.0,-0.05),fill=colors[1:5],legend=c("Connectivity","Core vs. non-core","Expression level","Expression variance","eGene vs. Non-eGene"),bty="n")
dev.off()

###tajima's D
cor.test(scores[,1],eqtls$tajD,method="spearman")   #p-value < 2.2e-16
cor.test(scores[,2],eqtls$tajD,method="spearman")   #p-value < 2.2e-16
cor.test(scores[,3],eqtls$tajD,method="spearman")   #p-value = 0.3015
cor.test(scores[,4],eqtls$tajD,method="spearman")   #p-value < 2.2e-16
cor.test(scores[,5],eqtls$tajD,method="spearman")   #p-value = 1.468e-12

tajD_cor=c(0.2352,0.1421,-0.0072,-0.0778,0.0496)
pca_variance_tajD=matrix(0,ncol=5,nrow=5)
for (i in c(1:5))
{
  pca_variance_tajD[,i]=pca_variance[,i]*tajD_cor[i]
}
colnames(pca_variance_tajD)=c("PC1","PC2","PC3","PC4","PC5")
##P<0.0001***
png("pca.eqtls.tajD.png",width=4,height=4,units='in',res=400)
par(mar=c(4,5,2,1))
barplot(pca_variance_tajD,horiz=F,col=colors[1:5],axes=T,ylim=c(-0.1,0.3),ylab="Rank correlation with Tajima's D")
text(0.7,0.29,"***")
text(1.9,0.29,"***")
text(3.1,0.29,"ns")
text(4.3,0.29,"***")
text(5.5,0.29,"***")

#axis(3)
#legend("bottomright",inset=c(-0.0,-0.05),fill=colors[1:5],legend=c("Connectivity","Core vs. non-core","Expression level","Expression variance","eGene vs. Non-eGene"),bty="n")
dev.off()

####tP_0_4
cor.test(scores[,1],eqtls$tP_0_4,method="spearman")   #p-value < 2.2e-16
cor.test(scores[,2],eqtls$tP_0_4,method="spearman")   #p-value = 0.1425
cor.test(scores[,3],eqtls$tP_0_4,method="spearman")   # p-value < 2.2e-16
cor.test(scores[,4],eqtls$tP_0_4,method="spearman")   #p-value < 2.2e-16
cor.test(scores[,5],eqtls$tP_0_4,method="spearman")   #p-value = 2.338e-06

tP_0_4_cor=c(0.2065,-0.0106,0.1520,-0.1022,-0.0340)  
pca_variance_tP_0_4=matrix(0,ncol=5,nrow=5)
for (i in c(1:5))
{
  pca_variance_tP_0_4[,i]=pca_variance[,i]*tP_0_4_cor[i]
}
colnames(pca_variance_tP_0_4)=c("PC1","PC2","PC3","PC4","PC5")
##P<0.0001***
png("pca.eqtls.tP_0_4.png",width=4,height=4,units='in',res=400)
par(mar=c(4,5,2,1))
barplot(pca_variance_tP_0_4,horiz=F,col=colors[1:5],axes=T,ylim=c(-0.1,0.3),ylab=expression(paste("Rank correlation with ",theta[0-fold]/theta[4-fold],sep="")))
text(0.7,0.29,"***")
text(1.9,0.29,"ns")
text(3.1,0.29,"***")
text(4.3,0.29,"***")
text(5.5,0.29,"***")

#axis(3)
#legend("bottomright",inset=c(-0.0,-0.05),fill=colors[1:5],legend=c("Connectivity","Core vs. non-core","Expression level","Expression variance","eGene vs. Non-eGene"),bty="n")
dev.off()


###dnds
cor.test(scores[,1],eqtls$kaks,method="spearman")  #p-value < 2.2e-16
cor.test(scores[,2],eqtls$kaks,method="spearman")  #p-value = 6.672e-08
cor.test(scores[,3],eqtls$kaks,method="spearman")  #p-value < 2.2e-16
cor.test(scores[,4],eqtls$kaks,method="spearman")  #p-value < 2.2e-16
cor.test(scores[,5],eqtls$kaks,method="spearman")  #p-value = 0.0004713

dnds_cor=c(0.2438,-0.0370,0.2272,-0.1251,-0.0240)  
pca_variance_dnds=matrix(0,ncol=5,nrow=5)
for (i in c(1:5))
{
  pca_variance_dnds[,i]=pca_variance[,i]*dnds_cor[i]
}
colnames(pca_variance_dnds)=c("PC1","PC2","PC3","PC4","PC5")
##P<0.0001***
png("pca.eqtls.dnds.png",width=4,height=4,units='in',res=400)
par(mar=c(4,5,2,1))
barplot(pca_variance_dnds,horiz=F,col=colors[1:5],axes=T,ylim=c(-0.1,0.3),ylab=expression(paste("Rank correlation with ",d[N]/d[S],sep="")))
text(0.7,0.29,"***")
text(1.9,0.29,"***")
text(3.1,0.29,"***")
text(4.3,0.29,"***")
text(5.5,0.29,"**")
#axis(3)
#legend("bottomright",inset=c(-0.0,-0.05),fill=colors[1:5],legend=c("Connectivity","Core vs. non-core","Expression level","Expression variance","eGene vs. Non-eGene"),bty="n")
dev.off()



