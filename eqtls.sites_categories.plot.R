library(RColorBrewer)
library(ggplot2)
library(data.table)
library(reshape2)

setwd("~/Dropbox/eQTL_paper/data/thetas")
eqtls=fread("eQTLs.gene.sites_category.summary.txt",header=T)

tP=eqtls[,c(1,2,3,4,5,6,8,10,12,14,16,18)]
tP_new=melt(tP,id.vars=c("gene","in_network","network_module","is_module_core_gene","is_egene"))

tajD=eqtls[,c(1,2,3,4,5,7,9,11,13,15,17,19)]
tajD_new=melt(tajD,id.vars=c("gene","in_network","network_module","is_module_core_gene","is_egene"))


###diversity plot
ggplot(tP_new)+geom_boxplot(aes(x=variable,y=value,fill=is_egene),outlier.size=NA)+labs(y=expression(theta[pi]),x="Sites categories")+scale_y_continuous(limits=c(0,0.05))+
#  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=0.5))+
  scale_x_discrete(labels=c("UTR3","UTR5","Intron","Upstream","Downstream","0-fold","4-fold"))
ggsave("tP.sites_categories.egene.png", width=6.5, height=4, dpi=300)

ggplot(tP_new)+geom_boxplot(aes(x=variable,y=value,fill=is_module_core_gene),outlier.size=NA)+labs(y=expression(theta[pi]),x="Sites categories")+scale_y_continuous(limits=c(0,0.05))+
scale_x_discrete(labels=c("UTR3","UTR5","Intron","Upstream","Downstream","0-fold","4-fold"))
ggsave("tP.sites_categories.core.png", width=7.5, height=4, dpi=300)

ggplot(tP_new)+geom_boxplot(aes(x=variable,y=value,fill=in_network),outlier.size=NA)+labs(y=expression(theta[pi]),x="Sites categories")+scale_y_continuous(limits=c(0,0.05))+
  #  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=0.5))+
  scale_x_discrete(labels=c("UTR3","UTR5","Intron","Upstream","Downstream","0-fold","4-fold"))
ggsave("tP.sites_categories.network.png", width=6.5, height=4, dpi=300)


###tajD
ggplot(tajD_new)+geom_boxplot(aes(x=variable,y=value,fill=is_egene),outlier.size=NA)+labs(y="Tajima's D",x="Sites categories")+scale_y_continuous(limits=c(-3,1.5))+
  #  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=0.5))+
  scale_x_discrete(labels=c("UTR3","UTR5","Intron","Upstream","Downstream","0-fold","4-fold"))
ggsave("tajD.sites_categories.egene.png", width=6.5, height=4, dpi=300)

ggplot(tajD_new)+geom_boxplot(aes(x=variable,y=value,fill=is_module_core_gene),outlier.size=NA)+labs(y="Tajima's D",x="Sites categories")+scale_y_continuous(limits=c(-3,1.5))+
  scale_x_discrete(labels=c("UTR3","UTR5","Intron","Upstream","Downstream","0-fold","4-fold"))
ggsave("tajD.sites_categories.core.png", width=7.5, height=4, dpi=300)

ggplot(tajD_new)+geom_boxplot(aes(x=variable,y=value,fill=in_network),outlier.size=NA)+labs(y="Tajima's D",x="Sites categories")+scale_y_continuous(limits=c(-3,1.5))+
  scale_x_discrete(labels=c("UTR3","UTR5","Intron","Upstream","Downstream","0-fold","4-fold"))
ggsave("tajD.sites_categories.network.png", width=7.5, height=4, dpi=300)

