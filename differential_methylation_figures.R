#####################################################
##### # Differential methylation Figures
##### 1.) Input: M values, lme results for ALL and APL
##### 2.) Produce a heat map of probes in top 10 Islands
##### 3.) Produce a venn diagram of overlap of top islands
#####################################################

setwd("/Users/Rachel/Desktop/UBC/540/Stats540/Project/R_objects")
library(gplots)
library(IlluminaHumanMethylation450k.db)
library(VennDiagram)

### M Values normalized and filtered
  load('CPGI2Probe_MList.Rdata')
  cpgi.probes.M.dat <- with(cpgi.probes.M,cbind(ALL[-ncol(ALL)], APL[-ncol(APL)], CTRL))
  rm(cpgi.probes.M)

  Island<-as.data.frame(IlluminaHumanMethylation450kCPGINAME)

## HEAT ALL (all porbes in islands with t>20)
  lme_geneset_all<-read.table(file="lme.t_ml_all.tab")
  top_lme_all<-subset(lme_geneset_all, abs(q.val)<1e-92) # t value cutoff
  top_probes_all<-subset(Island, Island$cpgiview.ucscname %in% top_lme_all$cgi)
  topdat_all<-subset(cpgi.probes.M.dat[,c(1:33,44:52)], rownames(cpgi.probes.M.dat[,c(1:33,44:52)])%in% top_probes_all$cpgiview.Probe_ID)
  topdat_all<- as.matrix(t(topdat_all))
  col<-c(rep("darkgoldenrod1", each=29),rep("forestgreen",each=13))
  heatmap.2(as.matrix(topdat_all),col=redblue(256), RowSideColors=col, scale=c("column"), margins = c(2, 10),key=F,density.info="none", trace="none",cexCol=0.5, Rowv=T, Colv=T,labCol =F,labRow=F, dendrogram="row")
  f <- as.factor(c(rep("ALL", each=29),rep("HBC",each=13)))
  legend(x="topleft", legend=levels(f), col=c("darkgoldenrod1","forestgreen"), pch=15)

## HEAT APL (all porbes in islands with t>20)
  lme_geneset_apl<-read.table(file="lme.t_ml_plp.tab")
  top_lme_apl<-subset(lme_geneset_apl, abs(q.val)<1e-105) # t value cutoff
  top_probes_apl<-subset(Island, Island$cpgiview.ucscname %in% top_lme_apl$cgi)
  topdat_apl<-subset(cpgi.probes.M.dat[,c(30:41,44:52)], rownames(cpgi.probes.M.dat[,c(30:41,44:52)])%in% top_probes_apl$cpgiview.Probe_ID)
  topdat_apl<- as.matrix(t(topdat_apl))
  col<-c(rep("forestgreen", each=4),rep("darkorchid1", each=8),rep("forestgreen",each=9))
  heatmap.2(as.matrix(topdat_apl),col=redblue(256),RowSideColors=col, scale=c("column"), margins = c(2, 10),key=F,density.info="none", trace="none",cexCol=0.5, Rowv=T, Colv=T,labCol =F,labRow=F, dendrogram="row")
  f <- as.factor(c(rep("HBC", each=4),rep("APL", each=8),rep("HBC",each=9)))
  legend(x="topleft", legend=levels(f), col=c("darkorchid1","forestgreen"), pch=15)
  


## VENN 
  ALLvctrl<-(subset(lme_geneset_all, abs(q.val)<1e-25))$cgi
  APLvctrl<-(subset(lme_geneset_apl, abs(q.val)<1e-25))$cgi
  de.islands <- list(ALL = ALLvctrl, APL = APLvctrl)
  
  pdf('venn.pdf')
  plot.new()
  venn.plot <- venn.diagram(de.islands, filename = NULL, fill = c("darkgoldenrod1", "darkorchid1"),
                            main="Overlap of Differentially Methylated Islands")
  grid.draw(venn.plot)
  dev.off()


##########################################################################################################
