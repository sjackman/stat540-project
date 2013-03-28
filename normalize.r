###################################################################
##### Normalization of 3 separate data files of clean, raw data
##### Output: Normalizated data sets
####################################################################

### 1.) Plot distribution of beta-values
### 2.) Normalize with bmiq
### 3.) Re-plot distributions

setwd('~/UBC Stats/STAT540/Group Project/')
library(plyr)
library(ggplot2)

library(IlluminaHumanMethylation450k.db)

##### Load data and normalization script (SET PATH!!!)
load('Data/All_3_metasets.Rdata')
load('Data/All_3_sets.Rdata')
dataList <- list(ALL = ALL.dat, APL = APL.dat, CTRL = CTRL.dat)
source('Scripts/BMIQ.R')
### For this to work, one must first install the package "RPMM"


#########################################
##### Set up for BMIQ Normalization
##########################################


### Requires knowledge of Type I and Type II probes:

probe_design <- IlluminaHumanMethylation450k_getProbes()
#### Type I comes in R and G.

#### Also requires vector of 1s and 2s (from above) in order of probes

type1 <- with(probe_design$I, union(R$Probe_ID, G$Probe_ID))
type2 <- probe_design$II$Probe_ID

design <- rep(2, nrow(ALL.dat))
design[rownames(ALL.dat) %in% type1] <- 1
rm(type1, type2)
#############################################
#### BMIQ Normalization
#############################################

### Normalizes sample-by-sample WITH NO NA's!!!
### So do it like this:
### 1.) Remove NAs
### 2.) Normalize
### 3.) Plot nice density plots
### 4.) Save object for later since normalizing takes about 1 minute per sample (52 total)
### 5.) Put NAs BACK IN to original spots
### 6.) Decide what to do with them later
 
### List of 3 (data sets)
### Each is list of normalized vector (without NAs) and set of NA indices
data.norm.noNAs <- lapply(dataList, function(dat){
  lapply(dat, function(x){
    na.ix <- which(is.na(x))
    samp.norm <- BMIQ(x[-na.ix], design[-na.ix], th1.v = c(0.2, 0.7))
    list(pnbeta = samp.norm, na.ix = na.ix)
  })
})
save(data.norm.noNAs, file = 'Data/normDat_rawList.Rdata')


#################################################
######## Plot each in panels a la CheckBMIQ
#################################################

data.t2norm.densities <- lapply(data.norm.noNAs, function(dat){
  ldply(dat, function(.list) 
    with(density(.list$pnbeta$all[design[-.list$na.ix] == 2]), 
         cbind(x = x, y = y)))
})

data.raw.densities <- lapply(dataList, function(dat){
  ldply(dat, function(samp){
    na.ix <- which(is.na(samp))
    t1dens <- density(samp[-na.ix][design[-na.ix] == 1])
    t2dens <- density(samp[-na.ix][design[-na.ix] == 2])
    cbind(t1x = t1dens$x, t1y = t1dens$y, t2x = t2dens$x, t2y = t2dens$y)
  })
})

data.t2norm.densities <- with(data.t2norm.densities, rbind(ALL, APL, CTRL))
data.raw.densities <- with(data.raw.densities, rbind(ALL, APL, CTRL))
all(data.t2norm.densities$.id == data.raw.densities$.id)
data.norm.densities <- cbind(data.t2norm.densities, data.raw.densities[-1])

plotDat <- with(data.norm.densities,
                data.frame(Sample = .id, x = c(x, t1x, t2x), y = c(y, t1y, t2y),
                           Type = rep(c('T2-Norm', 'T1', 'T2'), 
                                      each = nrow(data.norm.densities))))

### Finally... plotting time
colors <- brewer.pal(3, 'Set1')
plot.settings <- list(superpose.line = list(col = colors))

pdf('Figures/Normalization_bySample.pdf', height = 11, width = 8)
xyplot(y ~ x|Sample, data = subset(plotDat, grepl('ALL', Sample)), groups = Type, 
       xlab = 'Beta', main = 'Density of beta-values', type = 'l', 
       par.settings = plot.settings, 
       auto.key = list(space = 'top', columns = 3, lines = T, points = F))
xyplot(y ~ x|Sample, data = subset(plotDat, !grepl('ALL', Sample)), groups = Type, 
       xlab = 'Beta', main = 'Density of beta-values', type = 'l', 
       par.settings = plot.settings, 
       auto.key = list(space = 'top', columns = 3, lines = T, points = F))
dev.off()

###################################################
##### Put NAs back and save
##################################################
rm(list = ls())
load('Data/normDat_rawList.Rdata')

### Function for putting NAs back into indices
naMergeIn <- function(vec, ix){
  foo <- rep(0, length(vec) + length(ix))
  foo[ix] <- NA
  foo[-ix] <- vec
  foo
}

### List of 3 data sets with NAs inserted where they should be
dataList.norm <- vector('list', 3); names(dataList.norm) <- names(data.norm.noNAs)
for(i in 1:3){
  print(i)
  tmp <- data.norm.noNAs[[i]]
  dataList.norm[[i]] <- llply(tmp, function(.list) 
    with(.list, naMergeIn(pnbeta$all, na.ix)), .progress = 'text')
}

dataList.norm <- lapply(dataList.norm, data.frame)
dataList.norm <- lapply(dataList.norm, function(x) {rownames(x) <- rownames(ALL.dat); x})
ALL.norm <- dataList.norm$ALL; APL.norm <- dataList.norm$APL; CTRL.norm <- dataList.norm$CTRL
save(ALL.norm, APL.norm, CTRL.norm, file = 'Data/All_3_sets_normalized.Rdata')

rm(data.norm.noNAs, i, tmp)

###################################################
##### Plot distributions of samples
###################################################
library(ggplot2)


### Load non-normalized data 
load("All_3_sets.Rdata")
GSE42865<-data.frame(beta=rowMeans(CTRL.dat), Series=as.factor("GSE42865"),norm=as.factor("Pre-Normalization"))
GSE42118<-data.frame(beta=rowMeans(APL.dat), Series=as.factor("GSE42118"),norm=as.factor("Pre-Normalization"))
GSE39141<-data.frame(beta=rowMeans(ALL.dat), Series=as.factor("GSE39141"),norm=as.factor("Pre-Normalization"))

## Rename normalized objects
GSE42865.n<-CTRL.norm
GSE42118.n<-APL.norm
GSE39141.n<-ALL.norm

Norm.Not<-rbind(GSE42118,GSE39141,GSE42865,GSE42118.n,GSE39141.n,GSE42865.n)

png('Density_3_series_beta_norm.png', height = 500, width = 600)
ggplot(Norm.Not, aes(beta, color=Series)) + geom_density()+
  labs(title = "Comparison of Beta Value Distribution", x="",y="Beta Frequency") +facet_grid(norm ~ .)+
  theme_bw()
dev.off()

