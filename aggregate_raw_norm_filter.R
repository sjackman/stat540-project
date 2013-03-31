
#############################################################
### Raw data
# Aggregate the beta values of the probes for each CpG island.
# Input: 1.) raw beta-value data sets
# 2.) metadata
# Output: 1.) beta-value data set with CPGI attached
# 2.) beta-value CPGI-level data (mean and median)
#############################################################

setwd("C://STAT540/project/")
library(IlluminaHumanMethylation450k.db)
library(plyr)
library(gtools)
library(lattice)
library(reshape2)

### Load in (EDIT PATH):
load('Data/All_3_metasets.Rdata')
load('Data/All_3_sets.Rdata')

dataList <- list(ALL = ALL.dat, APL = APL.dat, CTRL = CTRL.dat)
#load('~/UBC Stats/STAT540/Group Project/Data/GSE42865_matrix.R')

####### Extracting and cleaning map of probe ID to CpG islands
cpginame <- as.data.frame(IlluminaHumanMethylation450kCPGINAME)
colnames(cpginame) <- c('Probe_ID', 'cpginame')
cpginame$cpginame <- factor(cpginame$cpginame)

#### There are 309,465 probes (out of total 485,577) in 27,176 islands

### Not all probes that map passed through the previous filter:
table(cpginame$Probe_ID %in% rownames(dataList$ALL))
#TRUE 
#309465 

### Restrict mappings to exclude filtered-out probes
cpginame <- subset(cpginame, Probe_ID %in% rownames(dataList$ALL))

##### Restrict all data sets to probes in mapping:
cpgi.probes.Beta <- lapply(dataList, function(x) x[cpginame$Probe_ID,])

####################################################
#### Aggregation CpG --> Island
####################################################

#### Tack on to data in both lists: CPGI as a factor

cpgi.probes.Beta <- lapply(cpgi.probes.Beta, function(x)
  cbind(x, cpgi = cpginame$cpginame))

#### Then aggregate both sets by island --> means of both betas and Ms
#### Note: avoid NA's in mean calculation

cpgi.Beta.mean <- lapply(cpgi.probes.Beta, function(x)
  simplify2array(by(
    x[,-ncol(x)], list(x$cpgi), function(y) colMeans(y, na.rm=T))))


# start of median
#### Then aggregate belta by island --> meidan of betas
#### Note: avoid NA's in median calculation

cpgi.Beta.median <- lapply(cpgi.probes.Beta, function(x)
  simplify2array(by(
    x[,-ncol(x)], list(x$cpgi), function(y) apply(y, 2, function(z) median(z, na.rm=T)))))

# end of median

##### Save these data sets:

########### Edit the path!!!!!!!!!!!!!!!!!!!
save(cpgi.probes.Beta, file = 'Data/CPGI2Probe_betaList_raw.Rdata')
# rows = cpgi, columns = samples
save(cpgi.Beta.mean, file = 'Data/CPGI_betaMeanList__raw.Rdata')
save(cpgi.Beta.median, file = 'Data/CPGI_betaMedianList_raw.Rdata')

#************************************************************************
#############################################################
### Normalized data
# Aggregate the beta values of the probes for each CpG island.
# Input: 1.) normalized beta-value data sets
# 2.) metadata
# Output: 1.) beta-value data set with CPGI attached
# 2.) beta-value CPGI-level data (mean and median)
#############################################################

rm(list =ls())
setwd("C://STAT540/project/")
library(IlluminaHumanMethylation450k.db)
library(plyr)
library(gtools)
library(lattice)
library(reshape2)

### Load in (EDIT PATH):
load('Data/All_3_metasets.Rdata')
load('Data/All_3_sets_normalized.Rdata')

dataList <- list(ALL = ALL.norm, APL = APL.norm, CTRL = CTRL.norm)
#load('~/UBC Stats/STAT540/Group Project/Data/GSE42865_matrix.R')

####### Extracting and cleaning map of probe ID to CpG islands
cpginame <- as.data.frame(IlluminaHumanMethylation450kCPGINAME)
colnames(cpginame) <- c('Probe_ID', 'cpginame')
cpginame$cpginame <- factor(cpginame$cpginame)

#### There are 309,465 probes (out of total 485,577) in 27,176 islands

### Not all probes that map passed through the previous filter:
table(cpginame$Probe_ID %in% rownames(dataList$ALL))
#TRUE 
#309465 

### Restrict mappings to exclude filtered-out probes
cpginame <- subset(cpginame, Probe_ID %in% rownames(dataList$ALL))

##### Restrict all data sets to probes in mapping:
cpgi.probes.Beta <- lapply(dataList, function(x) x[cpginame$Probe_ID,])

####################################################
#### Aggregation CpG --> Island
####################################################

#### Tack on to data in both lists: CPGI as a factor

cpgi.probes.Beta <- lapply(cpgi.probes.Beta, function(x)
  cbind(x, cpgi = cpginame$cpginame))

#### Then aggregate both sets by island --> means of both betas and Ms
#### Note: avoid NA's in mean calculation

cpgi.Beta.mean <- lapply(cpgi.probes.Beta, function(x)
  simplify2array(by(
    x[,-ncol(x)], list(x$cpgi), function(y) colMeans(y, na.rm=T))))


# start of median
#### Then aggregate belta by island --> meidan of betas
#### Note: avoid NA's in median calculation

cpgi.Beta.median <- lapply(cpgi.probes.Beta, function(x)
  simplify2array(by(
    x[,-ncol(x)], list(x$cpgi), function(y) apply(y, 2, function(z) median(z, na.rm=T)))))

# end of median

##### Save these data sets:

########### Edit the path!!!!!!!!!!!!!!!!!!!
save(cpgi.probes.Beta, file = 'Data/CPGI2Probe_betaList_norm.Rdata')
# rows = cpgi, columns = samples
save(cpgi.Beta.mean, file = 'Data/CPGI_betaMeanList__norm.Rdata')
save(cpgi.Beta.median, file = 'Data/CPGI_betaMedianList_norm.Rdata')

