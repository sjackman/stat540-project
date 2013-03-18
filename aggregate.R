# Aggregate the beta values of the probes for each CpG island.
# Input: 1.) normalized beta-value data sets
#        2.) metadata
# Output: 1.) beta-value AND M-value data sets with CPGI attached
#         2.) beta-value AND M-value CPGI-level data (mean)

setwd('~/UBC Stats/STAT540/Group Project/')
library(IlluminaHumanMethylation450k.db)

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


##### First restrict all data sets to probes in mapping:
cpgi.probes.Beta <- lapply(dataList, function(x) x[cpginame$Probe_ID,])

##### Create M values:
cpgi.probes.M <- lapply(cpgi.probes.Beta, logit)

#### Tack on to data in both lists: CPGI as a factor
cpgi.probes.Beta <- lapply(cpgi.probes.Beta, function(x) 
  cbind(x, cpgi = cpginame$cpginame))
cpgi.probes.M <- lapply(cpgi.probes.M, function(x) 
  cbind(x, cpgi = cpginame$cpginame))

#### Then aggregate both sets by island --> means of both betas and Ms
#### Note: avoid NA's in mean calculation

cpgi.Beta <- lapply(cpgi.probes.Beta, function(x) 
  simplify2array(by(
    x[,-ncol(x)], list(x$cpgi), function(y) colMeans(y, na.rm=T))))

cpgi.M <- lapply(cpgi.probes.M, function(x) 
  simplify2array(
    by(x[,-ncol(x)], list(x$cpgi), function(y) colMeans(y, na.rm=T))))

#### Transpose: (want features in rows for most packages)
cpgi.Beta <- lapply(cpgi.Beta, t)
cpgi.M <- lapply(cpgi.M, t)
 

##### Save these data sets:

########### Edit the path!!!!!!!!!!!!!!!!!!!
save(cpgi.probes.Beta, file = 'Data/CPGI2Probe_betaList.Rdata')
save(cpgi.probes.M, file = 'Data/CPGI2Probe_MList.Rdata')
save(cpgi.Beta, file = 'Data/CPGI_betaList.Rdata')
save(cpgi.M, file = 'Data/CPGI_MList.Rdata')
