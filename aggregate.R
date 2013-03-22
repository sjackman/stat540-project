# Aggregate the beta values of the probes for each CpG island.
# Input: 1.) normalized beta-value data sets
#        2.) metadata
# Output: 1.) beta-value AND M-value data sets with CPGI attached
#         2.) beta-value AND M-value CPGI-level data (mean)

setwd('~/UBC Stats/STAT540/Group Project/')
library(IlluminaHumanMethylation450k.db)
library(plyr)
library(gtools)
library(lattice)
library(reshape2)

### Load in (EDIT PATH):
load('Data/All_3_metasets.Rdata')
load('Data/All_3_sets_normAndFilt.Rdata')

dataList <- list(ALL = ALL.dat, APL = APL.dat, CTRL = CTRL.dat)
#load('~/UBC Stats/STAT540/Group Project/Data/GSE42865_matrix.R')

####### Extracting and cleaning map of probe ID to CpG islands
cpginame <- as.data.frame(IlluminaHumanMethylation450kCPGINAME)
colnames(cpginame) <- c('Probe_ID', 'cpginame')
cpginame$cpginame <- factor(cpginame$cpginame)

#### There are 309,465 probes (out of total 485,577) in 27,176 islands

### Not all probes that map passed through the previous filter:
table(cpginame$Probe_ID %in% rownames(dataList$ALL))
# FALSE   TRUE 
# 19933 289532

### Restrict mappings to exclude filtered-out probes
cpginame <- subset(cpginame, Probe_ID %in% rownames(dataList$ALL))

##### Restrict all data sets to probes in mapping:
cpgi.probes.Beta <- lapply(dataList, function(x) x[cpginame$Probe_ID,])


###############################################
##### M-value transformation
###############################################

#### Ranges:
(ranges <- ldply(cpgi.probes.Beta, function(x) range(x, na.rm = T)))
#    .id           V1        V2
# 1  ALL 2.123553e-14 0.9999967
# 2  APL 1.292368e-02 0.9906188
# 3 CTRL 3.266047e-14 0.9994150
names(ranges) <- c('Set', 'Min', 'Max')
rownames(ranges) <- ranges$Set; ranges$Set <- NULL
ranges <- as.matrix(ranges)

#### Check out:
logit(ranges)
logit(ranges + 1e-6) ## This one seems reasonable! (somewhat symmetric on the tails)
### Add epsilon to Beta-values before transforming to M-values:

cpgi.probes.M <- lapply(cpgi.probes.Beta, function(x) logit(x + 1e-6))

#### Explore a bit:
whole.M <- with(cpgi.probes.M, cbind(ALL, APL, CTRL))
whole.M.tall <- melt(t(whole.M), value.name = 'M', varnames = c('Sample', 'Probe_ID'))
png('Figures/Normalized_M_CpG_allSamps.png', width = 2000, height = 400)
par(mar = c(8, 4, 4, 2))
bwplot(M~Sample, data = whole.M.tall, panel = panel.violin,
       scales = list(x = list(rot = 90)), xlab = 'Sample', ylab = 'M values',
       main = 'Distribution of M-values from normalized and filtered Beta-values')
dev.off()
rm(whole.M.tall)

####################################################
####  Aggregation CpG --> Island
####################################################

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
