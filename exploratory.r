

#############################################
###### Exploratory Analysis

###### 1.) "Dealing with" NAs
###### 2.) Determine threshold and counts
###### 3.) Count of methylated and unmethylated probes 
###### 4.) Compare control samples across assays (scatterplots of beta-values)
###### 5.) Identify probes with difference > 0.2 between controls



setwd('~/UBC Stats/STAT540/Group Project/')
library(plyr)
library(ggplot2)

##### Load data (SET PATH!!!)
load('Data/All_3_metasets.Rdata')
load('Data/All_3_sets.Rdata')
dataList <- list(ALL = ALL.dat, APL = APL.dat, CTRL = CTRL.dat)
##########################################
### Counting and plotting NAs per sample:
##########################################

na_per_samp <- lapply(dataList, function(x) apply(x, 2, function(y) sum(is.na(y))))
na_per_samp <- lapply(na_per_samp, sort)

### Plot NAs per sample:
png('Figures/NA_per_sample.png', height = 400, width = 600)
par(mai = c(2,1,1,0.5))
plot(unlist(na_per_samp), xaxt = 'n', ylab = '# NA readings (/485577)', 
     xlab = 'Sample', col = rep(c('red', 'green', 'blue'), 
                                laply(na_per_samp, length)), 
     main = 'Number of NA readings per sample')
axis(1, at = 1:length(unlist(na_per_samp)), 
     labels = do.call(c, lapply(na_per_samp, names)), las = 2, cex = 0.8)
legend('topleft', legend = names(na_per_samp), col = c('red', 'green', 'blue'))
dev.off()

## Count NAs per probe:
na_per_probe <- lapply(dataList, function(x) apply(x, 1, function(x) sum(is.na(x))))
ldply(na_per_probe, table)

# $ALL
# 
#      0      1      2      3      4      5      6      7      8      9     10 
# 462250  20220   2358    478    142     55     30     14      7      4      5 
#     11     12     13     14     16     18     20     30 
#      5      3      1      1      1      1      1      1 
# 
# $APL
# 
#      0      1      2      3      4      5      7 
# 482070   3440     58      2      3      3      1 
# 
# $CTRL
# 
#      0      1      2      3      4 
# 483502   2008     52     12      3 

### List problematic probes (more than x NA's) along with # of NA's:
(badProbe1 <- with(na_per_probe, ALL[ALL>8]))
# cg03429968   NA.73057  NA.115870  NA.116942  NA.159765  NA.169942  NA.178116 
# 12         10         14         20         11         10         10 
# NA.187460  NA.209940  NA.247471  NA.272104  NA.280289  NA.283613  NA.291451 
# 13          9         30         16         12         18         11 
# NA.296724  NA.307777  NA.312479  NA.317955  NA.338564  NA.341252  NA.373203 
# 11         11          9         10          9          9         11 
# NA.425525  NA.427054 
# 10         12
(badProbe2 <- with(na_per_probe, APL[APL>3]))
# NA.12427   NA.51759 cg04025049 cg06973463  NA.319802  NA.356032  NA.438443 
# 5          5          4          7          5          4          4  
(badProbe3 <- with(na_per_probe, CTRL[CTRL>3]))
# NA.247471 NA.256954 NA.299865 
# 4         4         4 

### Only one intersect: NA.247471 which is in badProbe1 and badProbe3

#######################
#### Justification of threshold (rewrite of Rachel's code): aka density plots
########################

library(directlabels) # allows for clearer labelling
library(lattice)

### Calculate means across samples for each probe, and reshape for plotting 
dat.probeMeans <- unlist(lapply(dataList, function(x) 
  apply(x, 1, function(y) mean(y, na.rm = T))))
plotDat <- data.frame(Beta = dat.probeMeans, Data = rep(c('ALL', 'APL', 'CTRL'), each = nrow(ALL.dat)))

png('Density_3_series_beta.png', width = 600, height = 400)
direct.label(densityplot(~Beta, data = plotDat, groups = Data, plot.points = F, 
                         lwd = 2, panel = function(...){
                           panel.densityplot(...)
                           panel.abline(v = 0.2, col = 'red')
                           panel.abline(v = 0.7, col = 'red')
                         }, xlab = 'Beta',
                         main = 'Density of Beta-Values for 3 arrays'))
dev.off()




######################
### Count of methylated probes (rewrite of Rachel's Code)
######################

### Percentages: (excluding na's)
meth_probes <- llply(dataList, function(x) 
  apply(x, 2, function(y) sum(y>0.8, na.rm=T)/length(y)))

## Plot these:
meth_probes <- llply(meth_probes, sort)
meth <- unlist(meth_probes)
names(meth) <- my.strsplit(names(meth), '\\.', i = 2)
colors <- rep(c('red', 'green', 'blue'), laply(meth_probes, length))

png('Figures/meth_per_sample.png', height = 400, width = 600)
par(mai = c(2,1,1,0.5))
barplot(meth, col = colors, las = 2, ylab = '% Meth', main = 'Proportion methylated (beta>0.8) probes per sample')
legend('topleft', legend = names(meth_probes), fill = c('red', 'green', 'blue'))
dev.off()

######################
### Count of UNmethylated probes
######################

### Percentages: (excluding na's)
unmeth_probes <- llply(dataList, function(x) 
  apply(x, 2, function(y) sum(y<0.2, na.rm=T)/length(y)))

## Plot these:
unmeth_probes <- llply(unmeth_probes, sort)
unmeth <- unlist(unmeth_probes)
names(unmeth) <- my.strsplit(names(unmeth), '\\.', i = 2)
colors <- rep(c('red', 'green', 'blue'), laply(unmeth_probes, length))

png('Figures/unmeth_per_sample.png', height = 400, width = 600)
par(mai = c(2,1,1,0.5))
barplot(unmeth, col = colors, las = 2, ylab = '% Unmeth', main = 'Proportion unmethylated (beta<0.2) probes per sample')
legend('topleft', legend = names(unmeth_probes), fill = c('red', 'green', 'blue'))
dev.off()

########################################################
##### CTRL PAIRS ANALYSIS FOR PROBE FILTERING
########################################################

my.pairs <- function(..., textCex = 2, textCol = 'red', method = c('pearson', 'spearman')){
  method <- match.arg(method)
  cor.panel <- function(x, y) 
    text(mean(range(x)), mean(range(y)), 
         labels = format(cor(x,y, method=method), digits = 3), 
         cex = textCex, col = textCol)
  pairs(..., upper.panel = cor.panel)
}

###### PLotting

### all 5 types of control samples from all data sets:
pbmc <- names(CTRL.dat)[1:3]
lcl <- names(CTRL.dat)[4:6]
nbc <- names(CTRL.dat)[7:9]
abc <- names(ALL.dat)[grep('HBC', names(ALL.dat))]
hbm <- names(APL.dat)[grep('HBM', names(APL.dat))]

ctrls.samps.list <- list(PBMC = pbmc, LCL = lcl, NBC = nbc, ABC = abc, HBM = hbm)
ctrls.samps <- unlist(ctrls.samps.list)

## Bind all data sets into one data frame:
whole.dat <- cbind(ALL.dat, APL.dat, CTRL.dat)

### Calculate average beta-value over all sets of controls:
AVG.CTRL.dat <- lapply(ctrls.samps.list, function(x) 
  apply(whole.dat[,x], 1, function(x) mean(x, na.rm = T)))
AVG.CTRL.dat <- do.call(cbind, AVG.CTRL.dat)

#### There are still na's (when I try to calculate correlation)
### Temporarily remove probes with any NAs for calculation of correlations:
AVG.noNA.dat <- AVG.CTRL.dat[apply(AVG.CTRL.dat, 1, function(x) !any(is.na(x))), ]

### Plot
png('Figures/CtrlPairs_Spearman.png', width = 1000, height = 1000)
my.pairs(AVG.noNA.dat, lower.panel = function(...){
  smoothScatter(..., add = T)
  abline(-0.2, 1, lty = 3, col = 'red')
  abline(0.2, 1, lty = 3, col = 'red')
}, method = 'spearman')
dev.off()

########################################
##### Probe exclusion by delta > 0.2
########################################

### Obtain list for each pair of controls from allNA-filtered data

### In list of vectors, gives returns another a vector of union of elements with number of sets included for each element
listHits <- function(.list){
  refset <- unique(unlist(.list))
  inclu <- do.call(rbind, lapply(.list, function(x) refset %in% x))
  foo <- apply(inclu, 2, sum); names(foo) <- refset
  sort(foo, decreasing = T)
}

## Set up pairs:
ctrlPairs <- combn(colnames(AVG.noNA.dat), 2)
ctrlPairs.dat <- alply(ctrlPairs, 2, function(x) AVG.noNA.dat[,x])
names(ctrlPairs.dat) <- alply(ctrlPairs, 2, function(x) paste(x, collapse = '.'))

## Bad (delta > 0.2) probes for each pair
badProbe.list <- llply(ctrlPairs.dat, function(x) abs(x[,2]-x[,1]) > 0.2)

badProbe.df <- do.call(cbind, badProbe.list)
aaply(badProbe.df, 2, sum)
# PBMC.LCL PBMC.NBC PBMC.ABC PBMC.HBM  LCL.NBC  LCL.ABC  LCL.HBM  NBC.ABC  NBC.HBM 
# 20851    24339    60082    71897    61552    88973    95164    27558    45762 
# ABC.HBM 
# 14996 

### List of common "misses"
table(listHits(lapply(badProbe.list, which)))
#     1     2     3     4     5     6     7     8     9 
# 21193 20702 17256 24387 14641 15344  9606  7995   310 

