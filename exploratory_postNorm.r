#####################################################
##### Exploratory Analysis after data normalization
##### 1.) Input: normalized data entitled "*.norm"
##### 2.) Produce plots of exploratory analysis
##### 3.) Find CpG probes to filter and remove from all data sets
##### 4.) Output normalized, filtered, CpG-Level data
#####     NB: in output, data files revert back to name "*.dat"
#####################################################

setwd('~/UBC Stats/STAT540/Group Project/')
library(plyr)
library(ggplot2)

##### Load data (SET PATH!!!)
load('Data/All_3_metasets.Rdata')
load('Data/All_3_sets_normalized.Rdata')
normList <- list(ALL = ALL.norm, APL = APL.norm, CTRL = CTRL.norm)

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
pbmc <- names(CTRL.norm)[1:3]
lcl <- names(CTRL.norm)[4:6]
nbc <- names(CTRL.norm)[7:9]
abc <- names(ALL.norm)[grep('HBC', names(ALL.norm))]
hbm <- names(APL.norm)[grep('HBM', names(APL.norm))]

ctrls.samps.list <- list(PBMC = pbmc, LCL = lcl, NBC = nbc, ABC = abc, HBM = hbm)
ctrls.samps <- unlist(ctrls.samps.list)

## Bind all data sets into one data frame:
whole.dat <- cbind(ALL.norm, APL.norm, CTRL.norm)

### Calculate average beta-value over all sets of controls:
AVG.CTRL.norm <- lapply(ctrls.samps.list, function(x) 
  apply(whole.dat[,x], 1, function(x) mean(x, na.rm = T)))
AVG.CTRL.norm <- do.call(cbind, AVG.CTRL.norm)

#### There are still na's (when I try to calculate correlation)
### Temporarily remove probes with any NAs for calculation of correlations:
AVG.noNA.norm <- AVG.CTRL.norm[apply(AVG.CTRL.norm, 1, function(x) !any(is.na(x))), ]
### How many probes were excluded in this way?
nrow(AVG.CTRL.norm) - nrow(AVG.noNA.norm) ## 9
### If excluded, it was NA across all samples in at least one of the 5 groups of ctrls


### Plot
png('Figures/CtrlPairs_norm_Spearman.png', width = 1000, height = 1000)
my.pairs(AVG.noNA.norm, lower.panel = function(...){
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
ctrlPairs <- combn(colnames(AVG.noNA.norm), 2)
ctrlPairs.norm <- alply(ctrlPairs, 2, function(x) AVG.noNA.norm[,x])
names(ctrlPairs.norm) <- alply(ctrlPairs, 2, function(x) paste(x, collapse = '.'))

## Bad (delta > 0.2) probes for each pair
badProbe.list <- llply(ctrlPairs.norm, function(x) abs(x[,2]-x[,1]) > 0.2)

badProbe.df <- do.call(cbind, badProbe.list)
aaply(badProbe.df, 2, sum)
# PBMC.LCL PBMC.NBC PBMC.ABC PBMC.HBM  LCL.NBC  LCL.ABC  LCL.HBM  NBC.ABC 
# 48838    13741    12134    32282    49820    48630    66399       66 
# NBC.HBM  ABC.HBM 
# 16970    18972

### List of common "misses"
table(listHits(lapply(badProbe.list, which)))
#     1     2     3     4     5     6     7     8     9 
# 15082 10823 14455 35594  5781  4332  2673  1266   183 

### How many "misses", not counting NA across any 
sum(apply(badProbe.df, 1, any)) ## 90189
sum(apply(badProbe.df, 1, function(x) sum(x)>=4)) ## 49829

#########################################
##### Get list of probe which "hit" fewer than 4 pairs of ctrl probes
##### This list will also exclude CpGs that have NA across a group of controls
#########################################

keepers <- rownames(badProbe.df)[apply(badProbe.df, 1, function(x) sum(x)<4)]
### 
dataList <- lapply(normList, function(x) x[keepers,])
### Count summary:
## 485577 CpG probes to begin with
## - 9 too many NAs
## -49829 from delta < 0.2 for 4 or more ctrl pairs
## --> 435739 left (to be aggregated into islands)

### Saving normalized and CpG-filtered data sets:
ALL.dat <- dataList$ALL; APL.dat <- dataList$APL; CTRL.dat <- dataList$CTRL
save(ALL.dat, APL.dat, CTRL.dat, file = 'Data/All_3_sets_normAndFilt.Rdata.Rdata')
