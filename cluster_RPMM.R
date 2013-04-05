#####################################################
# generate the RPMM fit data for the raw and data
#####################################################

setwd('/ubc/cs/research/irmtraud/people/jingyun/localFile_backup/Documents/stat540/project/')
library(RColorBrewer)
library(RPMM)

##### Load raw & filt data (modify the path accordingly)
load('/ubc/cs/research/irmtraud/people/jingyun/localFile_backup/Documents/stat540/project/Data/All_3_sets.Rdata')
dim(ALL.dat)
# [1] 485577     33

# merge 3 data sets and design
merged.dat <- cbind(APL.dat, ALL.dat, CTRL.dat)

# remove probes containing NA in any of the 52 samples
NAs <- apply(merged.dat, 1, function(x) any(is.na(x)))
sum(NAs)
# remove 27366 NA-containing probes, 458211 remains 
merged.dat.noNAs <- merged.dat[!NAs, ] 

# number = top 10%
top.threshold <- ceiling((dim(merged.dat)[1] - sum(NAs)) * 0.1)
# order the probes by variance across the 52 samples
probes.var <- apply(merged.dat.noNAs, 1, function(x) var(x, na.rm=TRUE))
probes.var <- sort(probes.var, decreasing=TRUE)
probes.top.10pct <- names(probes.var)[1:top.threshold]

# keep only data for the probes that rank in the top 45822 in terms of variation across the 52 samples
merged.dat.top10pct <- merged.dat[probes.top.10pct, ]
dim(merged.dat.top10pct)
head(merged.dat.top10pct)


# Fit RPMM Model on the top 10 pct of the NA-free data
mmfit.raw <- blcTree(t(as.matrix(merged.dat.top10pct)), verbose=1)

save(mmfit.raw, file = '/ubc/cs/research/irmtraud/people/jingyun/localFile_backup/Documents/stat540/project/Data/mmfit_raw.Rdata')

#####################################################
# generate the RPMM fit data for the normalized and filtered data
#####################################################

setwd('/ubc/cs/research/irmtraud/people/jingyun/localFile_backup/Documents/stat540/project/')
library(RColorBrewer)
library(RPMM)

##### Load norm & filt data (modify the path accordingly)
load('/ubc/cs/research/irmtraud/people/jingyun/localFile_backup/Documents/stat540/project/Data/All_3_sets_normAndFilt.Rdata')
dim(ALL.dat)
# [1] 435739     33

# merge 3 data sets and design
merged.dat <- cbind(APL.dat, ALL.dat, CTRL.dat)

# remove probes containing NA in any of the 52 samples
NAs <- apply(merged.dat, 1, function(x) any(is.na(x)))
# remove 25190 NA-containing probes, 410549 remains 
merged.dat.noNAs <- merged.dat[!NAs, ] 

# number = top 10%
top.threshold <- ceiling((dim(merged.dat)[1] - sum(NAs)) * 0.1)
# order the probes by variance across the 52 samples
probes.var <- apply(merged.dat.noNAs, 1, function(x) var(x, na.rm=TRUE))
probes.var <- sort(probes.var, decreasing=TRUE)
probes.top.10pct <- names(probes.var)[1:top.threshold]

# keep only data for the probes that rank in the top 41055 in terms of variation across the 52 samples
merged.dat.top10pct <- merged.dat[probes.top.10pct, ]
dim(merged.dat.top10pct)
head(merged.dat.top10pct)


# Fit RPMM Model on the top 10 pct of the NA-free data
mmfit.normFilt <- blcTree(t(as.matrix(merged.dat.top10pct)), verbose=1)

save(mmfit.normFilt, file = '/ubc/cs/research/irmtraud/people/jingyun/localFile_backup/Documents/stat540/project/Data/mmfit_normFilt.Rdata')
