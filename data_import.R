###############################################
### This file was created to detail the importation and cleaning of the 3 data sets
###  (with associated metadata), and saving them as R-data objects


#### Change this and make it your own!.. preferably should be a folder for data only
setwd('~/UBC Stats/STAT540/Group Project/Data')

### You might have to install these packages using two lines:
# source('http://bioconductor.org/biocLite.R')
# biocLite('packageName')
# If they are installed, load them:
library(GEOquery)
library(affy)
library(simpleaffy)

### Get each GSE series from GEO site (this will download the "series matrix file" 
#can take awhile so once downloaded save as .R object for later)

######################
### Get data from GEObase, extract expressions, (do once)
######################

GSE42118 <- getGEO('GSE42118') 
show(GSE42118) ## 8 APL and 2 healthy marrow
save(GSE42118, file="GSE42118_matrix.R")

GSE39141 <- getGEO('GSE39141') 
show(GSE39141) ## 33 samples (ALL and healthy B cells)
save(GSE39141, file="GSE39141_matrix.R")

GSE42865 <- getGEO('GSE42865') 
show(GSE42865) ## 16 healthy cells B cells
save(GSE42865, file="GSE42865_matrix.R")

### Extract expression matrices (turn into data frames at once) 
APL.dat <- as.data.frame(exprs(GSE42118[[1]]))
ALL.dat <- as.data.frame(exprs(GSE39141[[1]]))
CTRL.dat <- as.data.frame(exprs(GSE42865[[1]]))


### Obtain the meta-data for the samples and rename them perhaps?
APL.meta <- pData(phenoData(GSE42118[[1]]))
ALL.meta <- pData(phenoData(GSE39141[[1]]))
CTRL.meta <- pData(phenoData(GSE42865[[1]]))

#####################
### Rename samples
#####################

### OK based on the info in these tables I will create some labels:
APL.meta$description
APL.meta$Group <- c('PLR', 'PLP', rep('PLR', 3), rep('PLP', 3), rep('HBM', 2))
  ## PLR: Remission Bone Marrow; PLP: Primary Diagnosis, HM: Healthy Marrow
ALL.meta$description 
ALL.meta$Group<- c(rep('ALL', 29), rep('HBC', 4))
  ## ALL: Case; HBC: Healthy B Cells

b <- CTRL.meta$characteristics_ch1.1[nrow(CTRL.meta)] ## the label for "healthy"
### Subset both meta-data and data for control (healthy) donors
CTRL.meta <- subset(CTRL.meta, characteristics_ch1.1 == as.character(b)); rm(b)
CTRL.meta <- droplevels(CTRL.meta)

CTRL.dat <- subset(CTRL.dat, select = as.character(CTRL.meta$geo_accession)) 


### Rename column names for pete's sake!
colnames(APL.dat) <- paste(APL.meta$Group, substr(colnames(APL.dat), 4, nchar(colnames(APL.dat))), sep = '_')
colnames(ALL.dat) <- paste(ALL.meta$Group, substr(colnames(ALL.dat), 4, nchar(colnames(ALL.dat))), sep = '_')
colnames(CTRL.dat) <- paste('HBC', substr(colnames(CTRL.dat), 4, nchar(colnames(CTRL.dat))), sep = '_')

### Save the data sets
save(APL.dat, ALL.dat, CTRL.dat, file = 'All_3_sets.Rdata')
save(APL.meta, ALL.meta, CTRL.meta, file = 'All_3_metasets.Rdata')

