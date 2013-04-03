# Fit a linear mixed model

setwd('~/UBC Stats/STAT540/Group Project/')
load('Data/CPGI2Probe_MList.Rdata')
source('coord_to_gene.R')
library(doParallel)
library(foreach)
library(lme4)
library(plyr)
library(reshape2)
library(nlme)

# Enable parallel computation.
registerDoParallel()

# Bind the datasets into one data frame.
cpgi.probes.M.dat <- with(cpgi.probes.M,
	cbind(ALL[-ncol(ALL)], APL[-ncol(APL)], CTRL))
rm(cpgi.probes.M)

################################################################################
#### ALL vs. Control
################################################################################

# Take a subset of the data.
data <- t(subset(cpgi.probes.M.dat,
	select = grepl('HBC|ALL', colnames(cpgi.probes.M.dat))))
cpgi <- cpgi.probes.M.dat[, 'cpgi', drop=FALSE]
rm(cpgi.probes.M.dat)

# Specify the design matrix.
design <- data.frame(
	Group=relevel(ref='HBC', factor(substr(rownames(data), 1, 3))),
	row.names=rownames(data))

# Reshape the data.
tall <- melt(cbind(design, data), id.vars=colnames(design),
	variable.name = 'probe', value.name = 'M')
tall$cgi <- cpgi[tall$probe,'cpgi']

# Write the reshaped data to a file.
save(tall, file='Data/ALL_tall.Rdata')

# Fit a linear model.
system.time(
	lm.t <- ddply(tall, .(cgi), .parallel=TRUE, .progress='text',
		.fun=function(x)
			coef(summary(lm(M ~ Group, x)))['GroupALL', c('t value', 'Pr(>|t|)')]))
colnames(lm.t) <- c('cgi', 't', 'p')
lm.t$q <- p.adjust(lm.t$p, 'BH')

# Write the results to a file.
write.table(lm.t, 'Data/lm.tab')

# Write the gene set to a file.
lm.geneset <- coordToGene(subset(lm.t, subset=q<1e-5, select=cgi, drop=TRUE))
writeLines(lm.geneset, 'Data/lm-geneset.txt')

# Fit a linear mixed-effects model.
system.time(
	lme.t <- ddply(tall, .(cgi), .parallel=TRUE, .progress='text',
		.fun=function(x)
			coef(summary(lmer(M ~ Group + (1|probe), x)))['GroupALL', 't value']))
colnames(lme.t) <- c('cgi', 't')

# Write the results to a file.
write.table(lme.t_reml, 'Data/lme.tab')

# Another linear mixed-effects model, with 'nlme' package in order to apply ML
# and return p-values
### I found that with very large data sets, ddply gets very slow. After a quick test, I found that it looked to get exponentially slower with the size of the data set. So I will attempt to do this by chopping up the data frame into pieces (10, for ~2500 genes each) and running the same part on each
ids <- as.character(unique(cpgi[,1]))
ids_list <- split(ids, ceiling(1:length(ids)/2500))
res.list <- vector('list', length(ids_list))
for(i in 1:length(ids_list)){
  print(i)
  tall.sub <- subset(tall, cgi %in% ids_list[[i]])
  tall.sub <- droplevels(tall.sub)
  lme.t_ml <- ddply(tall.sub, .(cgi), .progress='text', .fun = function(x){
    res <- lme(M ~ Group, data = x, random = ~1 | probe, method = 'ML', na.action = 'na.omit')
    summary(res)$tTable['GroupALL', c('Value', 't-value', 'p-value')]
  })
  res.list[[i]] <- lme.t_ml
}
### And it turns out I am right... this was much faster. It took the same amount of time to run this than 10% of the whole data.

write.table(lme.t_ml, 'Data/lme_ml.tab')


# Another linear mixed-effects model, with 'nlme' package in order to apply ML
# and return p-values
### I found that with very large data sets, ddply gets very slow. After a quick test, I found that it looked to get exponentially slower with the size of the data set. So I will attempt to do this by chopping up the data frame into pieces (10, for ~2500 genes each) and running the same part on each
ids <- as.character(unique(cpgi[,1]))
ids_list <- split(ids, ceiling(1:length(ids)/2500))
res.list <- vector('list', length(ids_list))
for(i in 1:length(ids_list)){
  print(i)
  tall.sub <- subset(tall, cgi %in% ids_list[[i]])
  tall.sub <- droplevels(tall.sub)
  lme.t_ml <- ddply(tall.sub, .(cgi), .progress='text', .fun = function(x){
    res <- lme(M ~ Group, data = x, random = ~1 | probe, method = 'ML', na.action = 'na.omit')
    summary(res)$tTable['GroupALL', c('Value', 't-value', 'p-value')]
  })
  res.list[[i]] <- lme.t_ml
}


lme.t_ml <- do.call(rbind, res.list)
write.table(lme.t_ml, 'Data/lme_ml.tab')



# Write the gene set to a file.
lme.geneset <- coordToGene(
	subset(lme.t, subset = abs(t) > 10, select = cgi, drop = TRUE))
writeLines(lme.geneset, 'Data/lme-geneset.txt')
writeLines(con='Data/lme-down-geneset.txt',
	coordToGene(subset(lme.t, subset = t < -10, select = cgi, drop = TRUE)))
writeLines(con='Data/lme-up-geneset.txt',
	coordToGene(subset(lme.t, subset = t > 10, select = cgi, drop = TRUE)))
