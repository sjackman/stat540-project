# Fit a linear mixed model

setwd('~/UBC Stats/STAT540/Group Project/')
load('Data/CPGI2Probe_MList.Rdata')
source('coord_to_gene.R')
library(latticeExtra)
library(lme4)
library(plyr)
library(RColorBrewer)
library(reshape2)
library(VennDiagram)
library(nlme)

# Bind the datasets into one data frame.
cpgi.probes.M.dat <- with(cpgi.probes.M,
	cbind(ALL[-ncol(ALL)], APL[-ncol(APL)], CTRL))
rm(cpgi.probes.M)

################################################################################
#### APL vs. Control
################################################################################
#Rename rows (APL = PLP + PLR)
substr(names(cpgi.probes.M.dat)[grepl('PLP|PLR', names(cpgi.probes.M.dat))], 1, 3) <-
  'APL'
# Take a subset of the data.
data <- t(subset(cpgi.probes.M.dat, 
                 select = grepl('HBC|APL', names(cpgi.probes.M.dat))))
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

# Fit a linear model.
lm.t <- ddply(tall, .(cgi), .progress='text', .fun=function(x)
		coef(summary(lm(M ~ Group, x)))['GroupALL', c('t value', 'Pr(>|t|)')])
colnames(lm.t) <- c('cgi', 't', 'p')
lm.t$q <- p.adjust(lm.t$p, 'BH')

# Write the results to a file.
write.table(lm.t, 'Data/lm.tab')

# Write the gene set to a file.
lm.geneset <- unique(coordToGene(as.character(
	subset(lm.t, subset=q<1e-5, select=cgi, drop=TRUE))))
writeLines(lm.geneset, 'Data/lm-geneset.txt')

lme_REML_topTable <- function(tall.data, coefName, size.each = 1500){
  ids <- as.character(unique(cpgi[,1]))
  ids_list <- split(ids, ceiling(1:length(ids)/size.each))
  res.list <- vector('list', length(ids_list))
  for(i in 1:length(ids_list)){
    print(i)
    tall.sub <- subset(tall.data, cgi %in% ids_list[[i]])
    tall.sub <- droplevels(tall.sub)
    lme.t_reml <- ddply(tall.sub, .(cgi), .progress='text', .fun=function(x)
      coef(summary(lmer(M ~ Group + (1|probe), x)))[coefName, c('Estimate','t value')])
    colnames(lme.t_reml) <- c('cgi', 'Estim', 't')
    res.list[[i]] <- lme.t_reml
  }
  return(do.call(rbind, res.list))
}

# Fit a linear mixed-effects model.
lme.t_reml <- ddply(tall, .(cgi), .progress='text', .fun=function(x)
  coef(summary(lmer(M ~ Group + (1|probe), x)))['GroupALL', 't value'])
colnames(lme.t_reml) <- c('cgi', 't')

# Write the results to a file.
write.table(lme.t_reml, 'Data/lme.tab')

# Another linear mixed-effects model, with 'nlme' package in order to apply ML
# and return p-values
### I found that with very large data sets, ddply gets very slow. After a quick test, I found that it looked to get exponentially slower with the size of the data set. So I will attempt to do this by chopping up the data frame into pieces (10, for ~2500 genes each) and running the same part on each
### Use this function: it will do the above. You need to specify coefName as it would appear in the result summary, e.g. "GroupALL" or "GroupAPL"
lme_ML_topTable <- function(tall.data, coefName){
  ids <- as.character(unique(cpgi[,1]))
  ids_list <- split(ids, ceiling(1:length(ids)/2500))
  res.list <- vector('list', length(ids_list))
  for(i in 1:length(ids_list)){
    print(i)
    tall.sub <- subset(tall.data, cgi %in% ids_list[[i]])
    tall.sub <- droplevels(tall.sub)
    lme.t_ml <- ddply(tall.sub, .(cgi), .progress='text', .fun = function(x){
      res <- lme(M ~ Group, data = x, random = ~1 | probe, method = 'ML', na.action = 'na.omit')
      summary(res)$tTable[coefName, c('Value', 't-value', 'p-value')]
    })
    res.list[[i]] <- lme.t_ml
  }
  return(do.call(rbind, res.list))
}

# Fit the LME
lme.t_ml_apl <- lme_ML_topTable(tall, 'GroupAPL')
# Write to table
write.table(lme.t_ml_apl, 'Data/lme_ml_apl.tab')


# Write the gene set to a file.
lme.geneset <- unique(coordToGene(as.character(
	subset(lme.t, subset = abs(t) > 10, select = cgi, drop = TRUE))))
writeLines(lme.geneset, 'Data/lme-geneset.txt')

# Plot a Venn diagram of the overlap of the fixed and mixed models.
plot.new()
grid.draw(venn.diagram(list(
	LinearModel = lm.geneset,
	LinearMixedEffectsModel = lme.geneset),
	filename=NULL, fill=c('red', 'blue')))

# Plot the denisty of the q-values of the linear model.
densityplot(lm.t$q,
	main='Density of q-values of the linear model',
	xlab='q-value')

# Plot the denisty of the t-statistic of the linear mixed-effects model.
densityplot(lme.t$t,
	main='Density of the t-statistic of the linear mixed-effects model',
	xlab='t-statistic')

# Plot the M-values of a CpG island.
stripplot(M ~ Group | probe, tall, group = Group,
	auto.key = TRUE, type = c('p', 'a'),
	subset = cgi == rownames(lme.t)[1])

# A similar plot.
mycol <- brewer.pal(7, 'Set1')[c(1,5,2,7)]
mycol <- rgb(t(col2rgb(mycol)), alpha = 180, maxColorValue=255)
my.par <- list(superpose.symbol = list(col = mycol, pch = 16))
stripplot(M ~ Group | probe, tall, groups = Group,
	par.settings = my.par,
	auto.key = TRUE, jitter = TRUE,
	layout = c(length(unique(tall$probe)), 1),
	subset = cgi == rownames(lme.t)[1])