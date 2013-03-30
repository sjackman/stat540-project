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

# Bind the datasets into one data frame.
cpgi.probes.M.dat <- with(cpgi.probes.M,
	cbind(ALL[-ncol(ALL)], APL[-ncol(APL)], CTRL))
rm(cpgi.probes.M)

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

# Fit a linear model.
lm.t <- ddply(tall, .(cgi), .progress='text', .fun=function(x)
		coef(summary(lm(M ~ Group, x)))['GroupALL', c('t value', 'Pr(>|t|)')])
colnames(lm.t) <- c('cgi', 't', 'p')
lm.t$q <- p.adjust(lm.t$p, 'BH')

# Write the results to a file.
write.table(lm.t, 'Data/lm.tab')

# Write the gene set to a file.
lm.geneset <- coordToGene(subset(lm.t, subset=q<1e-5, select=cgi, drop=TRUE))))
writeLines(lm.geneset, 'Data/lm-geneset.txt')

# Fit a linear mixed-effects model.
lme.t <- ddply(tall, .(cgi), .progress='text', .fun=function(x)
	coef(summary(lmer(M ~ Group + (1|probe), x)))['GroupALL', 't value'])
colnames(lme.t) <- c('cgi', 't')

# Write the results to a file.
write.table(lme.t, 'Data/lme.tab')

# Write the gene set to a file.
lme.geneset <- coordToGene(
	subset(lme.t, subset = abs(t) > 10, select = cgi, drop = TRUE))
writeLines(lme.geneset, 'Data/lme-geneset.txt')
writeLines(con='Data/lme-down-geneset.txt',
	coordToGene(subset(lme.t, subset = t < -10, select = cgi, drop = TRUE)))
writeLines(con='Data/lme-up-geneset.txt',
	coordToGene(subset(lme.t, subset = t > 10, select = cgi, drop = TRUE)))

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
