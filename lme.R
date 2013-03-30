# Fit a linear mixed model

setwd('~/UBC Stats/STAT540/Group Project/')
load('Data/CPGI2Probe_MList.Rdata')
library(latticeExtra)
library(lme4)
library(plyr)
library(RColorBrewer)
library(reshape2)

# Bind the datasets into one data frame.
cpgi.probes.M.dat <- with(cpgi.probes.M,
	cbind(ALL[-ncol(ALL)], APL[-ncol(APL)], CTRL))
rm(cpgi.probes.M)

# Take a subset of the data.
data <- t(subset(cpgi.probes.M.dat, select = grepl('HBC|ALL', rownames(data))))
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

# Write the gene set to a file.
source('coord_to_gene.R')
geneset <- coordToGene(as.character(subset(lm.t, subset=q<1e-5, select=cgi, drop=TRUE)))
writeLines(geneset, 'Data/lm-geneset.txt')

# Fit a linear mixed-effects model.
lme.t <- dlply(tall, .(cgi), .progress='text', .fun=function(x)
	coef(summary(lmer(M ~ Group + (1|probe), x)))['GroupALL', 't value'])
colnames(lme.t) <- c('cgi', 't')

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
