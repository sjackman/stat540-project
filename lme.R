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

# Specify the design matrix.
data <- t(subset(cpgi.probes.M.dat, select=-cpgi))
design <- data.frame(
	Group=relevel(ref='HBC', factor(substr(rownames(data), 1, 3))),
	row.names=rownames(data))

# Reshape the data.
tall <- melt(cbind(design, data), id.vars=colnames(design),
	variable.name = 'probe', value.name = 'M')
tall$cgi <- cpgi.probes.M.dat[tall$probe,'cpgi']

# Extract a small subset of the CpG islands for two groups.
small.cgi <- head(names(table(tall$cgi)), n=12)
small <- droplevels(subset(tall,
	subset = cgi %in% small.cgi & Group %in% c('HBC', 'ALL')))

# Fit a linear model.
lm.fit <- dlply(small, .(cgi), .progress='text',
	.fun=function(x) lm(M ~ Group, x))
lm.t <- t(sapply(lm.fit,
	function(x) coef(summary(x))['GroupALL', c('t value', 'Pr(>|t|)')]))
colnames(lm.t) <- c('t', 'p')

# Fit a linear mixed-effects model.
lme.fit <- dlply(small, .(cgi), .progress='text',
	.fun=function(x) lmer(M ~ Group + (1|probe), x))
lme.t <- data.frame(t=sapply(lme.fit,
	function(x) coef(summary(x))['GroupALL', 't value']))

# Plot the M-values of a CpG island.
stripplot(M ~ Group | probe, small, group = Group,
	auto.key = TRUE, type = c('p', 'a'),
	subset = cgi == rownames(lme.t)[1])

# A similar plot.
mycol <- brewer.pal(7, 'Set1')[c(1,5,2,7)]
mycol <- rgb(t(col2rgb(mycol)), alpha = 180, maxColorValue=255)
my.par <- list(superpose.symbol = list(col = mycol, pch = 16))
stripplot(M ~ Group | probe, small, groups = Group,
	par.settings = my.par,
	auto.key = TRUE, jitter = TRUE,
	layout = c(length(unique(small$probe)), 1),
	subset = cgi == rownames(lme.t)[1])
