# Aggregate the beta values of the probes for each CpG island.

#load('~/UBC Stats/STAT540/Group Project/Data/GSE42865_matrix.R')
cpginame <- as.data.frame(IlluminaHumanMethylation450kCPGINAME)
colnames(cpginame) <- c('Probe_ID', 'cpginame')
rownames(cpginame) <- cpginame$Probe_ID
cpginame$cpginame <- factor(cpginame$cpginame)
CTRL.cpginame <- merge(CTRL.dat, cpginame, by='row.names')
sampleNames <- colnames(CTRL.dat)
CTRL.cpgislands <- simplify2array(by(
	CTRL.cpginame[,sampleNames],
	list(CTRL.cpginame$cpginame),
	colMeans))

#load('~/UBC Stats/STAT540/Group Project/Data/GSE39141_matrix.R')
ALL.cpginame <- merge(ALL.dat, cpginame, by='row.names')
ALL.cpgislands <- simplify2array(by(
	ALL.cpginame[,colnames(ALL.dat)],
	list(ALL.cpginame$cpginame),
	colMeans))

# Combine the CTRL and ALL data sets.
CTRL.meta$Group <- 'CTRL'
design <- data.frame(
	Group = relevel(factor(c(CTRL.meta$Group, ALL.meta$Group)), 'CTRL'),
	row.names = c(row.names(CTRL.meta), row.names(ALL.meta)))
data <- rbind(CTRL.cpgislands, ALL.cpgislands)

# Fit a linear model.
library(limma)
mm <- model.matrix(~Group, design)
fit <- eBayes(lmFit(t(data), mm))
tt <- topTable(fit, coef='GroupALL', n=12)

# Reshape the data.
library(reshape2)
tall <- melt(cbind(design, data), id.vars=colnames(design),
	variable.name = 'probe', value.name = 'beta')

# Plot the hits.
library(lattice)
stripplot(beta ~ Group | probe, tall, subset = probe %in% tt$ID, auto.key=T)
