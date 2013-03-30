# Identify differentially methylated CpG islands using limma

library(lattice)
library(limma)
library(reshape2)

# Load the CpG island M values.
load('Data/CPGI_MList.Rdata')

# Combine the CTRL and ALL data sets.
data <- rbind(
	do.call(cbind, cpgi.M[['CTRL']]),
	do.call(cbind, cpgi.M[['ALL']]))
design <- data.frame(
	Group=relevel(ref='HBC', factor(substr(rownames(data), 1, 3))),
	row.names=rownames(data))

# Fit a linear model.
mm <- model.matrix(~Group, design)
fit <- eBayes(lmFit(t(data), mm))
tt <- topTable(fit, coef='GroupALL', number=Inf, p.value=1e-10)

# Reshape the data.
tall <- melt(cbind(design, data), id.vars=colnames(design),
  variable.name = 'probe', value.name = 'M')

# Plot the hits.
stripplot(M ~ Group | probe, tall,
	subset = probe %in% tt[1:12,'ID'], auto.key=T)

# Write the gene set to a file.
source('coord_to_gene.R')
writeLines(coordToGene(tt$ID), 'Data/limma-geneset.txt')
