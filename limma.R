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

# Get the coordinates of genes.
library(biomaRt)
library(GenomicRanges)
genes <- subset(subset = hgnc_symbol != '',
	getBM(
		c('hgnc_symbol', 'ensembl_gene_id', 'chromosome_name', 'start_position', 'end_position'),
		mart=useMart('ensembl', 'hsapiens_gene_ensembl')))
genes.gr <- GRanges(
	seqnames = Rle(genes$chromosome_name),
	ranges = IRanges(start = genes$start_position, end = genes$end_position),
	names = genes$hgnc_symbol)

# Convert the coordinates of CpG islands to gene names.
cpgi <- data.frame(t(simplify2array(strsplit(tt$ID, '[:-]'))),
	row.names=tt$ID,
	stringsAsFactors=FALSE)
colnames(cpgi) <- c('chr', 'start', 'end')
cpgi$chr <- factor(sub('chr', '', cpgi$chr))
cpgi$start <- as.numeric(cpgi$start)
cpgi$end <- as.numeric(cpgi$end)
cpgi.gr <- GRanges(
	seqnames = Rle(cpgi$chr),
	ranges = IRanges(start=cpgi$start, end=cpgi$end),
	names = rownames(cpgi))
overlaps <- findOverlaps(cpgi.gr, genes.gr)
geneset <- genes[subjectHits(overlaps), 'hgnc_symbol']

# Write the gene set to a file.
fd <- file('Data/limma-geneset.txt')
writeLines(geneset, fd)
close(fd)
