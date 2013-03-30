# Convert a set of coordinates to a set of gene names

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

coordToGene <- function(x) {
# Convert the coordinates of CpG islands to gene names.
	cgi <- data.frame(t(simplify2array(strsplit(x, '[:-]'))),
		row.names=x,
		stringsAsFactors=FALSE)
	colnames(cgi) <- c('chr', 'start', 'end')
	cgi$chr <- factor(sub('chr', '', cgi$chr))
	cgi$start <- as.numeric(cgi$start)
	cgi$end <- as.numeric(cgi$end)
	cgi.gr <- GRanges(
		seqnames = Rle(cgi$chr),
		ranges = IRanges(start=cgi$start, end=cgi$end),
		names = rownames(cgi))
	overlaps <- findOverlaps(cgi.gr, genes.gr)
	return(genes[subjectHits(overlaps), 'hgnc_symbol'])
}
