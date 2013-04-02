# Run the complete pipeline of analyses.
# Import the data, normalize, filter, perform hierarchical clustering,
# fit a linear mixed-effects model, and identify differentially methylated
# genes.

# Import the data from GEO
source('data_import.R')

# Perform exploratory analysis of the raw data
source('exploratory.r')

# Normalize the data using BMIQ
source('normalize.r')

# Filter the data
source('exploratory_postNorm.r')

# Group probes into CpG islands before and after normalizing and
# filtering the data
source('aggregate_raw_norm_filter.R')

# Group probes into CpG islands after normalizing and filtering the data
source('aggregate.R')

# Perform hierarchical clustering of the samples
source('hcluster.R')

# Fit a linear mixed-effects model, and identify differentially methylated
# CpG islands for ALL vs control
source('lme_all.R')

# Fit a linear mixed-effects model, and identify differentially methylated
# CpG islands for APL vs control
source('lme_plp.R')

# Produce a heat map of probes in top 10 Islands and a Venn diagram of
# overlap of top islands
source('differential_methylation_figures.R')

# Identify enriched Gene Ontology terms
source('topGO.R')
