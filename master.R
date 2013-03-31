# Run the complete pipeline of analyses.
# Import the data, normalize, filter, perform hierarchical clustering,
# fit a linear mixed-effects model, and identify differentially methylated
# genes.

# Import the data from GEO
source('data_import.R')

# Perform some exploratory analysis of the raw data
source('exploratory.r')

# Normalize the data using BMIQ
source('normalize.r')

# Filter the data
source('exploratory_postNorm.r')

# Group probes into CpG islands for the raw data
source('aggregate_raw_norm_filter.R')

# Group probes into CpG islands for the normalized and filtered data
source('aggregate.R')

# Perform hierarchical clustering of the samples
source('hcluster.R')

# Fit a linear mixed-effects model, and identify differentially methylated
# CpG islands
source('lme.R')
