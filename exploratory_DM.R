#### Hits analysis

#### Take lists from all 3 analyses
### --> q < 1e-25
### 1.) Venn Diagram
### 

setwd('~/UBC Stats/STAT540/Group Project/')

library(plyr)
library(reshape2)
library(latticeExtra)
library(VennDiagram)

### top-tables from LME (with p-values)
topML.all <- read.table('Data/lme.t_ml_all.tab')
topML.plp <- read.table('Data/lme.t_ml_plp.tab')
topML.apl <- read.table('Data/lme.t_ml_apl.tab')
### top-tables from REML (lmer, no p-values)
topRE.all <- read.table('Data/lme.t_reml_all.tab')
topRE.apl <- read.table('Data/lme.t_reml_apl.tab')

###################################################################
### Compare MLE and REML in APL vs. control (most small q-values)
###################################################################
topBoth.apl <- merge(topML.apl, topRE.apl, by = 'cgi')
head(topBoth.apl)
with(topBoth.apl, table(Value - Estim < 1e-2))

### Difference in t-statistics:
with(topBoth.apl, plot(density(t.value - t), main = 'Difference in t-statistics between LME and REML per gene in APL vs. control', xlab = 'Difference'))

## Relationship between t-difference and q-value:
with(topBoth.apl, plot(x = t.value - t, y = q.val))

## Correlation:
xyplot(t.value ~ t, data = topBoth.apl)
with(topBoth.apl, cor(t.value, t, method = 'spearman'))
with(topBoth.apl, max(abs(t.value - t))) ## 0.223889 (something to worry about... I think not)


##########################################################################
### Compare PLP and APL
##########################################################################

# Merge into on data frame
topPL <- cbind(topML.plp, topML.apl[-1])
names(topPL)[2:5] <- jPaste(names(topPL)[2:5], 'PLP')
names(topPL)[6:9] <- jPaste(names(topPL)[6:9], 'APL')


# Overlap of "significant" hits
### With approximately equal-sized sets:
venn.list <- with(topPL, list(APL = which(q.valAPL < 1e-25), PLP = which(q.valPLP < 1e-18)))
plot.new()
grid.draw(venn.diagram(venn.list, filename = NULL))

## With the same q-cutoff:
plot.new()
grid.draw(venn.diagram(lapply(topPL[c('q.valPLP', 'q.valAPL')], function(x) which(x<1e-25)), filename = NULL))


#########################################################################
#### Strip plot of top hits
#########################################################################

### Get number of probes per island
load('Data/CPGI2Probe_MList.Rdata')
cpgi.lth <- cpgi.probes.M[[1]]
nopi <- with(cpgi.lth, tapply(ALL_956761, cpgi, length))

### Merge results from PLP and ALL analyses
topBoth <- cbind(topML.plp, topML.all[-1])
names(topBoth)[2:5] <- jPaste(names(topBoth)[2:5], 'PLP')
names(topBoth)[6:9] <- jPaste(names(topBoth)[6:9], 'ALL')

topBoth$nopi <- nopi[as.character(topBoth$cgi)]

### Now choose those islands with 5 probes
chooseExample <- subset(topBoth, nopi == 5)

### And choose interesting islands 
### (these cutoffs were finagled to get the most exemplary gene per group)
whichIx <- with(chooseExample,
                      c(which.min(q.valPLP * q.valALL),
                        which(q.valPLP > 0.05 & q.valALL < 1e-15),
                        which(q.valPLP < 1e-24 & q.valALL > 0.05),
                        which.max(q.valPLP*q.valALL)))

cgi.examples <- as.character(chooseExample[whichIx, 'cgi'])

## good.. now take subset of M values for these islands and convert to tall:
cpgi.probes.M.dat <- with(cpgi.probes.M,
                          cbind(ALL[-ncol(ALL)], APL[-ncol(APL)], CTRL))
data <- subset(cpgi.probes.M.dat, cpgi %in% cgi.examples,
               select = c(grepl('HBC|ALL|PLP|cpgi', names(cpgi.probes.M.dat))))
data$probe <- rownames(data)
tall <- melt(data, id.vars = c('cpgi', 'probe'), value.name = 'M', variable.name = 'Sample')
tall$Group <- substr(tall$Sample, 1, 3)
tall <- droplevels(tall)

tall$probe <- factor(tall$probe)
tall$Group <- factor(tall$Group)

### Strip Plots
#my.par <- list(superpose.symbol = list(col = rgb(0, 0, 1, 0.7), pch = 16))

for(i in 1:length(levels(tall$cpgi))){
  CPGI <- levels(tall$cpgi)[i]
  #pdf(jPaste('Example', i, '.pdf'), width = 10, height = 6)
  tall.sub <- subset(tall, cpgi == CPGI)
  tall.sub <- droplevels(tall.sub)
  title <- paste(CPGI, paste(c('q.PLP', 'q.ALL'), '=', format(subset(chooseExample, cgi == CPGI)[1, c('q.valPLP', 'q.valALL')], digits = 4), collapse = ', '))
  print(stripplot(M ~ Group | probe, tall.sub, col = rgb(0, 0, 1, 0.6), pch = 16,
          type = c('p', 'a'), jitter = TRUE,
          layout = c(5, 1), main = title))
  #dev.off()
}

