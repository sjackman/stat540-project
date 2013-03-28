###########################################
### Script with examples for mixture models
##########################################
### Multiple versions of the same thing just to make sure it works!

## Mixture model packages
### nlme
### lme4
### CpGassoc

###
setwd('~/UBC Stats/STAT540/Group Project/')
load('Data/CPGI2Probe_betaList.Rdata')
load('Data/CPGI2Probe_MList.Rdata')
library(latticeExtra)
library(reshape2)
library(plyr)
### Bind data into one big ass data set: 
cpgi.probes.Beta.dat <- with(cpgi.probes.Beta, 
                             cbind(ALL[-ncol(ALL)], APL[-ncol(APL)], CTRL))
cpgi.probes.M.dat <- with(cpgi.probes.M, 
                             cbind(ALL[-ncol(ALL)], APL[-ncol(APL)], CTRL))

#### Dummy example of Mixed Effect Model (one Island):

### Pick one location to fit a linear mixed effect model on:
set.seed(127)
Loc1 <- sample(cpgi.probes.M.dat$cpgi, 1)
dummy <- subset(cpgi.probes.M.dat, cpgi == Loc1)

#### Specify groups for samples
Group <- c(rep('ALL', 29), rep('HBC', 4), rep('APL', 8), rep('HBM', 2), rep('HBC', 9))
Group.meta <- cbind(Sample = names(dummy)[-ncol(dummy)], Gp = Group)
dummy$Probe <- row.names(dummy)


#### Reshape data frame for easy formula call:
library(reshape2)
dummy.tall <- melt(dummy, id.vars = c('Probe', 'cpgi'), 
                   variable.name = 'Sample',
                   value.name = 'M')
dummy.tall <- merge(dummy.tall, Group.meta, by = 'Sample')
dummy.tall <- droplevels(dummy.tall)
dummy.tall$Probe <- factor(dummy.tall$Probe)

##########################################################################
#### Fitting linear mixed effects model with a random effect put on probe
##########################################################################

############################### nlme
library(nlme)
dummy.grouped <- groupedData(M ~ Gp|Probe, data = dummy.tall)

system.time(
  fit.nlme <- lme(M ~ Gp + Probe, data = dummy.grouped, 
                  random = ~Probe, na.action = na.omit)
)
  
## So slow!
# user  system elapsed 
# 368.57    0.59  377.67 


############################## lme4
library(lme4)
fit1.lme4 <- lmer(M ~ Gp + (1|Probe), dummy.tall)
fit2.lme4 <- lmer(M ~ Gp + (1|Probe) + (0+Gp|Probe), dummy.tall)
summary(fit1.lme4)
fit1.lme4

############################## CpGassoc
library(CpGassoc)
data(samplecpg,samplepheno,package="CpGassoc")
results <- cpg.assoc(samplecpg,samplepheno$weight,large.data=FALSE)
results


########### Plot to ensure results:
mycol <- brewer.pal(7, 'Set1')[c(1,5,2,7)]
mycol <- rgb(t(col2rgb(mycol)), alpha = 180, maxColorValue=255)
my.par <- list(superpose.symbol = list(col = mycol, pch = 16))
stripplot(M~Gp|Probe, data = dummy.tall, groups = Gp, par.settings = my.par,
          auto.key = list(space = 'top', columns = 4), jitter = T,
          layout = c(length(unique(dummy.tall$Probe)), 1))

#################################################
#### Differential Expression Analysis:
##############################################

### 1.) Determine groupings:


# CTRL.meta$Group <- 'CTRL'
# design <- data.frame(
#   Group = relevel(factor(c(CTRL.meta$Group, ALL.meta$Group)), 'CTRL'),
#   row.names = c(row.names(CTRL.meta), row.names(ALL.meta)))
# data <- rbind(CTRL.cpgislands, ALL.cpgislands)
# 
# # Fit a linear model.
# library(limma)
# mm <- model.matrix(~Group, design)
# fit <- eBayes(lmFit(t(data), mm))
# tt <- topTable(fit, coef='GroupALL', n=12)
# 
# # Reshape the data.
# library(reshape2)
# tall <- melt(cbind(design, data), id.vars=colnames(design),
#              variable.name = 'probe', value.name = 'beta')
# 
# # Plot the hits.
# library(lattice)
# stripplot(beta ~ Group | probe, tall, subset = probe %in% tt$ID, auto.key=T)