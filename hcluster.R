###########################################################################################################
# This script contains clustering procedure for:
# (A)pre-normalized 
# (B)post-normalized
# (C)post-normalized and filtered by pairwise ctrl group comparison
###########################################################################################################

#####################################################
##### Clustering for pre-normalized raw data:
##### 1.) Probe-based (485577 probes for the CPG sites)
##### 2.) CPG-island-based (median or mean)
##### 3.) (TBA)Gene-based ? DM-genes-based
#####################################################

setwd('/ubc/cs/research/irmtraud/people/jingyun/localFile_backup/Documents/stat540/project/')
library(RColorBrewer)



##### Load raw data (modify the path accordingly)
load('Data/All_3_metasets.Rdata')
load('Data/All_3_sets.Rdata')

#####################################################
##### (1) probe level: clustering by Belta values #####

# merge 3 data sets and design
merged.dat <- cbind(APL.dat, ALL.dat, CTRL.dat)
merged.des <- data.frame(Samples=colnames(merged.dat), Group=c(APL.meta$Group, ALL.meta$Group, rep('HBC', 9)))

# compute pairwise distances
pr.dis <- dist(t(merged.dat), method = "euclidean")

# sanity check between the sample name and group assigned
pr.nms <- with(merged.des, paste(Samples, Group, sep = "_"))
all(unlist(lapply(pr.nms, function(x) substr(x, 1, 3)==substr(x, nchar(x)-2, nchar(x))))) # TRUE

 
# Explore which clustering method performs the best

# compute hierarchical clustering using different linkage types
pr.hc.s <- hclust(pr.dis, method = "single")
pr.hc.c <- hclust(pr.dis, method = "complete")
pr.hc.a <- hclust(pr.dis, method = "average")
pr.hc.w <- hclust(pr.dis, method = "ward")

# plot the dendrograms 

# (I) overview of the 4 clustering methods w/o sample labels
pdf("Figures/clusterOverviewByProbe_4Methods_raw.pdf")
par(mar = c(0, 4, 4, 2))
par(mfrow = c(2, 2))
plot(pr.hc.s, labels = FALSE, main = "Single", xlab = "")
 
plot(pr.hc.c, labels = FALSE, main = "Complete", xlab = "")
 
plot(pr.hc.a, labels = FALSE, main = "Average", xlab = "")
 
plot(pr.hc.w, labels = FALSE, main = "Ward", xlab = "")
dev.off()

# (II) plot for each individual method with sample labels to judge which one groups the samples the best

# ward --> FINAL CHOICE
pdf("Figures/clusterByProbe_ward_raw.pdf")
par(mfrow = c(1, 1))
# decide how many clusters we should cut
unique(merged.des$Group) # 5 categories: ALL HBC HBM PLP PLR 

# label the leaf nodes by coloring them according to their groups
ward.dendro <- as.dendrogram(pr.hc.w)
# color code each node of the dendrogram according to their group
colors <- brewer.pal(5, 'Set1')
names(colors) <- c("ALL", "HBC", "HBM", "PLP", "PLR")

colLab <- function(n) {
    if(is.leaf(n)) {
      a <- attributes(n)
      # find group name
      a.group <- substr(a$label, 1, 3)
      # retrieve the corresponding color
      attr(n, "nodePar") <-
        c(a$nodePar, list(lab.col = colors[a.group],lab.bg='grey50',pch=sample(19:25,1)))
      attr(n, "frame.plot") <- TRUE
    }
    n
  }
clusDendro <- dendrapply(ward.dendro, colLab)
plot(clusDendro, cex = 0.6, main = "Ward showing 5 clusters")
legend("topright", legend =c("ALL", "HBC", "HBM", "PLP", "PLR"), fill = colors, title="Sample Groups", box.col="transparent")

rect.hclust(pr.hc.w, k = 5, border=c("cyan"))
dev.off()

# single
pdf("Figures/clusterByProbe_single_raw.pdf")
par(mfrow = c(1, 1))
plot(pr.hc.s, labels = merged.des$Samples, cex = 0.6, main = "single showing 5 clusters")
rect.hclust(pr.hc.s, k = 5)
dev.off()

# complete  --> Not as good as Ward, tho better than single and average
pdf("Figures/clusterByProbe_complete_raw.pdf")
par(mfrow = c(1, 1))

# color leaf node of the dendrogram
comp.dendro <- as.dendrogram(pr.hc.c)
clusDendro <- dendrapply(comp.dendro, colLab)
plot(clusDendro, cex = 0.6, main = "complete showing 5 clusters")
legend("topright", legend =c("ALL", "HBC", "HBM", "PLP", "PLR"), fill = colors, title="Sample Groups", box.col="transparent")

#plot(pr.hc.c, labels = merged.des$Samples, cex = 0.6, main = "complete showing 5 clusters")
rect.hclust(pr.hc.c, k = 5, border=c("cyan"))
dev.off()

# Average
pdf("Figures/clusterByProbe_Average_raw.pdf")
par(mfrow = c(1, 1))
plot(pr.hc.a, labels = merged.des$Samples, cex = 0.6, main = "Average showing 5 clusters")
rect.hclust(pr.hc.a, k = 5)
dev.off()


# (III) headmap with ward dendrogram
pdf("Figures/heatmap_ClusterByProbe_raw.pdf")
# scale
sMerged.Dat <- t(scale(t(merged.dat)))
heatmap(as.matrix(sMerged.Dat), Rowv = NA, Colv = NULL, hclustfun = function(x) hclust(x, 
method = "ward"), distfun = function(x) dist(x, method = "euclidean"), scale = "none", 
labCol = pr.nms, labRow = NA, margin = c(8, 1), ColSideColor = brewer.pal(11, 
"RdGy")[c(4, 7)][unclass(merged.des$Group)])
dev.off()

# (IV) another overview plot for ward vs. complete using the scaled data
pdf("Figures/clusterOverviewByProbe_wad_com_scaledRaw.pdf")
par(mar = c(0, 4, 4, 2))
par(mfrow = c(2, 1))

# recalculate pairwise distance
pr.dis.s <- dist(t(sMerged.Dat), method = "euclidean")
pr.hc.c.s <- hclust(pr.dis.s, method = "complete")
pr.hc.w.s <- hclust(pr.dis.s, method = "ward")

plot(pr.hc.c.s, labels = merged.des$Samples, cex = 0.6, main = "complete showing 5 clusters (scaled)")
rect.hclust(pr.hc.c.s, k = 5)

plot(pr.hc.w.s, labels = merged.des$Samples, cex = 0.6, main = "Ward showing 5 clusters (scaled)")
rect.hclust(pr.hc.w.s, k = 5)
dev.off()


# (v) another ward plot using the scaled data: the ALL samples have a better clustering here
pdf("Figures/clusterByProbe_ward_scaledRaw.pdf")
par(mfrow = c(1, 1))

# label the leaf nodes by coloring them according to their groups
ward.dendro.s <- as.dendrogram(pr.hc.w.s)
clusDendro <- dendrapply(ward.dendro.s, colLab)
plot(clusDendro, cex = 0.6, main = "Scaled Ward showing 5 clusters")
legend("topright", legend =c("ALL", "HBC", "HBM", "PLP", "PLR"), fill = colors, title="Sample Groups", box.col="transparent")

rect.hclust(pr.hc.w.s, k = 5, border=c("cyan"))
dev.off()

#####################################################
##### (2) cpgi level: clustering by Belta values #####

##### Load the raw probe data
#load('Data/CPGI2Probe_betaList_raw.Rdata')

### (I) ward plot for the mean cpgi belta values

# construct a data frame with rows = cpgi IDs, columns = mean belta value for each sample
load('Data/CPGI_betaMeanList_raw.Rdata')
#### Transpose to place the cpgi features in rows
cpgi.Beta.mean <- lapply(cpgi.Beta.mean, t)
cpgi.Beta.mean.whole <- do.call(cbind, cpgi.Beta.mean)
dim(cpgi.Beta.mean.whole)
# [1] 27176    52


# ward cluster based on the mean cpgi belta values
# compute pairwise distances
pr.dis.mean <- dist(t(cpgi.Beta.mean.whole), method = "euclidean")
pr.hc.w.mean <- hclust(pr.dis.mean, method = "ward")

pdf("Figures/clusterByMeanCPGI_ward_Raw.pdf")
par(mfrow = c(1, 1))

# label the leaf nodes by coloring them according to their groups
mean.ward.dendro <- as.dendrogram(pr.hc.w.mean)
clusDendro <- dendrapply(mean.ward.dendro, colLab)
plot(clusDendro, cex = 0.6, main = "CPGI-based Ward(mean) Clustering showing 5 clusters")
legend("topright", legend =c("ALL", "HBC", "HBM", "PLP", "PLR"), fill = colors, title="Sample Groups", box.col="transparent")

rect.hclust(pr.hc.w.mean, k = 5, border=c("cyan"))
dev.off()


# (II) ward plot for the median cpgi belta values
load('Data/CPGI_betaMedianList_raw.Rdata')
#### Transpose to place the cpgi features in rows
cpgi.Beta.median <- lapply(cpgi.Beta.median, t)
cpgi.Beta.median.whole <- do.call(cbind, cpgi.Beta.median)
dim(cpgi.Beta.median.whole)
# [1] 27176    52


# ward cluster based on the median cpgi belta values
# compute pairwise distances
pr.dis.median <- dist(t(cpgi.Beta.median.whole), method = "euclidean")
pr.hc.w.median <- hclust(pr.dis.median, method = "ward")

pdf("Figures/clusterByMedianCPGI_ward_Raw.pdf")
par(mfrow = c(1, 1))

# label the leaf nodes by coloring them according to their groups
median.ward.dendro <- as.dendrogram(pr.hc.w.median)
clusDendro <- dendrapply(median.ward.dendro, colLab)
plot(clusDendro, cex = 0.6, main = "CPGI-based Ward(median) Clustering showing 5 clusters")
legend("topright", legend =c("ALL", "HBC", "HBM", "PLP", "PLR"), fill = colors, title="Sample Groups", box.col="transparent")

rect.hclust(pr.hc.w.median, k = 5, border=c("cyan"))
dev.off()


#####################################################
##### Clustering for post-normalized data:
##### 1.) Probe-based (485577 probes for the CPG sites)
##### 2.) CPG-island-based
##### 3.) Gene-based ?DM-genes-based
#####################################################
# load the data
load('Data/All_3_sets_normalized.Rdata')
dim(ALL.norm)
# [1] 485577     33



#####################################################
##### (1) probe level: clustering by Belta values #####

# compute pairwise distances
merged.dat <- cbind(APL.norm, ALL.norm, CTRL.norm)
pr.dis <- dist(t(merged.dat), method = "euclidean")

# sanity check between the sample name and group assigned
merged.des <- data.frame(Samples=colnames(merged.dat), Group=c(APL.meta$Group, ALL.meta$Group, rep('HBC', 9)))
pr.nms <- with(merged.des, paste(Samples, Group, sep = "_"))
all(unlist(lapply(pr.nms, function(x) substr(x, 1, 3)==substr(x, nchar(x)-2, nchar(x))))) # TRUE

# compute hierarchical clustering using ward linkage type
pr.hc.w <- hclust(pr.dis, method = "ward")


# ward 
pdf("Figures/clusterByProbe_ward_norm.pdf")
par(mfrow = c(1, 1))

# label the leaf nodes by coloring them according to their groups
ward.dendro <- as.dendrogram(pr.hc.w)
# color code each node of the dendrogram according to their group
clusDendro <- dendrapply(ward.dendro, colLab)
plot(clusDendro, cex = 0.6, main = "Ward showing 5 clusters")
legend("topright", legend =c("ALL", "HBC", "HBM", "PLP", "PLR"), fill = colors, title="Sample Groups", box.col="transparent")

rect.hclust(pr.hc.w, k = 5, border=c("cyan"))
dev.off()




#####################################################
##### (2) cpgi level: clustering by Belta values #####
##### Load the raw probe data
#load('Data/CPGI2Probe_betaList_raw.Rdata')

### (I) ward plot for the mean cpgi belta values

# construct a data frame with rows = cpgi IDs, columns = mean belta value for each sample
load('Data/CPGI_betaMeanList_norm.Rdata')
dim(cpgi.Beta.mean$ALL)
# [1]    33 27176



#### Transpose to place the cpgi features in rows
cpgi.Beta.mean <- lapply(cpgi.Beta.mean, t)
cpgi.Beta.mean.whole <- do.call(cbind, cpgi.Beta.mean)
dim(cpgi.Beta.mean.whole)
# [1] 27176    52


# ward cluster based on the mean cpgi belta values
# compute pairwise distances
pr.dis.mean <- dist(t(cpgi.Beta.mean.whole), method = "euclidean")
pr.hc.w.mean <- hclust(pr.dis.mean, method = "ward")

pdf("Figures/clusterByMeanCPGI_ward_Norm.pdf")
par(mfrow = c(1, 1))

# label the leaf nodes by coloring them according to their groups
mean.ward.dendro <- as.dendrogram(pr.hc.w.mean)
clusDendro <- dendrapply(mean.ward.dendro, colLab)
plot(clusDendro, cex = 0.6, main = "CPGI-based Ward(mean) Clustering showing 5 clusters")
legend("topright", legend =c("ALL", "HBC", "HBM", "PLP", "PLR"), fill = colors, title="Sample Groups", box.col="transparent")

rect.hclust(pr.hc.w.mean, k = 5, border=c("cyan"))
dev.off()


# (II) ward plot for the median cpgi belta values
load('Data/CPGI_betaMedianList_raw.Rdata')
#### Transpose to place the cpgi features in rows
cpgi.Beta.median <- lapply(cpgi.Beta.median, t)
cpgi.Beta.median.whole <- do.call(cbind, cpgi.Beta.median)
dim(cpgi.Beta.median.whole)
# [1] 27176    52


# ward cluster based on the median cpgi belta values
# compute pairwise distances
pr.dis.median <- dist(t(cpgi.Beta.median.whole), method = "euclidean")
pr.hc.w.median <- hclust(pr.dis.median, method = "ward")

pdf("Figures/clusterByMedianCPGI_ward_Norm.pdf")
par(mfrow = c(1, 1))

# label the leaf nodes by coloring them according to their groups
median.ward.dendro <- as.dendrogram(pr.hc.w.median)
clusDendro <- dendrapply(median.ward.dendro, colLab)
plot(clusDendro, cex = 0.6, main = "CPGI-based Ward(median) Clustering showing 5 clusters")
legend("topright", legend =c("ALL", "HBC", "HBM", "PLP", "PLR"), fill = colors, title="Sample Groups", box.col="transparent")

rect.hclust(pr.hc.w.median, k = 5, border=c("cyan"))
dev.off()


#####################################################
##### Clustering for post-normalized and filtered data:
##### 1.) Probe-based (435739 probes for the CPG sites)
##### 2.) CPG-island-based
##### 3.) Gene-based ?DM-genes-based
#####################################################
# load the data
load('Data/All_3_sets_normAndFilt.Rdata')
dim(ALL.dat)
# [1] 435739     33


#####################################################
##### (1) probe level: clustering by Belta values #####

# compute pairwise distances
merged.dat <- cbind(APL.dat, ALL.dat, CTRL.dat)
pr.dis <- dist(t(merged.dat), method = "euclidean")

# sanity check between the sample name and group assigned
merged.des <- data.frame(Samples=colnames(merged.dat), Group=c(APL.meta$Group, ALL.meta$Group, rep('HBC', 9)))
pr.nms <- with(merged.des, paste(Samples, Group, sep = "_"))
all(unlist(lapply(pr.nms, function(x) substr(x, 1, 3)==substr(x, nchar(x)-2, nchar(x))))) # TRUE

# compute hierarchical clustering using ward linkage type
pr.hc.w <- hclust(pr.dis, method = "ward")


# ward 
pdf("Figures/clusterByProbe_ward_normFilter.pdf")
par(mfrow = c(1, 1))

# label the leaf nodes by coloring them according to their groups
ward.dendro <- as.dendrogram(pr.hc.w)
# color code each node of the dendrogram according to their group
clusDendro <- dendrapply(ward.dendro, colLab)
plot(clusDendro, cex = 0.6, main = "Ward showing 5 clusters")
legend("topright", legend =c("ALL", "HBC", "HBM", "PLP", "PLR"), fill = colors, title="Sample Groups", box.col="transparent")

rect.hclust(pr.hc.w, k = 5, border=c("cyan"))
dev.off()




##################################################### IGNORE THE CODE BELOW. THERE IS A BUG, I WILL FIX IT LATER
##### (2) cpgi level: clustering by Belta values #####
##### Load the raw probe data
#load('Data/CPGI2Probe_betaList_raw.Rdata')

### (I) ward plot for the mean cpgi belta values

# construct a data frame with rows = cpgi IDs, columns = mean belta value for each sample
load('Data/CPGI_betaMeanList_normFilt.Rdata')
length(cpgi.Beta.mean$ALL)
[1] 27176


cpgi.Beta.mean.whole <- do.call(cbind, cpgi.Beta.mean)
dim(cpgi.Beta.mean.whole)
# [1] 27176    3

library (plyr)
mean.whole.expand <- cbind(ldply(cpgi.Beta.mean.whole[1, 1], unlist), ldply(cpgi.Beta.mean.whole[1, 2], unlist), ldply(cpgi.Beta.mean.whole[1, 3], unlist))
for(i in 2: dim(cpgi.Beta.mean.whole)[1]) {
 print(i)
 new.row <- cbind(ldply(cpgi.Beta.mean.whole[i, 1], unlist), ldply(cpgi.Beta.mean.whole[i, 2], unlist), ldply(cpgi.Beta.mean.whole[i, 3], unlist))
 mean.whole.expand <- rbind(mean.whole.expand, new.row)
}

mean.whole.expand <- data.frame(mean.whole.expand)
rownames(mean.whole.expand) <- rownames(cpgi.Beta.mean.whole)


# ward cluster based on the mean cpgi belta values
# compute pairwise distances
pr.dis.mean <- dist(t(cpgi.Beta.mean.whole), method = "euclidean")
pr.hc.w.mean <- hclust(pr.dis.mean, method = "ward")

pdf("Figures/clusterByMeanCPGI_ward_NormFilt.pdf")
par(mfrow = c(1, 1))

# label the leaf nodes by coloring them according to their groups
mean.ward.dendro <- as.dendrogram(pr.hc.w.mean)
clusDendro <- dendrapply(mean.ward.dendro, colLab)
plot(clusDendro, cex = 0.6, main = "CPGI-based Ward(mean) Clustering showing 5 clusters")
legend("topright", legend =c("ALL", "HBC", "HBM", "PLP", "PLR"), fill = colors, title="Sample Groups", box.col="transparent")

rect.hclust(pr.hc.w.mean, k = 5, border=c("cyan"))
dev.off()


# (II) ward plot for the median cpgi belta values
load('Data/CPGI_betaMedianList_raw.Rdata')
#### Transpose to place the cpgi features in rows
cpgi.Beta.median <- lapply(cpgi.Beta.median, t)
cpgi.Beta.median.whole <- do.call(cbind, cpgi.Beta.median)
dim(cpgi.Beta.median.whole)
# [1] 27176    52


# ward cluster based on the median cpgi belta values
# compute pairwise distances
pr.dis.median <- dist(t(cpgi.Beta.median.whole), method = "euclidean")
pr.hc.w.median <- hclust(pr.dis.median, method = "ward")

pdf("Figures/clusterByMedianCPGI_ward_NormFilt.pdf")
par(mfrow = c(1, 1))

# label the leaf nodes by coloring them according to their groups
median.ward.dendro <- as.dendrogram(pr.hc.w.median)
clusDendro <- dendrapply(median.ward.dendro, colLab)
plot(clusDendro, cex = 0.6, main = "CPGI-based Ward(median) Clustering showing 5 clusters")
legend("topright", legend =c("ALL", "HBC", "HBM", "PLP", "PLR"), fill = colors, title="Sample Groups", box.col="transparent")

rect.hclust(pr.hc.w.median, k = 5, border=c("cyan"))
dev.off()
