###########################################################################################################
# This script contains clustering with AU&BP p-values procedure for:
# (A)pre-normalized 
# (B)post-normalized
# (C)post-normalized and filtered by pairwise ctrl group comparison
###########################################################################################################
setwd('/ubc/cs/research/irmtraud/people/jingyun/localFile_backup/Documents/stat540/project/')
library(RColorBrewer)
library(pvclust)



#####################################################
#*** Generate the pvclustering result by running the 
# following section of code on a cluster server 
# DON'T try run them as they take a whole day to finish
#####################################################

# (i) Clustering for pre-normalized raw data:
setwd('/ubc/cs/research/irmtraud/people/jingyun/localFile_backup/Documents/stat540/project/')

##### Load raw data (modify the path accordingly)
load('/ubc/cs/research/irmtraud/people/jingyun/localFile_backup/Documents/stat540/project/Data/All_3_sets.Rdata')

# merge 3 data sets and design
merged.dat <- cbind(APL.dat, ALL.dat, CTRL.dat)

# plot clustering using ward method
# cluster with bootstrapping
result.raw <- pvclust(merged.dat, method.dist="correlation", method.hclust="ward", nboot=100)
save(result.raw, file = '/ubc/cs/research/irmtraud/people/jingyun/localFile_backup/Documents/stat540/project/Data/result_raw.Rdata')



# (ii) Clustering for post-normalized raw data:
setwd('/ubc/cs/research/irmtraud/people/jingyun/localFile_backup/Documents/stat540/project/')
##### Load normalized data (modify the path accordingly)
load('/ubc/cs/research/irmtraud/people/jingyun/localFile_backup/Documents/stat540/project/Data/All_3_sets_normalized.Rdata')

# merge 3 data sets and design
merged.dat <- cbind(APL.norm, ALL.norm, CTRL.norm)

# plot clustering using ward method
# cluster with bootstrapping
result.norm <- pvclust(merged.dat, method.dist="correlation", method.hclust="ward", nboot=100)
save(result.norm, file = '/ubc/cs/research/irmtraud/people/jingyun/localFile_backup/Documents/stat540/project/Data/result_norm.Rdata')


# (iii) Clustering for post-normalized & filtered raw data:
setwd('/ubc/cs/research/irmtraud/people/jingyun/localFile_backup/Documents/stat540/project/')

##### Load raw data (modify the path accordingly)
load('/ubc/cs/research/irmtraud/people/jingyun/localFile_backup/Documents/stat540/project/Data/All_3_sets_normAndFilt.Rdata')

# merge 3 data sets and design
merged.dat <- cbind(APL.dat, ALL.dat, CTRL.dat)

# plot clustering using ward method
# cluster with bootstrapping
result.normFilt <- pvclust(merged.dat, method.dist="correlation", method.hclust="ward", nboot=100)
save(result.normFilt, file = '/ubc/cs/research/irmtraud/people/jingyun/localFile_backup/Documents/stat540/project/Data/result_normFilt.Rdata')


#####################################################
#*** Clustering result for pre-normalized raw data:
#*** 1.) Probe-based (485577 probes for the CPG sites)
#####################################################

### load result from clustering with bootstrapping
load('Data/result_raw.Rdata')

# plot clustering using ward method
pdf("Figures/pvclusterByProbe_raw.pdf", width=30, height=20)
par(mfrow = c(2, 1))

# label the leaf nodes by coloring them according to their groups
result_dendrogram <- as.dendrogram(result.raw$hclust)
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

# show leaf nodes colors
clusDendro <- dendrapply(result_dendrogram, colLab)
plot(clusDendro,  main = "Ward Clustering showing 5 clusters")
legend("topright", legend =c("ALL", "HBC", "HBM", "PLP", "PLR"), fill = colors, title="Sample Groups", box.col="transparent")
rect.hclust(result.raw$hclust, k = 5, border=c("cyan"))



# show AU & BP values
plot(result.raw,  main = "Ward Clustering for raw data with AU/Bp p-values")
rect.hclust(result.raw$hclust, k = 5, border=c("cyan"))
dev.off()

#*********************************************************************************

#####################################################
#*** Clustering result for post-normalized norm data:
#*** 1.) Probe-based (485577 probes for the CPG sites)
#####################################################

### load result from clustering with bootstrapping
load('Data/result_norm.Rdata')

# plot clustering using ward method
pdf("Figures/pvclusterByProbe_norm.pdf", width=30, height=20)
par(mfrow = c(2, 1))

# label the leaf nodes by coloring them according to their groups
result_dendrogram <- as.dendrogram(result.norm$hclust)
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

# show leaf nodes colors
clusDendro <- dendrapply(result_dendrogram, colLab)
plot(clusDendro,  main = "Ward Clustering showing 5 clusters")
legend("topright", legend =c("ALL", "HBC", "HBM", "PLP", "PLR"), fill = colors, title="Sample Groups", box.col="transparent")
rect.hclust(result.norm$hclust, k = 5, border=c("cyan"))



# show AU & BP values
plot(result.norm,  main = "Ward Clustering for normalized data with AU/Bp p-values")
rect.hclust(result.norm$hclust, k = 5, border=c("cyan"))
dev.off()

#*********************************************************************************

#####################################################
#*** Clustering result for pre-normalized and Filtered data:
#*** 1.) Probe-based (435739 probes for the CPG sites)
#####################################################

### load result from clustering with bootstrapping
load('Data/result_normFilt.Rdata')

# plot clustering using ward method
pdf("Figures/pvclusterByProbe_normFilt.pdf", width=30, height=20)
par(mfrow = c(2, 1))

# label the leaf nodes by coloring them according to their groups
result_dendrogram <- as.dendrogram(result.normFilt$hclust)
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

# show leaf nodes colors
clusDendro <- dendrapply(result_dendrogram, colLab)
plot(clusDendro,  main = "Ward Clustering showing 5 clusters")
legend("topright", legend =c("ALL", "HBC", "HBM", "PLP", "PLR"), fill = colors, title="Sample Groups", box.col="transparent")
rect.hclust(result.normFilt$hclust, k = 5, border=c("cyan"))



# show AU & BP values
plot(result.normFilt,  main = "Ward Clustering for normalized and filtered with AU/Bp p-values")
rect.hclust(result.normFilt$hclust, k = 5, border=c("cyan"))
dev.off()

#*********************************************************************************
