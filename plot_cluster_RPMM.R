setwd('/ubc/cs/research/irmtraud/people/jingyun/localFile_backup/Documents/stat540/project/')
library(RColorBrewer)
library(RPMM)


#######################################################################################
#####  RAW DATA
#####
#######################################################################################



##### Load raw data and the mmfit-fited result(modify the path accordingly)
load('/ubc/cs/research/irmtraud/people/jingyun/localFile_backup/Documents/stat540/project/Data/mmfit_raw.Rdata')
load('/ubc/cs/research/irmtraud/people/jingyun/localFile_backup/Documents/stat540/project/Data/All_3_sets.Rdata')

# merge 3 data sets and design
merged.dat <- cbind(APL.dat, ALL.dat, CTRL.dat)
sample.names <- colnames(merged.dat)


# find how many groups resulted from the RPMM clustering 
class.raw <- blcTreeLeafMatrix(mmfit.raw)
dim(class.raw) # 9 groups for these 52 samples. 

# record the membership of the samples for the 9 groups
membership <- data.frame(class=colnames(class.raw),  members=vector(mode="character", length=length(colnames(class.raw))), stringsAsFactors = FALSE)
rownames(membership) <- colnames(class.raw)



for(i in 1:length(colnames(class.raw))) {
  class.name <- membership[i, "class"]
	member.samples <- sample.names[which(class.raw[ ,class.name]==1)]
	#print(as.character(do.call(paste, c(as.list(member.samples), sep=" "))))
	membership[class.name,"members"] <- as.character(do.call(paste, c(as.list(member.samples), sep=",  ")))

}


# output the membership table
pdf("Figures/RPMM_cluster_raw.pdf", width=25, height=20)
par(mfrow = c(3, 1))
# plot the tree
plotTree.blcTree(mmfit.raw, )

plot.new()

plot.new()
library(gridExtra)
rownames(membership) <- NULL
grid.table(membership)
dev.off()



#######################################################################################
#####  Normalized and Filtered DATA
#####
#######################################################################################
##### Load normFilt data and the mmfit-fited result(modify the path accordingly)
load('/ubc/cs/research/irmtraud/people/jingyun/localFile_backup/Documents/stat540/project/Data/mmfit_normFilt.Rdata')
load('/ubc/cs/research/irmtraud/people/jingyun/localFile_backup/Documents/stat540/project/Data/All_3_sets_normAndFilt.Rdata')


# merge 3 data sets and design
merged.dat <- cbind(APL.dat, ALL.dat, CTRL.dat)
sample.names <- colnames(merged.dat)


# find how many groups resulted from the RPMM clustering 
class.normFilt <- blcTreeLeafMatrix(mmfit.normFilt)
dim(class.normFilt) # 9 groups for these 52 samples. 

# record the membership of the samples for the 9 groups
membership <- data.frame(class=colnames(class.normFilt),  members=vector(mode="character", length=length(colnames(class.normFilt))), stringsAsFactors = FALSE)
rownames(membership) <- colnames(class.normFilt)



for(i in 1:length(colnames(class.normFilt))) {
	class.name <- membership[i, "class"]
	member.samples <- sample.names[which(class.normFilt[ ,class.name]==1)]
	#print(as.character(do.call(paste, c(as.list(member.samples), sep=" "))))
	membership[class.name,"members"] <- as.character(do.call(paste, c(as.list(member.samples), sep=",  ")))

}


# output the membership table
pdf("Figures/RPMM_cluster_normFilt.pdf", width=25, height=20)
par(mfrow = c(3, 1))
# plot the tree
plotTree.blcTree(mmfit.normFilt, )

plot.new()

plot.new()
library(gridExtra)
rownames(membership) <- NULL
grid.table(membership)
dev.off()
     
