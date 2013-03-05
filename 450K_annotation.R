## 450K annotation package ##
source("http://bioconductor.org/biocLite.R")
biocLite('AnnotationForge')
biocLite("IlluminaHumanMethylation450k.db")

library(IlluminaHumanMethylation450k.db)



### 450k CPG COORDINATE  (probe location in genome)
coor <- IlluminaHumanMethylation450kCPGCOORDINATE
mapped_probes <- mappedkeys(coor)
coor<- as.data.frame(coor[mapped_probes])
head(coor)
nrow(coor)

# Probe chromosome (could be useful way to break data down)
CHR <- as.data.frame(IlluminaHumanMethylation450kCHR36)
head(CHR)
nrow(CHR)




#### HOW MANY ISLANDS? ####
###  27176 with and average of 11 probes in an island, 
#could serve as or regions to look for differential methylation
#but if we use this strategy we will lose 176112 probes which are not in islands
# have to look at these numbers again after initial QC steps



#UCSC CpG island (Range that Illumina has defined at CpG island)
Island<-as.data.frame(IlluminaHumanMethylation450kCPGINAME)
head(Island)
nrow(Island) #309465/485577 probes are in islands (including shelfs and shores)

isl.range<-Island$cpgiview.ucscname
# islands with multiple probes
multiple<-intersect(isl.range,isl.range)
length(multiple) # 27176 islands

## data frame of number of probes in each island
probes<-as.data.frame(table(Island[2]))
head(probes)
summary(probes) # Max 117 probes/island average 11 probes/island

# cpg island relation (what part of island is probe in?)
relation<-as.data.frame(IlluminaHumanMethylation450kCPGIRELATION)
head(relation)
nrow(relation)

# Aggregate the beta values of the probes for each CpG island.
cpginame <- as.data.frame(IlluminaHumanMethylation450kCPGINAME)
colnames(cpginame) <- c('Probe_ID', 'cpginame')
rownames(cpginame) <- cpginame$Probe_ID
cpginame$cpginame <- factor(cpginame$cpginame)
CTRL.cpginame <- merge(CTRL.dat, cpginame, by='row.names')
sampleNames <- colnames(CTRL.dat)
CTRL.cpgislands <- simplify2array(by(
	CTRL.cpginame[,sampleNames],
	list(CTRL.cpginame$cpginame),
	colMeans))
