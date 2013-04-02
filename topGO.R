#####################################################
##### # Geneset Enrichment Analysis
##### 1.) Input: Differentially Methylated Island List entitle "myInterestingIslands
##### 2.) Test for enriched GO groups in the list (Fisher and KS)
##### 3.) Convert Island name list to associated genename for input into pathway analysis
##### 4.) Call the genes in a GO group
#####################################################

setwd("/Users/Rachel/Desktop/UBC/540/Stats540/Project/R_objects")
library(IlluminaHumanMethylation450k.db)
library(topGO)

#### ALL ISLANDS 
#GO annotations for all Islands on the 450K
  GO <- as.list(IlluminaHumanMethylation450kGO2PROBE)
  head(GO) #GO<-GO2probe, want probe to GO then Island to GO
  probe2GO<-inverseList(GO)

#Summaize GO groups associated with each island (object for topGO function)
  Island<-as.data.frame(IlluminaHumanMethylation450kCPGINAME)
  #lookup all the GO of each Island and store as list
    islGO<-function(x) probe2GO[[Island[x,1]]]
    isl<-as.list(1:nrow(Island))
  #Island GO data (from probe GO data)
    Island2GO<-lapply(isl,islGO)
    names(Island2GO)<-Island$cpgiview.ucscname
    save(Island2GO, file="island2go.R")
  
    load(file="island2go.R")

### Produce topGO object
  ###list of all Islands
    islandNames<-Island$cpgiview.ucscname

  ### Top Differentially Methylated Islands from linear mixed-effects model
    lme_geneset_apl<-read.table(file="lme_ml.tab")
    top_lme<-subset(lme_geneset_apl, abs(t.value)>15) # t value cutoff
    myInterestingIslands<-as.character(top_lme$cgi)

    islandList <- factor(as.integer(islandNames %in% myInterestingIslands))
    names(islandList) <- islandNames
  
  ### topGO object
  GOdata <- new("topGOdata", ontology = "MF", allGenes = islandList, annot = annFUN.gene2GO, gene2GO = Island2GO)
      

### TESTS
  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  resultFisher
  resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
  resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")
  ###Summary of top GO groups
  allRes_APL <- GenTable(GOdata, classicFisher = resultFisher,
                     classicKS = resultKS, elimKS = resultKS.elim,orderBy = "elimKS", 
                     ranksOf = "classicFisher", topNodes = 10)
save(allRes_ALL, file="allRes_ALL.R") 

write.table(allRes_ALL, file="allRes_ALL.txt", sep="\t") 
write.table(allRes_APL, file="allRes_APL.txt", sep="\t") 
load(file="allRes_ALL.R")
################################################################

## Genes in top Islands
  x <- IlluminaHumanMethylation450kSYMBOL
  # Get the probe identifiers that are mapped to a gene symbol
  mapped_probes <- mappedkeys(x)
  xx <- as.data.frame(x[mapped_probes])
  gen.isl<-merge(Island, xx, by.x="cpgiview.Probe_ID", by.y="probe_id")
  # 219813 probes have Islands and genes
  # dont need probe Id anymore
  gen.isl[1]<-NULL
  gen.isl<-unique(gen.isl) #21263 Islands associated with 14770 genes

# function to pull out genes associated with top islands
int.genes<-gen.isl[gen.isl$cpgiview.ucscname %in% myInterestingIslands, 2]

lapply(int.genes, write, "intgenes.txt", append=TRUE)
# feed top genes into pathway analysis tool

################################################################
# Genes in top GOs

prbs<-GO[["GO:0046332"]]
gene<-subset(xx,probe_id %in%prbs)
gogene<-unique(gene$symbol)

lapply(gogene, write, "gogenes.txt", append=TRUE)
