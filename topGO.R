#####################################################
##### # Geneset Enrichment Analysis
##### 1.) Input: Differentially Methylated Island List entitle "myInterestingIslands
##### 2.) Test for enriched GO groups in the list (Fisher and KS)
##### 3.) Convert Island name list to associated genename for input into pathway analysis
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
  ###Diffrentially Methylated Island List (random sample of islands now)
  ###
    # remove this random sampling line when have actual list (change to loading list)
    myInterestingIslands <- sample(names(Island2GO), length(names(Island2GO)) / 1000) 
  ###
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
  allRes <- GenTable(GOdata, classicFisher = resultFisher,
                     classicKS = resultKS, elimKS = resultKS.elim,orderBy = "elimKS", 
                     ranksOf = "classicFisher", topNodes = 10)
          

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


