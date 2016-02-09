codeDir <- get(load("codeDir.RData"))
source(file.path(codeDir, "intSiteLogic.R"))

sampleID <- 5
sampleID <- 81
message(sampleID)

completeMetadata <- get(load("completeMetadata.RData"))[sampleID,]
print(t(completeMetadata), quote=FALSE)  

workingDir=completeMetadata$alias
minPercentIdentity=completeMetadata$minPctIdent
maxAlignStart=completeMetadata$maxAlignStart
maxLength=2500
refGenome="hg18"
refGenome="mm9"

## processAlignments <- function(workingDir, minPercentIdentity, maxAlignStart, maxLength, refGenome)

codeDir <- get(load("codeDir.RData"))
source(file.path(codeDir, "programFlow.R"))#for get_reference_genome function

setwd(workingDir)
message("Entering ", workingDir)

load("keys.RData")
keys$uPID <- paste(keys$R2, keys$R1, sep=":")
keys <- dplyr::group_by(keys, uPID) %>% dplyr::mutate(count=n()) 
ukeys <- (dplyr::group_by(keys, uPID) %>%
          dplyr::mutate(count=n(), uid=1:n()) %>%
          dplyr::top_n(n=1, uid) )


processPsl <- function(algns, from, keys){
    
    stopifnot(from == "R1" | from == "R2")
    stopifnot(c("R2", "R1", "names", "uPID", "count") %in% colnames(keys))
    
    algns$from <- from
    algns <- merge(algns, keys, by.x="qName", by.y=from)
    
    algns <- (dplyr::mutate(algns,  POI=100*(matches+repMatches)/qSize) %>%
              dplyr::filter(strand=="+" | strand=="-") %>%
              dplyr::select(uPID, count, tName, strand, tStart, tEnd, from, POI, qStart))
    
    return(algns)
}

psl.R1 <- read.psl(list.files(".", "R1.*.fa.psl.gz"),
                   bestScoring=F, removeFile=F)
psl.R1 <- processPsl(psl.R1, from="R1", keys=ukeys)
psl.R1 <- dplyr::filter(psl.R1, POI>minPercentIdentity)

## complement R1 strand
strand <- c("+", "-")
psl.R1$Cstrand <- strand[3-match(psl.R1$strand, strand)]
psl.R1$strand <- psl.R1$Cstrand



psl.R2 <- read.psl(list.files(".", "R2.*.fa.psl.gz"),
                   bestScoring=F, removeFile=F)
psl.R2 <- processPsl(psl.R2, from="R2", keys=ukeys)
psl.R2 <- dplyr::filter(psl.R2, POI>minPercentIdentity &
                        qStart<=maxAlignStart)

psl.Pair <- merge(psl.R2,
                  psl.R1,
                  by=c("uPID", "tName", "strand"),
                  suffixes = c(".R2",".R1"))

rm(psl.R1, psl.R2)
gc()

psl.Pair <- dplyr::mutate(psl.Pair,
                          position=ifelse(strand=="+", tStart.R2, tEnd.R2),
                          breakpoint=ifelse(strand=="+", tEnd.R1, tStart.R1))

psl.Pair <- dplyr::filter(psl.Pair, abs(breakpoint-position)<maxLength)

psl.Pair <- dplyr::select(psl.Pair,
                          uPID,
                          chr=tName, strand, position, breakpoint,
                          count=count.R2)

psl.Pair <- (dplyr::group_by(psl.Pair, uPID) %>%
             dplyr::mutate(hits=n()))




#' clean up alignments and prepare for int site calling
#'
#' @param algns df with: 21 columns from BLAT (psl)
#' @param from is "R1" or "R2"  
#' @return Granges object

processBLATData <- function(algns, from){
stopifnot(from == "R1" | from == "R2")
algns$from <- from
algns <- merge(algns, keys[c(from, "names")], by.x="qName", by.y=from)
algns.gr <- GRanges(seqnames=Rle(algns$tName),
                    ranges=IRanges(start=algns$tStart, end=algns$tEnd),
                    strand=Rle(algns$strand),
                    seqinfo=seqinfo(get_reference_genome(refGenome)))

names(algns.gr) <- algns[,"names"]
mcols(algns.gr) <- algns[,c("matches", "repMatches", "misMatches", "qStart", "qEnd", "qSize", "tBaseInsert", "from")]
rm(algns)
algns.gr
}

load("keys.RData")

psl.R2 <- list.files(".", pattern="R2.*.fa.psl.gz")
message("R2 psl:\n", paste(psl.R2, collapse="\n"),"\n")
hits.R2 <- processBLATData(read.psl(psl.R2, bestScoring=F, removeFile=F),
"R2")
save(hits.R2, file="hits.R2.RData")

psl.R1 <- list.files(".", pattern="R1.*.fa.psl.gz")
message("R1 psl:\n", paste(psl.R1, collapse="\n"),"\n")
hits.R1 <- processBLATData(read.psl(psl.R1, bestScoring=F, removeFile=F),
"R1")
save(hits.R1, file="hits.R1.RData")

load("stats.RData")

                                        #no more '.p' or '.l' nomenclature here
  #we're combining alignments from both sides of the read

  
  allAlignments <- append(hits.R1, hits.R2)
  stopifnot(!any(strand(allAlignments)=="*"))
  #TODO: star strand is impossible '*' 
  
  readsAligning <- length(unique(names(allAlignments)))
  
  stats <- cbind(stats, readsAligning)
  cat("readsAligning:\t", readsAligning, "\n")
  
  allAlignments$percIdent <- 100*(allAlignments$matches +
                                  allAlignments$repMatches)/allAlignments$qSize
  
  allAlignments <- subset(allAlignments,
                          allAlignments$percIdent >= minPercentIdentity &
                          allAlignments$qStart <= maxAlignStart &
                          allAlignments$tBaseInsert <= 5)
  
  #even if a single block spans the vast majority of the qSize, it's NOT ok to
  #accept the alignment as it will give a spurrious integration site/breakpoint
  #that's a few dozen nt away from the real answer.  That won't be caught by our
  #collapsing algorithms
  
  readsWithGoodAlgnmts <- length(unique(names(allAlignments)))
  
  stats <- cbind(stats, readsWithGoodAlgnmts)
  cat("readsWithGoodAlgnmts:\t", readsWithGoodAlgnmts, "\n")
  
  # allAlignments now are GRangesList
  # and later only keep concordant pairs
  allAlignments <- split(allAlignments, names(allAlignments))
 
  # need alignemnent when R1 is on one strand and R2 is on the opposite 
  numStrands <- unname(table(strand(allAlignments)))
  
  #just a quick pre-filter to reduce amount of work in future steps
  allAlignments <- subset(allAlignments, numStrands[,1]>=1 & numStrands[,2]>=1)
  
  ######## REDUCE ALIGNMENTS INTO POTENTIAL SITES AT THE READ-LEVEL ###########  
  #private method stuff is so that we can maintain the revmap
  #doing sapply(split(allAlignments, names(allAlignments)), reduce, min.gapwidth=maxLength) is incredibly slow
  #might be faster with dplyr or data.table?
  pairedAlignments <- GenomicRanges:::deconstructGRLintoGR(flank(allAlignments, -1, start=T))
  pairedAlignments <- reduce(pairedAlignments, min.gapwidth=maxLength, with.revmap=TRUE, ignore.strand=TRUE)
  pairedAlignments <- GenomicRanges:::reconstructGRLfromGR(pairedAlignments, flank(allAlignments, -1, start=T))
  pairedAlignments <- unlist(pairedAlignments)
  # strand can be '*" here because paired has + and -

  #names are no longer unique identifier
  #candidate sites are either a unique site, a member of a multihit cluster, or a member of a chimera
  pairedAlignments$pairingID <- seq(pairedAlignments)

  
  ########## IDENTIFY PROPERLY-PAIRED READS ##########
  #properly-paired reads will have one representitive each from R1 and R2
  #there can be multiple candidate sites per read
  alignmentsPerPairing <- sapply(pairedAlignments$revmap, length)
  
  allAlignments <- unlist(allAlignments, use.names=FALSE)
  alignmentSources <- split(allAlignments$from[unlist(pairedAlignments$revmap)],
                            as.vector(Rle(pairedAlignments$pairingID, alignmentsPerPairing)))
  # alignmentSources have values of only "R1" and "R2" 
  
  
  R1Counts <- sapply(alignmentSources, function(x){sum(x=="R1")})
  R2Counts <- sapply(alignmentSources, function(x){sum(x=="R2")})
  oneEach <- R1Counts==1 & R2Counts==1
  
  sitesFrom2Alignments <- pairedAlignments[alignmentsPerPairing==2]
  properlyPairedAlignments <- sitesFrom2Alignments[oneEach[sitesFrom2Alignments$pairingID]]
  #still can be with multihits
  
  #assign strand to be whatever was seen on the LTR read (i.e. R2) in allAlignments
  allPairedSingleAlignments <- allAlignments[unlist(properlyPairedAlignments$revmap)]
  #guaranteed to have R1 and R2 at this point
  R1s <- allPairedSingleAlignments[allPairedSingleAlignments$from=="R1"]
  R2s <- allPairedSingleAlignments[allPairedSingleAlignments$from=="R2"]
  strand(properlyPairedAlignments) <- strand(R2s)

  #need to kick out properlyPaired Alignments that are 'outward facing'
  #these will have union widths that are far greater than the reduction widths
  unionWidths <- width(punion(R1s, R2s, ignore.strand=T, fill.gap=T))

  properlyPairedAlignments <- properlyPairedAlignments[abs(unionWidths-width(properlyPairedAlignments)) < 5]

  numProperlyPairedAlignments <- length(unique(names(properlyPairedAlignments)))

  stats <- cbind(stats, numProperlyPairedAlignments)
  cat("numProperlyPairedAlignments:\t", numProperlyPairedAlignments, "\n")
  
  ########## IDENTIFY MULTIPLY-PAIRED READS (multihits) ##########  
  properlyPairedAlignments$sampleName <- sapply(strsplit(names(properlyPairedAlignments), "%"), "[[", 1)
  properlyPairedAlignments$ID <- sapply(strsplit(names(properlyPairedAlignments), "%"), "[[", 2)

  multihitNames <- unique(names(properlyPairedAlignments[duplicated(properlyPairedAlignments$ID)]))
  unclusteredMultihits <- subset(properlyPairedAlignments, names(properlyPairedAlignments) %in% multihitNames)
  ##unclusteredMultihits <- standardizeSites(unclusteredMultihits) #not sure if this is required anymore

  ########## IDENTIFY UNIQUELY-PAIRED READS (real sites) ##########  
  allSites <- properlyPairedAlignments[!properlyPairedAlignments$ID %in% unclusteredMultihits$ID]
  
  save(allSites, file="rawSites.RData")
  
  ## standardizing sites will be done in reportMaker
  ##allSites <- standardizeSites(allSites)
  sites.final <- dereplicateSites(allSites)
  
  if(length(sites.final)>0){
    sites.final$sampleName <- allSites[1]$sampleName
    sites.final$posid <- paste0(as.character(seqnames(sites.final)),
                                as.character(strand(sites.final)),
                                start(flank(sites.final, width=-1, start=TRUE)))
    }
  
  save(sites.final, file="sites.final.RData")
  save(allSites, file="allSites.RData")

  numAllSingleReads <- length(allSites)
  stats <- cbind(stats, numAllSingleReads)
  cat("numAllSingleReads:\t", numAllSingleReads, "\n")
  numAllSingleSonicLengths <- 0
  if( length(sites.final)>0 ) {
        numAllSingleSonicLengths <- length(unlist(sapply(1:length(sites.final), function(i){
        unique(width(allSites[sites.final$revmap[[i]]]))})))
  }
  stats <- cbind(stats, numAllSingleSonicLengths)
  cat("numAllSingleSonicLengths:\t", numAllSingleSonicLengths, "\n")
  numUniqueSites <- length(sites.final)
  stats <- cbind(stats, numUniqueSites)
  cat("numUniqueSites:\t", numUniqueSites, "\n")

########## IDENTIFY IMPROPERLY-PAIRED READS (chimeras) ##########
  singletonAlignments <- pairedAlignments[alignmentsPerPairing==1]
  strand(singletonAlignments) <- strand(allAlignments[unlist(singletonAlignments$revmap)])
  t <- table(names(singletonAlignments))
  chimeras <- subset(singletonAlignments, names(singletonAlignments) %in% 
                       names(subset(t, t==2))) #should be >=?
  #not an already-assigned read
  chimeras <- chimeras[!names(chimeras) %in% names(properlyPairedAlignments)]
  chimeras <- split(chimeras, names(chimeras))
  
  chimeras <- subset(chimeras,
                     sapply(chimeras, function(x){sum(R1Counts[x$pairingID]) ==
                                                    sum(R2Counts[x$pairingID])}))
  
  dereplicatedChimeras <- dereplicateSites(unlist(chimeras, use.names=FALSE))
  
  chimera <- length(dereplicatedChimeras)
  
  stats <- cbind(stats, chimera)
  cat("chimera:\t", chimera, "\n")
  chimeraData <- list("chimeras"=chimeras, "dereplicatedChimeras"=dereplicatedChimeras)
  save(chimeraData, file="chimeraData.RData")
  save(stats, file="stats.RData")
  

  ########## IDENTIFY MULTIPLY-PAIRED READS (multihits) ##########  
  clusteredMultihitPositions <- GRangesList()
  clusteredMultihitLengths <- list()
  
  if(length(unclusteredMultihits) > 0){
      ##library("igraph") #taken care of earlier on in the code
      
      ##medians are based on all the potential sites for a given read, so we want these counts preserved pre-condensentaiton
    multihits.medians <- round(median(width(split(unclusteredMultihits, unclusteredMultihits$ID))))

    #condense identical R1/R2 pairs (as determined by keys.RData) and add the new 
    #primary key (the combination of R1's key ID + R2's key ID -> readPairKey)
    keys$names <- sapply(strsplit(as.character(keys$names), "%"), "[[", 2)
    keys$readPairKey <- paste0(keys$R1, "_", keys$R2)
    mcols(unclusteredMultihits)$readPairKey <- keys[match(mcols(unclusteredMultihits)$ID, keys$names), "readPairKey"] #merge takes too much memory
    
    multihits.split <- unique(split(unclusteredMultihits, unclusteredMultihits$readPairKey))
    multihits.split <- flank(multihits.split, -1, start=T) #now just care about solostart
    
    overlaps <- findOverlaps(multihits.split, multihits.split, maxgap=5)
    edgelist <- matrix(c(queryHits(overlaps), subjectHits(overlaps)), ncol=2)
    
    clusteredMultihitData <- clusters(graph.edgelist(edgelist, directed=F))
    clusteredMultihitNames <- split(names(multihits.split), clusteredMultihitData$membership)
    clusteredMultihitPositions <- GRangesList(lapply(clusteredMultihitNames, function(x){
      unname(granges(unique(unlist(multihits.split[x]))))
    }))
    clusteredMultihitLengths <- lapply(clusteredMultihitNames, function(x){
      #retrieve pre-condensed read IDs, then query for median fragment length
      readIDs <- unique(unclusteredMultihits[unclusteredMultihits$readPairKey %in% x]$ID)
      data.frame(table(multihits.medians[readIDs]))
    })
  }
  stopifnot(length(clusteredMultihitPositions)==length(clusteredMultihitLengths))
  multihitData <- list(unclusteredMultihits, clusteredMultihitPositions, clusteredMultihitLengths)
  names(multihitData) <- c("unclusteredMultihits", "clusteredMultihitPositions", "clusteredMultihitLengths")

  save(multihitData, file="multihitData.RData")
  
  #making new variable multihitReads so that the naming in stats is nice
  ##multihitReads <- length(multihitNames) #multihit names is already unique
  multihitReads <- length(unique(multihitData$unclusteredMultihits$ID))
  stats <- cbind(stats, multihitReads)
  cat("multihitReads:\t", multihitReads, "\n")
  multihitSonicLengths <- 0
  if( length(multihitData$clusteredMultihitLengths)>0 ) {
        multihitSonicLengths <- sum(sapply(multihitData$clusteredMultihitLengths, nrow))
  }
  stats <- cbind(stats, multihitSonicLengths) 
  cat("multihitSonicLengths:\t", multihitSonicLengths, "\n")
  multihitClusters <- length(multihitData$clusteredMultihitPositions) #
  stats <- cbind(stats, multihitClusters)
  cat("multihitClusters:\t", multihitClusters, "\n")
  
  totalSonicLengths <- numAllSingleSonicLengths + multihitSonicLengths
  stats <- cbind(stats, totalSonicLengths)
  cat("totalSonicLengths:\t", totalSonicLengths, "\n")
  totalEvents <- numUniqueSites + multihitClusters
  stats <- cbind(stats, totalEvents)
  cat("totalEvents:\t", totalEvents, "\n")
  save(stats, file="stats.RData")
  
}
