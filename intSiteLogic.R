## load hiReadsProcessor.R
libs <- c("plyr", "BiocParallel", "Biostrings", "GenomicAlignments" ,"hiAnnotator" ,"sonicLength", "GenomicRanges", "BiocGenerics")
junk <- sapply(libs, require, character.only=TRUE)
if( any(!junk) ) {
    message("Libs not loaded:")
    print(data.frame(Loaded=junk[!junk]))
    stop()
}
codeDir <- get(load("codeDir.RData"))
stopifnot(file.exists(file.path(codeDir, "hiReadsProcessor.R")))
source(file.path(codeDir, "hiReadsProcessor.R"))
source(file.path(codeDir, "standardization_based_on_clustering.R"))

## nesessary libraries
stopifnot(require("ShortRead"))
stopifnot(require("GenomicRanges"))
stopifnot(require("igraph"))


#' find reads originating from vector
#' @param vectorSeq vector sequence fasta file
#' @param primerLTR primer and LTR sequence
#' @param reads.p DNAStringSet, reads on primer side
#' @param reads.l DNAStringSet, reads on linker side
#' @return character, qNames for the vector reads
#' @example findVectorReads(vectorSeq, reads.p, reads.l)
findVectorReads <- function(vectorSeq, primerLTR="GAAAATCTCTAGCA",
                            reads.p, reads.l,
                            debug=FALSE) {
    ##require(stringr)
    require(Biostrings)
    
    Vector <- readDNAStringSet(vectorSeq)
    
    ltrLoci <- stringr::str_locate_all(Vector, primerLTR)[[1]][, "start"]
    if( length(ltrLoci)>2 ) stop("Expecting 2 LTR regions, got\n",
                                 paste(ltrLoci, collapse=", "))
    if( length(ltrLoci)<1 ) stop("Expecting 2 LTR regions, got\n",
                                 paste(ltrLoci, collapse=", "))
    ltrpos <- ltrLoci[1]
    
    globalIdentity <- 0.75
    blatParameters <- c(minIdentity=70, minScore=15, stepSize=3, 
                        tileSize=8, repMatch=112312, dots=1000, 
                        q="dna", t="dna", out="psl")
    
    
    hits.v.p <- try(read.psl(blatSeqs(query=reads.p, subject=Vector,     
                                      blatParameters=blatParameters, parallel=F),
                             bestScoring=F) )
    if( class(hits.v.p) == "try-error" ) hits.v.p <- data.frame()
    if ( debug ) save(hits.v.p, file="hits.v.p.RData")    
    
    hits.v.l <- try(read.psl(blatSeqs(query=reads.l, subject=Vector, 
                                      blatParameters=blatParameters, parallel=F),
                             bestScoring=F) )
    if( class(hits.v.l) == "try-error" ) hits.v.l <- data.frame()
    if ( debug ) save(hits.v.l, file="hits.v.l.RData")    
    
    hits.v.p <- dplyr::filter(hits.v.p, tStart  > ltrpos &
                                        ##tStart  < ltrpos+nchar(primerLTR)+10 &
                                        matches > globalIdentity*qSize &
                                        ##strand == "+" &
                                        qStart  <= 5) 
    hits.v.l <- dplyr::filter(hits.v.l, matches>globalIdentity*qSize )
                                        ##strand=="-") 
    hits.v <- try(merge(hits.v.p[, c("qName", "tStart")],
                        hits.v.l[, c("qName", "tStart")],
                        by="qName")
                 ,silent = TRUE)
    if( class(hits.v) == "try-error" ) hits.v <- data.frame()
    
    ##hits.v <- dplyr::filter(hits.v, tStart.y >= tStart.x &
    ##                                tStart.y <= tStart.x+2000)
    
    if ( debug ) {
        save(reads.p, file="reads.p.RData")
        save(reads.l, file="reads.l.RData")
    }
    
    vqName <- unique(hits.v$qName)
    
    message(length(vqName), " vector sequences found")
    return(vqName)
}
## vqName <- findVectorReads(vectorSeq, reads.p, reads.l)

getTrimmedSeqs <- function(qualityThreshold, badQuality, qualityWindow, primer,
                           ltrbit, largeLTRFrag, linker, linker_common, mingDNA,
                           read1, read2, alias, vectorSeq){
  
  ##### Load libraries #####
  ##library("hiReadsProcessor")
  ##library("ShortRead")
  
  stats <- data.frame()
  message(alias)
  workingDir <- alias
  suppressWarnings(dir.create(workingDir, recursive=TRUE))
  setwd(workingDir)
  
  stats.bore <- data.frame(sample=alias)
  message("\t read data and make quality heatmaps")  
  
  filenames <- list(read1, read2)
  
  reads <- lapply(filenames, function(x) {
    sapply(x, readFastq)
  })
  
  stats.bore$Reads.l.beforeTrim <- sum(sapply(reads[[1]], length))
  stats.bore$Reads.p.beforeTrim <- sum(sapply(reads[[2]], length))
  
  r <- lapply(reads, function(x){
    seqs <- x[[1]]
    if(length(seqs) > 0){
      #remove anything after 5 bases under Q30 in 10bp window
      trimmed <- trimTailw(seqs, badQuality, qualityThreshold,
                           round(qualityWindow/2))
      #get rid of anything that lost more than half the bases
      trimmed <- trimmed[width(trimmed) > 65]
      if(length(trimmed) > 0){
        trimmedSeqs <- sread(trimmed)
        trimmedqSeqs <- quality(quality(trimmed))
        names(trimmedSeqs) <- names(trimmedqSeqs) <- 
          sapply(sub("(.+) .+","\\1",ShortRead::id(trimmed)),
                 function(z){paste0(alias, "%", strsplit(z, "-")[[1]][2])})
      }
    }
    list(trimmedSeqs, trimmedqSeqs)
  })
  
  reads <- sapply(r, "[[", 1)
  qualities <- sapply(r, "[[", 2)
  #this is needed for primerID quality scores later on
  R1Quality <- qualities[[1]]
  
  stats.bore$Reads.p.afterTrim <- length(reads[[2]])
  stats.bore$Reads.l.afterTrim <- length(reads[[1]])
  print(stats.bore) 
  
  message("\t trim adaptors")
  
  #'.p' suffix signifies the 'primer' side of the amplicon (i.e. read2)
  #'.l' suffix indicates the 'liner' side of the amplicon (i.e. read1)
  
  res.p <- pairwiseAlignSeqs(reads[[2]], patternSeq=primer,
                             qualityThreshold=1, doRC=F)
  
  reads.p <- reads[[2]]
  if(length(res.p) > 0){
    reads.p <- trimSeqs(reads[[2]], res.p, side='left', offBy=1)
  }
  
  stats.bore$primed <- length(reads.p)
  
  res.ltr <- pairwiseAlignSeqs(reads.p, patternSeq=ltrbit, 
                               qualityThreshold=1, doRC=F)
  
  if(length(res.ltr) > 0 ){
    reads.p <- trimSeqs(reads.p, res.ltr, side='left', offBy=1)
  }
  
  stats.bore$LTRed <- length(reads.p)
  
  if(grepl("N", linker)){
    res <- primerIDAlignSeqs(subjectSeqs=reads[[1]], patternSeq=linker,
                             doAnchored=T, qualityThreshold1=1, 
                             qualityThreshold2=1, doRC=F)
    res.l <- res[["hits"]]
    res.pID <- res[["primerIDs"]]
  }else{
    res.l <- pairwiseAlignSeqs(subjectSeqs=reads[[1]], patternSeq=linker,
                               qualityThreshold=.95, doRC=F, side="middle")
    start(res.l) <- 1
  }
  
  reads.l <- reads[[1]]
  if(length(res.l) > 0 ){
    reads.l <- trimSeqs(reads[[1]], res.l, side='left', offBy=1)
    if(grepl("N", linker)){ #i.e. contains a primerID
      R1Quality <- R1Quality[match(names(res.pID), names(R1Quality))]
      primerIDs <- trimSeqs(reads[[1]], res.pID, side="middle")
      primerIDQuality <- subseq(R1Quality, start=start(res.pID),
                                end=end(res.pID))
      primerIDData <- list(primerIDs, primerIDQuality)
      
      save(primerIDData, file="primerIDData.RData")
    }
  }
  
  stats.bore$linkered <- length(reads.l)
  
  print(stats.bore) 
  
  ## check if reads were sequenced all the way by checking for opposite adaptor
  message("\t trim opposite side adaptors")
  
  res.p <- NULL
  tryCatch(res.p <- pairwiseAlignSeqs(reads.p, linker_common,
                                      qualityThreshold=.55, side='middle',
                                      doRC=F),
           error=function(e){print(paste0("Caught ERROR in intSiteLogic: ",
                                          e$message))})
  
  if(!is.null(res.p)){
    #if we see the common sequence, pitch the rest
    end(res.p) <- width(reads.p[names(res.p)]) + 1
    if(length(res.p) > 0 ){
      reads.p <- c(reads.p[!names(reads.p) %in% names(res.p)],
                   trimSeqs(reads.p, res.p, side='right', offBy=1))
    }
  }
  
  res.l <- NULL
  tryCatch(res.l <- pairwiseAlignSeqs(reads.l, largeLTRFrag,
                                      qualityThreshold=.55, side='middle', doRC=F),
           error=function(e){print(paste0("Caught ERROR in intSiteLogic: ",
                                          e$message))})
  
  if(!is.null(res.l)){
    end(res.l) <- width(reads.l[names(res.l)]) + 1
    if(length(res.l) > 0 ){
      reads.l <- c(reads.l[!names(reads.l) %in% names(res.l)],
                   trimSeqs(reads.l, res.l, side='right', offBy=1))
    }
  }
  
  
  reads.p <- subset(reads.p, width(reads.p) > mingDNA)
  reads.l <- subset(reads.l, width(reads.l) > mingDNA)
  
  stats.bore$reads.p_afterTrim <- length(reads.p)
  stats.bore$reads.l_afterTrim <- length(reads.l)
  
  
  message("\t trim vector") 
  vqName <- findVectorReads(file.path("..", vectorSeq),
                            paste0(primer, ltrbit),
                            reads.p, reads.l,
                            debug=TRUE)
  
  reads.p <- reads.p[!names(reads.p) %in% vqName]
  reads.l <- reads.l[!names(reads.l) %in% vqName]
  
  
  stats.bore$reads.p_afterVTrim <- length(reads.p)
  stats.bore$reads.l_afterVTrim <- length(reads.l)
  
  toload <- intersect(names(reads.p), names(reads.l))
  
  stats.bore$reads.lLength <- as.integer(mean(width(reads.l)))  
  stats.bore$reads.pLength <- as.integer(mean(width(reads.p)))
  
  stats.bore$curated <- length(toload)
  
  print(stats.bore)
  stats <- rbind(stats, stats.bore)
  
  reads.p <- reads.p[toload]
  reads.l <- reads.l[toload]
  
  #dereplicate seqs for faster alignments
  #this is re-expand at the beginning of callSeqs
  reads.p.u <- unique(reads.p)
  reads.l.u <- unique(reads.l)
  
  names(reads.p.u) <- seq(reads.p.u)
  names(reads.l.u) <- seq(reads.l.u)
  
  keys <- data.frame("R2"=match(reads.p, reads.p.u),
                     "R1"=match(reads.l, reads.l.u), "names"=toload)
  
  save(keys, file="keys.RData")
  
  if(length(toload) > 0){
    #cap number of reads per thread-we care about speed rather than # of procs
    chunks.p <- split(seq_along(reads.p.u), ceiling(seq_along(reads.p.u)/30000))
    for(i in c(1:length(chunks.p))){
      writeXStringSet(reads.p.u[chunks.p[[i]]], file=paste0("R2-", i, ".fa"),
                      append=TRUE)
    }
    
    chunks.l <- split(seq_along(reads.l.u), ceiling(seq_along(reads.l.u)/30000))
    for(i in c(1:length(chunks.l))){    
      writeXStringSet(reads.l.u[chunks.l[[i]]], file=paste0("R1-", i, ".fa"),
                      append=TRUE)
    }
    
    save(stats, file="stats.RData")
    alias #return 'value' which ultimately gets saved as trimStatus.RData
  }else{
    stop("error - no curated reads")
  }
}

processAlignments <- function(workingDir, minPercentIdentity, maxAlignStart, maxLength, refGenome){
  
  codeDir <- get(load("codeDir.RData"))
  source(paste0(codeDir, "/programFlow.R"))#for get_reference_genome function
  
  setwd(workingDir)
  cat(workingDir,"\n")

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
    mcols(algns.gr) <- algns[,c("matches", "qStart", "qSize", "tBaseInsert", "from")]
    rm(algns)
    algns.gr
  }
  
  load("keys.RData")
  
  hits.R2 <- processBLATData(read.psl(system("ls R2*.fa.psl.gz", intern=T), bestScoring=F, removeFile=F), "R2")
  save(hits.R2, file="hits.R2.RData")
  
  hits.R1 <- processBLATData(read.psl(system("ls R1*.fa.psl.gz", intern=T), bestScoring=F, removeFile=F), "R1")
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
  
  allAlignments$percIdent <- 100 * allAlignments$matches/allAlignments$qSize
  
  #doing this first subset speeds up the next steps
  allAlignments <- subset(allAlignments,
                          allAlignments$percIdent >= minPercentIdentity
                          & allAlignments$qStart <= maxAlignStart
                          & allAlignments$tBaseInsert <= 5)
  
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
  unclusteredMultihits <- standardizeSites(unclusteredMultihits) #not sure if this is required anymore

  ########## IDENTIFY UNIQUELY-PAIRED READS (real sites) ##########  
  allSites <- properlyPairedAlignments[!properlyPairedAlignments$ID %in% unclusteredMultihits$ID]
  
  save(allSites, file="rawSites.RData")
  
  allSites <- standardizeSites(allSites)
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

    #medians are based on all the potential sites for a given read, so we want these counts preserved pre-condensentaiton
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
