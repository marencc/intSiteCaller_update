getTrimmedSeqs <- function(qualityThreshold, badQuality, qualityWindow, primer,
                          ltrbit, largeLTRFrag, linker, linker_common, mingDNA,
                          read1, read2, alias, vectorSeq){
  
  ##### Load libraries #####
  libs <- c("hiAnnotator", "hiReadsProcessor", "plyr", "ShortRead",
            "BiocParallel")
  sapply(libs, require, character.only=TRUE)
  
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
           error=function(e){print(paste0("Caught ERROR in pairedIntSites: ",
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
           error=function(e){print(paste0("Caught ERROR in pairedIntSites: ",
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
  #we want to do this at the end so that we don't have to worry about partial
  #vector alignments secondary to incomplete trimming of long reads
  
  Vector <- readDNAStringSet(vectorSeq)
  
  blatParameters <- c(minIdentity=70, minScore=5, stepSize=3, 
                      tileSize=8, repMatch=112312, dots=1000, 
                      q="dna", t="dna", out="psl")
  
  findAndRemoveVector <- function(reads, Vector, blatParameters, minLength=10){
    
    hits.v <- read.psl(blatSeqs(query=reads, subject=Vector, 
                                blatParameters=blatParameters, parallel=F),
                       bestScoring=F)
    
    #collapse instances where a single read has multiple vector alignments
    hits.v <- reduce(GRanges(seqnames=hits.v$qName, IRanges(hits.v$qStart,
                                                            hits.v$qEnd)),
                     min.gapwidth=1200)
    names(hits.v) <- as.character(seqnames(hits.v))
    
    hits.v <- hits.v[start(hits.v)<=5 & width(hits.v)>minLength]
    
    reads[!names(reads) %in% names(hits.v)]
    
  }
  
  tryCatch(reads.p <- findAndRemoveVector(reads.p, Vector,
                                          blatParameters=blatParameters),
           error=function(e){print(paste0("Caught ERROR in pairedIntSites: ",
                                          e$message))})
  
  tryCatch(reads.l <- findAndRemoveVector(reads.l, Vector,
                                          blatParameters=blatParameters),
           error=function(e){print(paste0("Caught ERROR in pairedIntSites: ",
                                          e$message))})
  
  stats.bore$reads.p_afterVTrim <- length(reads.p)
  stats.bore$reads.l_afterVTrim <- length(reads.l)
  
  toload <- intersect(names(reads.p), names(reads.l))
  
  stats.bore$reads.lLength <- mean(width(reads.l))  
  stats.bore$reads.pLength <- mean(width(reads.p))
  
  stats.bore$curated <- length(toload)
  
  print(stats.bore)
  stats <- rbind(stats, stats.bore)
  
  #this could probably be cleaner with sapplys
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
    chunks.p <- split(seq_along(reads.p.u), ceiling(seq_along(reads.p.u)/20000))
    for(i in c(1:length(chunks.p))){
      writeXStringSet(reads.p.u[chunks.p[[i]]], file=paste0("R2-", i, ".fa"),
                      append=TRUE)
    }
    
    chunks.l <- split(seq_along(reads.l.u), ceiling(seq_along(reads.l.u)/20000))
    for(i in c(1:length(chunks.l))){    
      writeXStringSet(reads.l.u[chunks.l[[i]]], file=paste0("R1-", i, ".fa"),
                      append=TRUE)
    }
    
    save(stats, file="stats.RData")
    getwd()
  }else{
    stop("error - no curated reads")
  }
}

processAlignments <- function(workingDir, minPercentIdentity, maxAlignStart, maxLength, refGenome){
  
  ##### Load libraries #####
  libs <- c("hiAnnotator", "hiReadsProcessor", "plyr", "GenomicRanges", "sonicLength")
  sapply(libs, require, character.only=TRUE)
  
  setwd(workingDir)
  
  dereplicateSites <- function(uniqueSites){
    if(length(uniqueSites)>0){
      sites.reduced <- reduce(flank(uniqueSites, 5, both=TRUE), min.gapwidth=5, with.revmap=T)
      sites.reduced$counts <- sapply(sites.reduced$revmap, length)
      
      allSites <- uniqueSites[unlist(sites.reduced$revmap)]
      allSites <- split(allSites, Rle(values=c(1:length(sites.reduced)), lengths=sites.reduced$counts))
      allSites <- unlist(reduce(allSites, min.gapwidth=5))
      mcols(allSites) <- mcols(sites.reduced)
      allSites
    }else{
      uniqueSites
    }
  }
  
  get_reference_genome <- function(reference_genome) {
    pattern <- paste0("\\.", reference_genome, "$")
    match_index <- which(grepl(pattern, installed.genomes()))
    if (length(match_index) != 1) {
      write("Installed genomes are:", stderr())
      write(installed.genomes(), stderr())
      stop(paste("Cannot find unique genome for", reference_genome))
    }
    BS_genome_full_name <- installed.genomes()[match_index]
    library(BS_genome_full_name, character.only=T)
    get(BS_genome_full_name)
  }

  # clean up alignments and prepare for int site calling
  processBLATData <- function(algns, from){
    algns$from <- from
    algns <- merge(algns, keys[c(from, "names")], by.x="qName", by.y=from)
    algns.gr <- GRanges(seqnames=Rle(algns$tName),
                       ranges=IRanges(start=algns$tStart, end=algns$tEnd),
                       strand=Rle(algns$strand),
                       seqinfo=seqinfo(get_reference_genome(refGenome)))
    
    names(algns.gr) <- algns[,"names"]
    mcols(algns.gr) <- algns[,c("matches", "qStart", "qSize", "tBaseInsert", "blockSizes", "from")]
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
  
  readsAligning <- length(unique(names(allAlignments)))
  
  stats <- cbind(stats, readsAligning)
  
  allAlignments$percIdent <- 100 * allAlignments$matches/allAlignments$qSize

  #doing this first subset speeds up the next steps
  allAlignments <- subset(allAlignments, allAlignments$percIdent >= minPercentIdentity
                          & allAlignments$qStart <= maxAlignStart)

  #if a single block spans the vast majority of the qSize, then it's ok if there's a
  #big gap before a few straggling basepairs that were mismapped
  minGoodBlockSize <- floor(allAlignments$qSize*(minPercentIdentity/100))
  blockSizes <- sapply(strsplit(allAlignments$blockSizes, ","), as.integer)
  blockSizeOverride <- sapply(seq(length(blockSizes)), function(x){
    any(blockSizes[[x]] >= minGoodBlockSize[x])
  })

  allAlignments = allAlignments[(allAlignments$tBaseInsert <= 5) | blockSizeOverride]
  
  readsWithGoodAlgnmts <- length(unique(names(allAlignments)))
  
  stats <- cbind(stats, readsWithGoodAlgnmts)
  
  #turn allAlignments into soloStarts - makes reductions easier and more robust
  allAlignments <- flank(allAlignments, -1)
  
  allAlignments <- split(allAlignments, names(allAlignments))
  
  numStrands <- unname(table(strand(allAlignments)))
  
  #just a quick pre-filter to reduce amount of work in future steps
  allAlignments <- subset(allAlignments, numStrands[,1]>=1 & numStrands[,2]>=1)
  
  
  ######## REDUCE ALIGNMENTS INTO POTENTIAL SITES AT THE READ-LEVEL ###########
  #private method stuff is so that we can quickly do a read-by-read reduction
  #doing sapply(split(allAlignments, names(allAlignments)), reduce, min.gapwidth=maxLength) is incredibly slow
  #might be faster with dplyr or data.table?
  pairedAlignments <- GenomicRanges:::deconstructGRLintoGR(allAlignments)
  #reduce() doesn't like different strands - set as "*" now and add strand info back in later
  strand(pairedAlignments) <- "*"
  pairedAlignments <- reduce(pairedAlignments, min.gapwidth=maxLength, with.revmap=TRUE)
  pairedAlignments <- GenomicRanges:::reconstructGRLfromGR(pairedAlignments, allAlignments)
  
  pairedAlignments <- unlist(pairedAlignments)
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
  
  R1Counts <- sapply(alignmentSources, function(x){sum(x=="R1")})
  R2Counts <- sapply(alignmentSources, function(x){sum(x=="R2")})
  oneEach <- R1Counts==1 & R2Counts==1
  
  sitesFrom2Alignments <- pairedAlignments[alignmentsPerPairing==2]
  properlyPairedAlignments <- sitesFrom2Alignments[oneEach[sitesFrom2Alignments$pairingID]]
  
  #assign strand to be whatever was seen on the LTR read (i.e. R2) in allAlignments
  strandDonor <- allAlignments[unlist(properlyPairedAlignments$revmap)]
  strandDonor <- strandDonor[strandDonor$from=="R2"]
  strand(properlyPairedAlignments) <- strand(strandDonor)
  
  numProperlyPairedAlignments <- length(unique(names(properlyPairedAlignments)))
  
  stats <- cbind(stats, numProperlyPairedAlignments)
  
  
  ########## IDENTIFY MULTIPLY-PAIRED READS (multihits) ##########  
  properlyPairedAlignments$clone <- sapply(strsplit(names(properlyPairedAlignments), "%"), "[[", 1)
  properlyPairedAlignments$ID <- sapply(strsplit(names(properlyPairedAlignments), "%"), "[[", 2)
  
  multihitNames <- unique(names(properlyPairedAlignments[duplicated(properlyPairedAlignments$ID)]))
  multihits <- subset(properlyPairedAlignments, names(properlyPairedAlignments) %in% multihitNames)
  dereplicatedMultihits <- dereplicateSites(multihits)
  
  #####
  #CLUSTER MULTIHITS HERE! (save as clusteredMultihits)
  #####
  clusteredMultihits<-"clusteredMultihits"
  
  multihitData <- list(multihits, dereplicatedMultihits, clusteredMultihits)
  
  save(multihitData, file="multihitData.RData")
  
  #making new variable multihitReads so that the naming in stats is nice
  multihitReads <- length(multihitNames) #multihit names is already unique
  stats <- cbind(stats, multihitReads)
  
  ########## IDENTIFY UNIQUELY-PAIRED READS (real sites) ##########  
  allSites <- properlyPairedAlignments[!properlyPairedAlignments$ID %in% multihits$ID]
  
  allSites$breakpoint <- end(flank(allSites, width=-1, both=FALSE,
                                     start=FALSE))
  
  sites.final <- dereplicateSites(allSites)
  
  if(length(sites.final)>0){
    sites.final$clone <- strsplit(names(allSites[1]), "%")[[1]][1]
    sites.final$intLoc <- start(flank(sites.final, width=-1, start=TRUE, both=FALSE))
    sites.final$posid <- paste0(as.character(seqnames(sites.final)),
                               as.character(strand(sites.final)),
                               sites.final$intLoc)
  }
  save(sites.final, file="sites.final.RData")
  save(allSites, file="allSites.RData")
  
  
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
  
  chimeraData <- list("chimeras"=chimeras, "dereplicatedChimeras"=dereplicatedChimeras)
  save(chimeraData, file="chimeraData.RData")
  save(stats, file="stats.RData")
  
}