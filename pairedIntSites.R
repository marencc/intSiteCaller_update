getTrimmedSeqs = function(qualityThreshold, badQuality, qualityWindow, primer, ltrbit, largeLTRFrag, linker, linker_common, mingDNA, read1, read2, alias, vectorSeq){
  
  ##### Load libraries #####
  libs <- c("hiAnnotator", "hiReadsProcessor", "plyr", "ShortRead", "BiocParallel")
  sapply(libs, require, character.only=TRUE)
  #.jinit(parameters="-Xrs")
  
  #source("/home/nerv/Rpackage/processing/hiReadsProcessor.R")
  stats <- data.frame()
  message(alias)
  #workingDir = paste0(masterDir, "/", alias)
  workingDir = alias
  suppressWarnings(dir.create(workingDir, recursive=TRUE))
  setwd(workingDir)
  
  stats.bore = data.frame(sample=alias)
  message("\t read data and make quality heatmaps")  
  
  filenames = list(read1, read2)
  
  reads <- lapply(filenames, function(x) {
    sapply(x, readFastq)
  })
  
  stats.bore$Reads.l.beforeTrim = sum(sapply(reads[[1]], length))#length(reads[[2]])
  stats.bore$Reads.p.beforeTrim = sum(sapply(reads[[2]], length))#(reads[[1]])
  
  r = lapply(reads, function(x){
    seqs = x[[1]]
    if(length(seqs) > 0){
      trimmed = trimTailw(seqs, badQuality, qualityThreshold, round(qualityWindow/2)) #remove anything after 5 bases under Q30 in 10bp window
      trimmed = trimmed[width(trimmed) > 65] #get rid of anything that lost more than half the bases
      if(length(trimmed) > 0){
        trimmedSeqs = sread(trimmed)
        trimmedqSeqs = quality(quality(trimmed))
        names(trimmedSeqs) = names(trimmedqSeqs) = sapply(sub("(.+) .+","\\1",ShortRead::id(trimmed)), function(z){paste0(alias, "%", strsplit(z, "-")[[1]][2])})
      }
    }
    list(trimmedSeqs, trimmedqSeqs)
  })
  
  reads = sapply(r, "[[", 1)
  qualities = sapply(r, "[[", 2)
  R1Quality = qualities[[1]] #this is needed for primerID quality scores later on
  
  #stopifnot(length(reads[[1]]) > 0, length(reads[[2]]) > 0)
  
  stats.bore$Reads.p.afterTrim = length(reads[[2]])
  stats.bore$Reads.l.afterTrim = length(reads[[1]])
  print(stats.bore) 
  
  message("\t trim adaptors")
  
  res.p <- pairwiseAlignSeqs(reads[[2]], patternSeq=primer, 
                             qualityThreshold=1, doRC=F)
  
  reads.p = reads[[2]]
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
    res = primerIDAlignSeqs(subjectSeqs = reads[[1]], patternSeq=linker,
                            doAnchored=T, qualityThreshold1=1, 
                            qualityThreshold2=1, doRC=F)
    res.l <- res[["hits"]]
    res.pID <- res[["primerIDs"]]
  }else{
    res.l <- pairwiseAlignSeqs(subjectSeqs = reads[[1]], patternSeq=linker,
                               qualityThreshold=.95, doRC=F, side="middle")
    start(res.l) = 1
  }
  
  reads.l = reads[[1]]
  if(length(res.l) > 0 ){
    reads.l <- trimSeqs(reads[[1]], res.l, side='left', offBy=1)
    if(grepl("N", linker)){
      R1Quality = R1Quality[match(names(res.pID), names(R1Quality))]
      primerIDs = trimSeqs(reads[[1]], res.pID, side="middle")
      primerIDQuality = subseq(R1Quality, start=start(res.pID), end=end(res.pID))
      primerIDData = list(primerIDs, primerIDQuality)
      
      save(primerIDData, file="primerIDData.RData")
      #primerIDs = trimSeqs(R1Quality, res.pID, side="middle") #would have to change trimSeqs checks
    }
  }
  
  stats.bore$linkered <- length(reads.l)
  
  print(stats.bore) 

  ## check if reads were sequenced all the way by checking for opposite adaptor ##
  message("\t trim opposite side adaptors") #this trimming can fail and we don't care all that much since the DNA molecule could just be super long and we never run into the other side
  
  res.p = NULL
  tryCatch(res.p <- pairwiseAlignSeqs(reads.p, linker_common, qualityThreshold=.55, side='middle', doRC=F), error=function(e){print(paste0("Caught ERROR in pairedIntSites: ", e$message))})
  
  if(!is.null(res.p)){
    end(res.p) = width(reads.p[names(res.p)]) + 1#if we see the common sequence, pitch the rest
    if(length(res.p) > 0 ){
      reads.p <- c(reads.p[!names(reads.p) %in% names(res.p)],
                   trimSeqs(reads.p, res.p, side='right', offBy=1))  
    }
  }
  
  res.l = NULL
  tryCatch(res.l <- pairwiseAlignSeqs(reads.l, largeLTRFrag, qualityThreshold=.55, side='middle', doRC=F), error=function(e){print(paste0("Caught ERROR in pairedIntSites: ", e$message))})
  
  if(!is.null(res.l)){
    end(res.l) = width(reads.l[names(res.l)]) + 1
    if(length(res.l) > 0 ){
      reads.l <- c(reads.l[!names(reads.l) %in% names(res.l)],
                   trimSeqs(reads.l, res.l, side='right', offBy=1))
    }
  }
  
  
  reads.p = subset(reads.p, width(reads.p) > mingDNA)
  reads.l = subset(reads.l, width(reads.l) > mingDNA)
  
  stats.bore$reads.p_afterTrim <- length(reads.p)
  stats.bore$reads.l_afterTrim <- length(reads.l)
  
  
  message("\t trim vector") 
  #we want to do this at the end so that we don't have to worry about partial vector alignments secondary to incomplete trimming of long reads
  
  Vector <- readDNAStringSet(vectorSeq)
  
  blatParameters <- c(minIdentity=70, minScore=5, stepSize=3, 
                      tileSize=8, repMatch=112312, dots=1000, 
                      q="dna", t="dna", out="psl")
  
  findAndRemoveVector <- function(reads, Vector, blatParameters, minLength=10){
    
    hits.p <- read.psl(blatSeqs(query=reads, subject=Vector, 
                                blatParameters=blatParameters, parallel = F), bestScoring=F)
    
    hits.p <- reduce(GRanges(seqnames=hits.p$qName, IRanges(hits.p$qStart, hits.p$qEnd)),
                     min.gapwidth=1200) #collapse instances where a single read has multiple vector alignments
    names(hits.p) = as.character(seqnames(hits.p))
       
    hits.p = hits.p[start(hits.p)<=5 & width(hits.p)>minLength]
    
    reads[!names(reads) %in% names(hits.p)]

  }
  
  tryCatch(reads.p <- findAndRemoveVector(reads.p, Vector, blatParameters=blatParameters), error=function(e){print(paste0("Caught ERROR in pairedIntSites: ", e$message))})
  tryCatch(reads.l <- findAndRemoveVector(reads.l, Vector, blatParameters=blatParameters), error=function(e){print(paste0("Caught ERROR in pairedIntSites: ", e$message))})
  
  stats.bore$reads.p_afterVTrim <- length(reads.p)
  stats.bore$reads.l_afterVTrim <- length(reads.l)
  
  toload = names(reads.p)
  toload <- intersect(names(reads.p), names(reads.l))
  
  stats.bore$reads.lLength = mean(width(reads.l))  
  stats.bore$reads.pLength = mean(width(reads.p))
  
  stats.bore$curated <- length(toload)
  
  
  print(stats.bore)
  stats <- rbind(stats, stats.bore)
  
  #list(reads.p, reads.l, stats)
  if(length(toload) > 0){
    chunks = split(toload, ceiling(seq_along(toload)/20000)) #if using PMACS, cap the number of reads per BLAT thread to 50K since we care about speed rather than number of processers
    for(chunk in names(chunks)){
      toWrite = chunks[[chunk]]
      writeXStringSet(reads.p[names(reads.p) %in% toWrite], file=paste0("p1-", chunk, ".fa"), append=TRUE)
      writeXStringSet(reads.l[names(reads.l) %in% toWrite], file=paste0("p2-", chunk, ".fa"), append=TRUE)
    }
    save(stats, file="stats.RData")
    getwd()
  }else{
    stop("error - no curated reads")
  }
}

processAlignments = function(workingDir, minPercentIdentity, maxAlignStart, maxLength){
  library("sonicLength")
  library("hiAnnotator")
  library("hiReadsProcessor")
  #library("igraph")
  
  #setwd(paste0(masterDir, "/", workingDir))
  setwd(workingDir)
  
  # clean up alignments and prepare for int site calling
  
  processBLATData = function(algns, from){
    names(algns) = algns$qName
    algns$repMatches = NULL
    algns$nCount = NULL
    algns$qNumInsert = NULL
    algns$tNumInsert = NULL
    algns$qName = NULL
    algns$tSize = NULL
    algns$blockCount = NULL
    algns$blockSizes = NULL
    algns$qStarts = NULL
    algns$tStarts = NULL
    algns$from = from
    algns
  }
  
  #source("~/EAS/PMACS_scripts/hiReadsProcessorTemp.R")
  library("BiocParallel")
  library("plyr")
  library("GenomicRanges")
  
  hits.R2 = processBLATData(read.psl(system("ls p1*.fa.psl.gz", intern=T), asGRanges=T, bestScoring=F, removeFile=F), "R2") 
  save(hits.R2, file="hits.R2.RData")
  
  hits.R1 = processBLATData(read.psl(system("ls p2*.fa.psl.gz", intern=T), asGRanges=T, bestScoring=F, removeFile=F), "R1") 
  save(hits.R1, file="hits.R1.RData")
  
  load("stats.RData")
  
  hits.p = append(hits.R1, hits.R2)
  
  readsAligning = length(unique(names(hits.p)))
  
  stats = cbind(stats, readsAligning)
  
  hits.pFilters = (hits.p$matches - hits.p$misMatches - hits.p$tBaseInsert - hits.p$qBaseInsert)/hits.p$qSize
  hits.p = subset(hits.p, hits.pFilters >=.9 & hits.pFilters <= 1.0)
  hits.p$percIdent = 100 * hits.p$matches/hits.p$qSize
  hits.p = subset(hits.p, hits.p$percIdent >= minPercentIdentity & hits.p$qStart <= maxAlignStart)
  
  readsWithGoodAlgnmts = length(unique(names(hits.p)))
  
  stats = cbind(stats, readsWithGoodAlgnmts)
  
  sites.final = GRanges()
  
  allSites = GRanges()
  
  hits.p = flank(hits.p, -1, start=T, both=F)
  
  hits.p = split(hits.p, names(hits.p))
  
  numStrands = as.data.frame.matrix(table(strand(hits.p)))
  
  #numReads = sapply(hits.p, function(x){x$from})
  
  #numStrands = numStrands$"-"==1 & numStrands$"+"==1 #does this also inadvertentaly filter out multihits?
  numStrands = numStrands$"-">=1 & numStrands$"+">=1 #just a quick pre-filter to reduce amount of work in future steps
  
  hits.p = subset(hits.p, numStrands)
  
  reduced = GenomicRanges:::deconstructGRLintoGR(hits.p)
  strand(reduced) = "*" #reduce() doesn't like different strands - we'll add strand info back in later
  reduced = reduce(reduced, min.gapwidth=maxLength, with.revmap=TRUE)
  reduced = GenomicRanges:::reconstructGRLfromGR(reduced, hits.p)
  
  reduced = unlist(reduced)
  reduced$reductionID = c(1:length(reduced)) #names are no longer unique identifier
  hits.p = unlist(hits.p, use.names = FALSE)
  
  lengths = sapply(reduced$revmap, length)
  
  alignmentsFromReads = split(hits.p[unlist(reduced$revmap)]$from, as.vector(Rle(reduced$reductionID, lengths)))
  R1Counts = sapply(alignmentsFromReads, function(x){sum(x=="R1")}) #ordered and named by reductionID
  R2Counts = sapply(alignmentsFromReads, function(x){sum(x=="R2")})
  oneEach = R1Counts==1 & R2Counts==1
  
  hits.pOLD = hits.p
  
  #double check that the subset of reduced really only has one from each read
  reduced2 = subset(reduced, lengths==2)
  reduced2 = subset(reduced2, oneEach[reduced2$reductionID])
  hits.reduced = hits.p[unlist(reduced2$revmap)]
  hits.reduced = subset(hits.reduced, hits.reduced$from=="R2")
  strand(reduced2) = strand(hits.reduced)
  
  #hits.p = reduced2
  
  good = length(unique(names(reduced2)))
  
  stats = cbind(stats, good)
  
  #=======
  
  nonOverlappingSingles = subset(reduced, lengths==1)
  strand(nonOverlappingSingles) = strand(hits.p[unlist(nonOverlappingSingles$revmap)])
  t = table(names(nonOverlappingSingles))
  chimeras = subset(nonOverlappingSingles, names(nonOverlappingSingles) %in% names(subset(t, t==2)))
  chimeras = chimeras[!names(chimeras) %in% names(reduced2)] #not an already-assigned read
  chimeras = split(chimeras, names(chimeras))
  foo = sapply(chimeras, function(x){sum(R1Counts[x$reductionID]) == sum(R2Counts[x$reductionID])})#SLOW!!
  chimeras = subset(chimeras, foo)
  chimera = length(unique(names(chimeras)))
  
  unlistedChimeras = unlist(chimeras, use.names=FALSE)
  reducedChimeras = reduce(unlistedChimeras, min.gapwidth=5)
  unlistedChimeras$chimeraID = subjectHits(findOverlaps(unlistedChimeras, reducedChimeras))
  
  chimeras = split(unlistedChimeras, names(unlistedChimeras))
  uniqueChimeraIDs = unique(sapply(chimeras, function(x){list(x$chimeraID)}))
  uniqueChimeras = length(unique(lapply(uniqueChimeraIDs, sort)))
  
  stats = cbind(stats, chimera)
  
  chimeraData = list("totalReads"=sum(chimera, good), "chimericReads"=chimera, "uniqueChimeras"=uniqueChimeras, "chimeras"=chimeras)
  save(chimeraData, file="chimeraData.RData")
  
  #hits.p = hits.p[lengths==1] #actually remove bad ones (i.e. not properly pairing)
  
  #hits.p = unlist(hits.p)
  
  #strand(hits.p) = as.character(strand(hits.R2[names(hits.p)])) #mark correct strand
  hits.p = reduced2 #from above
  
  hits.p$clone = sapply(strsplit(names(hits.p), "\\."), "[[", 1)
  hits.p$ID = sapply(strsplit(names(hits.p), "%"), "[[", 2)
  #Theoretically shouldn't have multihits across samples, so it's ok to do it here
  multihitNames = unique(names(hits.p[duplicated(hits.p$ID)]))
  multihits = subset(hits.p, names(hits.p) %in% multihitNames)
  multihitReads = length(multihitNames) #multihit names is already unique
  multihitData = list(multihitReads, multihits)
  
  save(multihitData, file="multihitData.RData")
  
  stats = cbind(stats, multihitReads)
  
  uniqueSites = hits.p[!hits.p$ID %in% multihits$ID] #just throwing away multihits for now
  
  #   load("~/Illumina/badSite.RData")
  #   badSites = subsetByOverlaps(uniqueSites, badSite)
  #   save(badSites, file="badSites.RData")
  #   uniqueSites = uniqueSites[!uniqueSites$ID %in% subsetByOverlaps(uniqueSites, badSite)$ID]
  #allSites = uniqueSites #allSites now has reps!
  uniqueSites$breakpoint = end(flank(uniqueSites, width=-1, both=FALSE, start=FALSE))
  
  if(length(uniqueSites)>0){
    allSitesSolostart = uniqueSites
    #end(allSitesSolostart) = start(uniqueSites)
    allSitesSolostart = flank(allSitesSolostart, 5, both=TRUE)
    #sites.final$intLoc = start(flank(sites.final, width=-1, start=TRUE, both=FALSE))
    
    sites.reduced = reduce(allSitesSolostart, min.gapwidth=5, with.revmap=T)
    sites.reduced$counts = sapply(sites.reduced$revmap, length)
    
    allSites = uniqueSites[unlist(sites.reduced$revmap)]
    allSites = split(allSites, Rle(values=c(1:length(sites.reduced)), lengths=sites.reduced$counts))
    allSites = unlist(reduce(allSites, min.gapwidth=5))
    mcols(allSites) = mcols(sites.reduced)
  }
  sites.reduced = allSites
  allSites = uniqueSites
  
  #sites.reduced$counts <- countOverlaps(sites.reduced, uniqueSites)
  sites.final = sites.reduced
  if(length(sites.reduced)>0){
    sites.reduced$clone = strsplit(names(uniqueSites[1]), "\\.")[[1]][1]
    sites.final = sites.reduced
    sites.final$intLoc = start(flank(sites.final, width=-1, start=TRUE, both=FALSE))
    #sites.final$Aliasposid = paste0(sites.final$clone, "_", as.character(seqnames(sites.final)), as.character(strand(sites.final)), sites.final$intLoc)
    sites.final$posid = paste0(as.character(seqnames(sites.final)), as.character(strand(sites.final)), sites.final$intLoc)
    sites.final$estAbund = NA #tmp for now, will be overwritten if sonicAbund is calculated
  }
  save(sites.final, file="sites.final.RData")
  save(stats, file="stats.RData")
  save(allSites, file="allSites.RData")  
}