getTrimmedSeqs = function(qualityThreshold, badQuality, qualityWindow, primer, ltrbit, largeLTRFrag, linker, hasPrimerID, linker_common, mingDNA, read1, read2, paired, alias, vectorSeq, PMACS){
  
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
  
  filenames = list(read2)
  
  if(paired){
    filenames = list(read1, read2)
  }
  
  reads <- lapply(filenames, function(x) {
    sapply(x, readFastq)
      #bore <- readFastq(y)
      
      #qual <- quality(bore)
      #qual <- FastqQuality(srsort(qual[sample(length(qual),100)]))
      #pheatmap(as(qual,'matrix'), cluster_rows=F, cluster_cols=F, 
      #         filename=paste0(gsub("*.fastq.*","",x),".qualityBeforeTrimming.pdf"), width=11, height=8.5)  
      
      #rm(qual)
      
      #bore
    #})
  })
  
  if(paired){
    stats.bore$Reads.l.beforeTrim = sum(sapply(reads[[1]], length))#length(reads[[2]])
    stats.bore$Reads.p.beforeTrim = sum(sapply(reads[[2]], length))#(reads[[1]])
  }else{
    stats.bore$Reads.p.beforeTrim = sum(sapply(reads[[1]], length))#(reads[[1]]) #hasn't been moved to reads[[2]] yet 
  }
  
  r = lapply(reads, function(x){
    res = DNAStringSet()
    qres = BStringSet()
    for(y in 1:length(x)){
      if(length(x[[y]]) > 0){
        trimmed = trimTailw(x[[y]], badQuality, qualityThreshold, round(qualityWindow/2)) #remove anything after 5 bases under Q30 in 10bp window
        trimmed = trimmed[width(trimmed) > (max(width(x[[y]]))/2)] #get rid of anything that lost more than half the bases
        if(length(trimmed) > 0){
          trimmedSeqs = sread(trimmed)
          trimmedqSeqs = quality(quality(trimmed))
          n = sapply(sub("(.+) .+","\\1",ShortRead::id(trimmed)), function(z){paste0(alias, ".", y, "%", strsplit(z, "-")[[1]][2])})
          names(trimmedSeqs) = names(trimmedqSeqs) = n
        #  names(trimme)
          res = append(res, trimmedSeqs)
          qres = append(qres, trimmedqSeqs)
        }
      }
    }
    list(res, qres)
  })
  
  reads = sapply(r, "[[", 1)
  qualities = sapply(r, "[[", 2)
  R1Quality = qualities[[1]]
  
  if(!paired){
    reads = list(NULL, reads[[1]]) #later we refer to reads[[2]] as the primed (i.e. LTR read) - easier to just modify here
  }
  
  #stopifnot(length(reads[[1]]) > 0, length(reads[[2]]) > 0)
  
  stats.bore$Reads.p.afterTrim = length(reads[[2]])
  if(paired){
    stats.bore$Reads.l.afterTrim = length(reads[[1]])
  }
  print(stats.bore) 
  
  message("\t trim adaptors") #if we don't see any of the correct adaptors, we should crash
  #turn these into sapplys?
  res.p = NULL
  for(x in primer){
    res.p <- append(res.p, pairwiseAlignSeqs(reads[[2]], patternSeq=x, 
                              qualityThreshold=1, doRC=F))#, mc.cores=8)
  }

  reads.p = reads[[2]]
  if(length(res.p) > 0){
    reads.p <- trimSeqs(reads[[2]], res.p, side='left', offBy=1)
  }
  
  stats.bore$primed <- length(reads.p)

  res.ltr = NULL
  for(x in ltrbit){
    res.ltr <- append(res.ltr, pairwiseAlignSeqs(reads.p, patternSeq=x, 
                                 qualityThreshold=1, doRC=F))
  }
  if(length(res.ltr) > 0 ){
    reads.p <- trimSeqs(reads.p, res.ltr, side='left', offBy=1)
  }
  
  stats.bore$LTRed <- length(reads.p)
  
  if(paired){
    res.l = NULL
    res.pID = NULL
    for(x in linker){
      if(hasPrimerID){
        res = primerIDAlignSeqs(subjectSeqs = reads[[1]], patternSeq=x,
                                              doAnchored=T, qualityThreshold1=1, 
                                              qualityThreshold2=1, doRC=F)
        res.l <- append(res.l, res[["hits"]])
        res.pID <- append(res.pID, res[["primerIDs"]])
      }else{
        res.l <- append(res.l, pairwiseAlignSeqs(subjectSeqs = reads[[1]], patternSeq=x,
                                                 qualityThreshold=.95, doRC=F, side="middle"))
        start(res.l) = 1
      }
    }
    reads.l = reads[[1]]
    if(length(res.l) > 0 ){
      reads.l <- trimSeqs(reads[[1]], res.l, side='left', offBy=1)
      if(hasPrimerID){
        R1Quality = R1Quality[match(names(res.pID), names(R1Quality))]
        primerIDs = trimSeqs(reads[[1]], res.pID, side="middle")
        primerIDQuality = subseq(R1Quality, start=start(res.pID), end=end(res.pID))
        primerIDData = list(primerIDs, primerIDQuality)
        
        save(primerIDData, file="primerIDData.RData")
        #primerIDs = trimSeqs(R1Quality, res.pID, side="middle") #would have to change trimSeqs checks
      }
    }

    stats.bore$linkered <- length(reads.l)
  }
  
  print(stats.bore) 
  
  message("\t trim vector") #this can fail and we can probably limp through the rest of the analysis
  
  Vector <- readDNAStringSet(vectorSeq)
  
  blatParameters <- c(minIdentity=70, minScore=5, stepSize=3, 
                      tileSize=8, repMatch=112312, dots=1000, 
                      q="dna", t="dna", out="psl")

  findAndRemoveVector <- function(reads, Vector, blatParameters, read,
                                  minLength=10) {

    hits.p <- read.psl(blatSeqs(query=reads, subject=Vector, 
                                blatParameters=blatParameters, parallel = F), bestScoring=F)

    hits.p <- reduce(GRanges(seqnames=hits.p$qName, IRanges(hits.p$qStart, hits.p$qEnd)),
                     min.gapwidth=1200)
    
    hits.p.r <- ranges(hits.p)
    names(hits.p.r) <- as.character(seqnames(hits.p))
    
    minAlignLength = minLength
    if(read=="l"){
      minAlignLength = 34
    }
    
    hits.p.r = hits.p.r[start(hits.p.r)<=5 & width(hits.p.r)>minAlignLength] #anything that doesn't fit these params
    #is likely hitting the LTR read-through (maybe change this number to mingDNA?) and should be 
    #filtered out later
    reads[!names(reads) %in% names(hits.p.r)]
    #trimmed <- trimSeqs(reads, hits.p.r, side=side, offBy = 0)
    #trimmed <- trimmed[width(trimmed)>=minLength]
    #c(reads[!names(reads) %in% names(hits.p.r)], trimmed)
  }
  
  tryCatch(reads.p <- findAndRemoveVector(reads.p, Vector, blatParameters=blatParameters, "p"), error=function(e){print(paste0("Caught ERROR in pairedIntSites: ", e$message))})
  
  if(paired){
    tryCatch(reads.l <- findAndRemoveVector(reads.l, Vector, blatParameters=blatParameters, "l"), error=function(e){print(paste0("Caught ERROR in pairedIntSites: ", e$message))})
  }
  stats.bore$reads.p_afterVTrim <- length(reads.p)
  if(paired){
    stats.bore$reads.l_afterVTrim <- length(reads.l)
  }
  
  ## check if reads were sequenced all the way by checking for opposite adaptor ##
  message("\t trim opposite side adaptors") #this trimming can fail and we don't care all that much since the DNA molecule could just be super long and we never run into the other side
  
  res.p = NULL
  for(x in linker_common){
    tryCatch(res.p <- append(res.p, pairwiseAlignSeqs(reads.p, x, qualityThreshold=.55, side='middle', doRC=F)), error=function(e){print(paste0("Caught ERROR in pairedIntSites: ", e$message))})
  }  
  if(!is.null(res.p)){
    end(res.p) = width(reads.p[names(res.p)]) + 1#if we see the common sequence, pitch the rest
    if(length(res.p) > 0 ){
      reads.p <- c(reads.p[!names(reads.p) %in% names(res.p)],
                   trimSeqs(reads.p, res.p, side='right', offBy=1))  
    }
  }
  
  if(paired){
    res.l = NULL
    for(x in largeLTRFrag){
      tryCatch(res.l <- append(res.l, pairwiseAlignSeqs(reads.l, x, qualityThreshold=.55, side='middle', doRC=F)), error=function(e){print(paste0("Caught ERROR in pairedIntSites: ", e$message))})
    }
    if(!is.null(res.l)){
      end(res.l) = width(reads.l[names(res.l)]) + 1
      if(length(res.l) > 0 ){
        reads.l <- c(reads.l[!names(reads.l) %in% names(res.l)],
                     trimSeqs(reads.l, res.l, side='right', offBy=1))
      }
    }
  }
  
  reads.p = subset(reads.p, width(reads.p) > mingDNA)
  stats.bore$reads.p_afterTrim <- length(reads.p)
  toload = names(reads.p)
  
  if(paired){
    reads.l = subset(reads.l, width(reads.l) > mingDNA)
    stats.bore$reads.l_afterTrim <- length(reads.l)
    toload <- intersect(names(reads.p), names(reads.l))
    stats.bore$reads.lLength = mean(width(reads.l))
  }
  
  stats.bore$curated <- length(toload)
  stats.bore$reads.pLength = mean(width(reads.p))
  
  print(stats.bore)
  stats <- rbind(stats, stats.bore)
  
  #list(reads.p, reads.l, stats)
  if(length(toload) > 0){
    if(PMACS){
      chunks = split(toload, ceiling(seq_along(toload)/20000)) #if using PMACS, cap the number of reads per BLAT thread to 50K since we care about speed rather than number of processers
      for(chunk in names(chunks)){
        toWrite = chunks[[chunk]]
        writeXStringSet(reads.p[names(reads.p) %in% toWrite], file=paste0("p1-", chunk, ".fa"), append=TRUE)
        if(paired){
          writeXStringSet(reads.l[names(reads.l) %in% toWrite], file=paste0("p2-", chunk, ".fa"), append=TRUE)
        }
      }
    }else{
      writeXStringSet(reads.p[names(reads.p) %in% toload], file="p1.fa", append=TRUE)
      if(paired){
        writeXStringSet(reads.l[names(reads.l) %in% toload], file="p2.fa", append=TRUE)
      }
    }
    save(stats, file="stats.RData")
    getwd()
  }else{
    stop("error - no curated reads")
  }
}

alignTrimmedSeqs = function(alignMethod, workingDir, indexPath, aliases, mingDNA, paired, maxCores){
  library("Rsubread")
  library("Rsamtools")
  library("doMC")
  library("hiReadsProcessor")
  message("\t align pairs to hg18")
  setwd(workingDir)
  
  scanParams = ScanBamParam(flag = scanBamFlag(isUnmappedQuery=FALSE), what="qwidth")

  mdata = list("p1")
  conversion = list("R2")
  
  if(any(paired)){
    mdata = append(mdata, "p2")
    conversion = append(conversion, "R1")
  }
  
  names(conversion) = unlist(mdata)
  
  ports = list(5561, 5562)
  names(ports) = unlist(mdata)
  
  alignments = mclapply(mdata, function(x){
    if(alignMethod=="RSubread"){
      Rsubread::align(index=indexPath, readfile1=paste0(x, ".fa"), 
            output_file=paste0(x, ".alignment.bam"), nthreads=2, nBestLocations=10,
            minFragLength=0, maxFragLength=1200, output_format="BAM", unique=FALSE, tieBreakHamming=FALSE, tieBreakQS=FALSE)
      sortBam(paste0(x, ".alignment.bam"), paste0(x, ".alignment.sorted"), maxMemory=2048) # sort entries in bam file
      indexBam(paste0(x, ".alignment.sorted.bam"))
      algns = readGAlignmentsFromBam(paste0(x, ".alignment.sorted.bam"), index=paste0(x, ".alignment.sorted.bam"), use.names=TRUE, param=scanParams)
      algns.gr = granges(algns)
      algns.gr$from = conversion[[x]]
      algns.gr$cigar = cigar(algns)
      algns.gr$qWidth = qwidth(algns)
      algns = algns.gr
    }
    if(alignMethod=="BLAT"){
      registerDoMC(4) #server chokes on more than 4 cores
      algns = read.psl(blatSeqs(query=paste0(x, ".fa"), subject="~/Indexes/hg18.2bit", standaloneBlat=F, host="localhost", port=ports[[x]], blatParameters=c(minIdentity = 85, minScore = round(.9*min(mingDNA)), tileSize = 11, repMatch = 112312, dots = 1000, maxDnaHits="10", q = "dna", t = "dna", out = "psl")), bestScoring=FALSE, asGRanges=TRUE, removeFile=FALSE)
      #hg18 = readDNAStringSet("~/Indexes/hg18.fa")
      #algns = read.psl(blatSeqs(query=paste0(x, ".fa"), subject=hg18, blatParameters=c(minIdentity = 80, minScore = round(.9*min(mingDNA)), stepSize = 5, tileSize = 11, repMatch = 112312, dots = 1000, q = "dna", t = "dna", out = "psl")), bestScoring=FALSE, asGRanges=TRUE)
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
      algns$from = conversion[[x]]
    }
    algns = split(algns, sapply(strsplit(names(algns), "\\."), "[[", 1))#split by clone, not rep - reps will be identified in processAlignment
  }, mc.cores=2)
  
  mclapply(aliases, function(x){
    #hits = alignments[[1]][[x]]
    hits.R2 = alignments[[1]][[x]]
    #f = paste0(masterDir, "/", names(alignments[[1]][x]), "/alignment.RData")
    f = paste0(x, "/alignment.RData")
        
    hits.R1 = NULL
    if(paired[match(x, aliases)]){
      hits.R1 = alignments[[2]][[x]]
    }
    #hits = append(hits.R1, hits.R2)
    
    save(hits.R1, hits.R2, file=f)
    
    #names(alignments[[1]][x]) #return names of aligned samples
    x
  }, mc.cores=min(6, length(aliases)))
  
}

processAlignments = function(alignMethod, PMACS, workingDir, minPercentIdentity, maxAlignStart, maxLength, hasPrimerID, paired, removeMultihits, csvOut){
  library("sonicLength")
  library("hiAnnotator")
  library("hiReadsProcessor")
  #library("igraph")
  
  #setwd(paste0(masterDir, "/", workingDir))
  setwd(workingDir)
  if(PMACS){
    
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
    
    if(paired){
      hits.R1 = processBLATData(read.psl(system("ls p2*.fa.psl.gz", intern=T), asGRanges=T, bestScoring=F, removeFile=F), "R1") 
      save(hits.R1, file="hits.R1.RData")
    }
  }else{
    load("alignment.RData")
  }
  load("stats.RData")
  
  hits.p = hits.R2
  
  if(paired){
    hits.p = append(hits.R1, hits.p)
  }
  
  readsAligning = length(unique(names(hits.p)))
  
  stats = cbind(stats, readsAligning)
  
  if(alignMethod=="RSubread"){
    hits.p$percIdent = 100 * cigarOpTable(hits.p$cigar)[,"M"]/hits.p$qWidth
    
    alignStarts = sapply(explodeCigarOpLengths(hits.p$cigar), "[[", 1)
    alignStarts[sapply(explodeCigarOps(hits.p$cigar), "[[", 1) != "S"] = 0
    hits.p$alignStart = alignStarts
    
    hits.p = subset(hits.p, percIdent >= minPercentIdentity & alignStart <= maxAlignStart)
  }
  if(alignMethod=="BLAT"){
    hits.pFilters = (hits.p$matches - hits.p$misMatches - hits.p$tBaseInsert - hits.p$qBaseInsert)/hits.p$qSize
    hits.p = subset(hits.p, hits.pFilters >=.9 & hits.pFilters <= 1.0)
    hits.p$percIdent = 100 * hits.p$matches/hits.p$qSize
    hits.p = subset(hits.p, hits.p$percIdent >= minPercentIdentity & hits.p$qStart <= maxAlignStart)
  }
  
  readsWithGoodAlgnmts = length(unique(names(hits.p)))
  
  stats = cbind(stats, readsWithGoodAlgnmts)
  
  sites.final = GRanges()
  
  allSites = GRanges()
  
  if(paired){
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
  }
  
  #hits.p = hits.p[lengths==1] #actually remove bad ones (i.e. not properly pairing)
  
  #hits.p = unlist(hits.p)
  
  #strand(hits.p) = as.character(strand(hits.R2[names(hits.p)])) #mark correct strand
  if(paired){
    hits.p = reduced2 #from above
  }
  hits.p$clone = sapply(strsplit(names(hits.p), "\\."), "[[", 1)
  hits.p$rep = sapply(strsplit(sapply(strsplit(names(hits.p), "%"), "[[", 1), "\\."), "[[", 2)
  hits.p$ID = sapply(strsplit(names(hits.p), "%"), "[[", 2)
  #Theoretically shouldn't have multihits across samples, so it's ok to do it here
  multihitNames = unique(names(hits.p[duplicated(hits.p$ID)]))
  multihits = subset(hits.p, names(hits.p) %in% multihitNames)
  multihitReads = length(multihitNames) #multihit names is already unique
  multihitData = list(multihitReads, multihits)
  
  save(multihitData, file="multihitData.RData")

  stats = cbind(stats, multihitReads)
  
  uniqueSites = hits.p
  if(removeMultihits){
    uniqueSites = hits.p[!hits.p$ID %in% multihits$ID] #just throwing away multihits for now
  }
  
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
  
#   if(hasPrimerID & paired & length(uniqueSites)>0){
#     load("primerIDData.RData")
#     
#     primerIDs = primerIDData[[1]]
#     primerIDquals = primerIDData[[2]]
#     #mat <- nucleotideSubstitutionMatrix(match = 1, mismatch = 0, baseOnly = TRUE)
#     matchCoords = match(names(uniqueSites), names(primerIDs))
#     primerIDs = primerIDs[matchCoords]
#     #primerIDquals = primerIDquals[matchCoords]
#     
#     splitGrp = rep(0, length(primerIDs))
#     for(i in c(1:length(sites.reduced))){
#       splitGrp[sites.reduced$mapping[[i]]] = i
#     }
#     
#     splitPrimerIDs = split(primerIDs, splitGrp)
#     #splitPrimerIDquals = split(primerIDquals, splitGrp)
#     
#     dfr = NULL
#     for(i in 1:length(splitPrimerIDs)){ #i is siteID
#       x = splitPrimerIDs[[i]]
#       len = length(x)
#       if(len>1){
#         mat = as.matrix(stringDist(x))
#         mat[lower.tri(mat, diag=TRUE)] = 12
#         toCheck = which(mat<=1, arr.ind=TRUE)
#         if(nrow(toCheck)>0){
#           rownames(toCheck) = NULL
#           seqNames = names(x) #order preserved in mat
#           net = graph.edgelist(toCheck, directed=FALSE)
#           clusts = clusters(net)
#           clusterIDs = clusts[["membership"]]
#           if(length(clusterIDs) < length(seqNames)){ #last seqname wasn't involved in any cluster, so it was never inferred
#             clusterIDs = append(clusterIDs, c((max(clusterIDs)+1):(max(clusterIDs)+length(seqNames)-length(clusterIDs))))
#           }
#         }
#         else{clusterIDs = c(1:nrow(mat))} #just fudge it so that each item is it's own 'cluster'
#         splitLengths = split(width(allSites[seqNames]), clusterIDs)
#         tmp = data.frame("clusterID"=paste0(i, ".", as.character(Rle(values=names(splitLengths), lengths=sapply(splitLengths, length)))), "width"=unlist(splitLengths, use.names=FALSE))
#         dfr = rbind(dfr, droplevels(unique(tmp)))
#       }
#       else{dfr=rbind(dfr, data.frame("clusterID"=paste0(i,".",1), "width"=width(allSites[names(x)])))}
#     }
#     rownames(dfr)=NULL
#     res = estAbund(dfr$clusterID, dfr$width)$theta
#     splitBy = as.integer(sapply(strsplit(names(res), "\\."), "[[", 1))
#     splitAbunds = split(res, splitBy)
#     pidAbunds = round(sapply(splitAbunds, sum))
#     sites.final$myPidAbund = pidAbunds
#   }
    
  
#   sites.reduced$cdHitspidAbund = sapply(splitPrimerIDs, function(x){
#     if(length(x)>1){ #hopfully save some computation?
#       #print("a")
#       #allPIDs = primerIDs[names(uniqueSites[x])]
#       #allPIDQs = primerIDquals[names(uniqueSites[x])]
#       #writeFasta(allPIDs, paste0(workingDir, "/in.fasta"))
#       writeFasta(x, paste0(workingDir, "/in.fasta"))
#       
#       #system("./cd-hit -i ./in.fasta -o ./out.fasta -T 7 -c .85 &> /dev/null")
#       system("cd-hit-est -i in.fasta -o out.fasta -c .85", ignore.stdout=TRUE, wait=TRUE)
#       length(readFasta(paste0(workingDir, "/out.fasta")))
#     }
#     else{
#       #print("b")
#       1 }
#   })
#   
  
#   sites.reduced$mypidAbund = sapply(c(1:length(sites.reduced)), function(x){
#     allPIDs = splitPrimerIDs[[x]]
#     allPIDQs = splitPrimerIDquals[[x]]
#     #PIDdata = lapply(c(1:length(allPIDs)), function(y){list(allPIDs[[y]], allPIDQs[[y]])})
#     clusters = list()
#     while(length(allPIDs) > 0){
#       seed = allPIDs[1]
#       seedQ = allPIDQs[1]
#       if(length(allPIDs)>1){ #if it's =0, it won't get through the while loop, if it's =1, then we just added last one to cluster
#         allPIDs = allPIDs[2:length(allPIDs)]
#         allPIDQs = allPIDQs[2:length(allPIDQs)]
#         
#         algnmnt = pairwiseAlignment(allPIDs, seed, PhredQuality(allPIDQs), PhredQuality(seedQ), scoreOnly=TRUE) > 10
#         seed = append(seed, allPIDs[algnmnt])
#         allPIDs = allPIDs[!algnmnt]
#         allPIDQs = allPIDQs[!algnmnt]
#       }
#       else{allPIDs = NULL} #it's finished, and the remaining value was stored in seed - this acts as a 'break' but allows seed to be added to clusters (next line)
#       clusters = append(clusters, list(seed))
#     }
#     length(clusters)
#   })
  
  #sites.reduced$counts <- countOverlaps(sites.reduced, uniqueSites)
  sites.final = sites.reduced
  if(length(sites.reduced)>0){
    sites.reduced$clone = strsplit(names(uniqueSites[1]), "\\.")[[1]][1]
    sites.final = sites.reduced
    sites.final$intLoc = start(flank(sites.final, width=-1, start=TRUE, both=FALSE))
    #sites.final$Aliasposid = paste0(sites.final$clone, "_", as.character(seqnames(sites.final)), as.character(strand(sites.final)), sites.final$intLoc)
    sites.final$posid = paste0(as.character(seqnames(sites.final)), as.character(strand(sites.final)), sites.final$intLoc)
    sites.final$estAbund = NA #tmp for now, will be overwritten if sonicAbund is calculated

#========SONIC LENGTH BELOW============
#     
    if(paired & removeMultihits){ #will get lots of ugly NA's on non-multihit-filtered because we can't gaurantee that 1 alignment = 1 read
            
#       allSites = allSites[unlist(sites.final$mapping)]
#       
#       allSites$realStart = rep(sites.final$intLoc, sites.final$counts)
      allSites$realStart[unlist(sites.final$revmap)] = rep(sites.final$intLoc, sites.final$counts)
      
      if(workingDir != "~/Illumina/finalOpt/Lib4_4"){
      #if(FALSE){
      
        sonicSites = allSites
        
        pos = paste0(as.character(seqnames(sonicSites)), as.character(strand(sonicSites)), sonicSites$realStart)
        
        if(any(allSites$rep>1)){ #has reps
          dfr = data.frame("posid"=pos, "qEnd"=width(sonicSites), "rep"=sonicSites$rep)
        }  else{
          dfr = data.frame("posid"=pos, "qEnd"=width(sonicSites))
        }
        
        
        #dfr$Aliasposid = paste0(dfr$Alias, "_", as.character(seqnames(sonicSites)), as.character(strand(sonicSites)), pos)
        #dfr$posid = paste0(as.character(seqnames(sonicSites)), as.character(strand(sonicSites)), sonicSites$realStart)
        dfr = droplevels(unique(dfr))
        
        if(any(allSites$rep>1)){ #has reps
          res = estAbund(dfr$posid, dfr$qEnd, dfr$rep)
        }  else{
          res = estAbund(dfr$posid, dfr$qEnd)
        }
        
        allAbund = round(res$theta)
        allAbundVar = res$var.theta
        
        #sites.final = unlist(sites.final, use.names=FALSE)
        #sites.final$estAbund = allAbund[sites.final$Aliasposid]
        #sites.final$estAbund = allAbund[sites.final$Aliasposid]
        sites.final$estAbund = allAbund[sites.final$posid]
        #sites.final$varEstAbund = allAbundVar[sites.final$Aliasposid]
      }
      
    }


#===========END SONIC ABUND==============
  }
  save(sites.final, file="sites.final.RData")
  save(stats, file="stats.RData")
  save(allSites, file="allSites.RData")
  
  if(csvOut & length(sites.final) > 0){
    genes <- dbGetQuery(connectToIntSitesDB(), "select distinct geneName, Chrom, strand, txStart, txEnd from hg18.refflat")
    genes <- makeGRanges(genes)
    sites.annotated = sites.final
    sites.annotated$revmap = NULL
    sites.annotated <- getSitesInFeature(sites.annotated, genes, colnam="inRefGene")
    sites.annotated <- getNearestFeature(sites.annotated, genes, colnam="nrstRefGene")
    sites.annotated$geneType <- with(mcols(sites.annotated),
                                     ifelse(inRefGene=="FALSE",
                                            paste0("~",as.character(nrstRefGene)),
                                            as.character(inRefGene)))
    
    out.df = as.data.frame(sites.annotated)
    out.df$strand = NULL
    write.table(out.df, file="annotation.csv", sep=",", row.names=F) 
  }
  
}

getIntSitesFromSeqs = function(parameters, alignMethod, alignDir, maxCores){
  library(parallel)
  library(ShortRead)
  
  suppressWarnings(dir.create(alignDir, recursive=TRUE))

  PMACS = F
  
  aliases = unlist(mclapply(parameters, function(x){
    tryCatch(eval(as.call(append(getTrimmedSeqs, append(x[1:15], PMACS))), error=function(e){print(paste0("ERROR: ", e$message))}))
  }, mc.cores = min(maxCores, length(parameters))))
  
  aliasErrors = aliases[grepl("ERROR: ", aliases)]
  aliases = aliases[!grepl("ERROR: ", aliases)]
  stopifnot(length(aliases) > 0)
  save(aliases, file=paste0(alignDir, "/aliases.RData"))
  save(aliasErrors, file=paste0(alignDir, "/aliasErrors.RData"))

  setwd(alignDir)
#  load("aliases.RData")
#  load("aliasErrors.RData")
  source("/media/8TB_PLAYGROUND/home/ericsherman/Illumina/alignSites.R")

  alignment = unlist(NewalignTrimmedSeqs(alignMethod, alignDir, "/home/ericsherman/Indexes/subread/hg18_noRandom/hg18_noRandom", aliases, sapply(parameters[match(aliases, unlist(sapply(parameters, "[[", 14)))], "[[", 10), sapply(parameters[match(aliases, unlist(sapply(parameters, "[[", 14)))], "[[", 13), maxCores))

  save(alignment, file=paste0(alignDir, "/alignment.RData"))
  
  foo = mclapply(aliases, function(x){
    eval(as.call(append(list(processAlignments, alignMethod, PMACS), parameters[[match(x, unlist(sapply(parameters, "[[", 14)))]][c(14,16:18,8,13,19,20)])))
  }, mc.cores = min(maxCores, length(alignment)))
  
  save(foo, file=paste0(alignDir, "/foo.RData"))
  
  return(NULL)
}


condense = function(workingDir, aliases, condensedAlias){
  library("sonicLength")
  library("GenomicRanges")
  library("ShortRead")
  suppressWarnings(dir.create(paste0(workingDir, "/condensed/")))
  as = GRanges()
  sites.final = GRanges()
  stats = NULL
  pidData = list(DNAStringSet(), BStringSet())
  for(alias in aliases){
    rep = strsplit(alias, "-")[[1]][[length(strsplit(alias, "-")[[1]])]]
    if(file.exists(paste0(workingDir, "/", alias, "/allSites.RData")) & length(get(load(paste0(workingDir, "/", alias, "/allSites.RData")))) > 1){
      allSites = get(load(paste0(workingDir, "/", alias, "/allSites.RData")))
      allSites$rep = rep
      as = append(as, allSites)
      primerIDData = get(load(paste0(workingDir, "/", alias, "/primerIDData.RData")))
      pidData[[1]] = append(pidData[[1]], primerIDData[[1]])
      pidData[[2]] = append(pidData[[2]], primerIDData[[2]])
      stats = rbind(stats, get(load(paste0(workingDir, "/", alias, "/stats.RData"))))
    }
  }
  suppressWarnings(dir.create(paste0(workingDir, "/condensed/", condensedAlias), recursive=FALSE))
  allSites = as
  save(allSites, file=paste0(workingDir, "/condensed/", condensedAlias, "/allSites.RData"))
  primerIDData = pidData
  save(primerIDData, file=paste0(workingDir, "/condensed/", condensedAlias, "/primerIDData.RData"))
  stats = colSums(Filter(is.numeric, stats))
  save(stats, file=paste0(workingDir, "/condensed/", condensedAlias, "/stats.RData"))
  
  uniqueSites = allSites
  allSitesSolostart = flank(uniqueSites, -1, start=T, both=F)
  allSitesSolostart = flank(allSitesSolostart, 5, both=TRUE)      
  sites.reduced = reduce(allSitesSolostart, min.gapwidth=5, with.revmap=T)
  sites.reduced$counts = sapply(sites.reduced$revmap, length)  
  
  allSites = uniqueSites[unlist(sites.reduced$revmap)]
  allSites = split(allSites, Rle(values=c(1:length(sites.reduced)), lengths=sites.reduced$counts))
  allSites = unlist(reduce(allSites, min.gapwidth=5))
  mcols(allSites) = mcols(sites.reduced)
  sites.final = allSites
  allSites = uniqueSites
  
  sites.final$intLoc = start(flank(sites.final, width=-1, start=TRUE, both=FALSE))
  #sites.final$Aliasposid = paste0(sites.final$clone, "_", as.character(seqnames(sites.final)), as.character(strand(sites.final)), sites.final$intLoc)
  sites.final$posid = paste0(as.character(seqnames(sites.final)), as.character(strand(sites.final)), sites.final$intLoc)
  sites.final$estAbund = NA #tmp for now, will be overwritten if sonicAbund is calculated
  
  allSites$realStart[unlist(sites.final$revmap)] = rep(sites.final$intLoc, sites.final$counts)

  #if(workingDir != "~/Illumina/finalOpt/Lib4_4"){
  if(TRUE){
    
    sonicSites = allSites
    
    pos = paste0(as.character(seqnames(sonicSites)), as.character(strand(sonicSites)), sonicSites$realStart)
    
    if(any(allSites$rep>1)){ #has reps
      dfr = data.frame("posid"=pos, "qEnd"=width(sonicSites), "rep"=sonicSites$rep)
    }  else{
      dfr = data.frame("posid"=pos, "qEnd"=width(sonicSites))
    }
    
    
    #dfr$Aliasposid = paste0(dfr$Alias, "_", as.character(seqnames(sonicSites)), as.character(strand(sonicSites)), pos)
    dfr$posid = paste0(as.character(seqnames(sonicSites)), as.character(strand(sonicSites)), sonicSites$realStart)
    dfr = droplevels(unique(dfr))
    
   if(any(allSites$rep>1)){ #has reps
     res = estAbund(dfr$posid, dfr$qEnd, dfr$rep)
   }  else{
     res = estAbund(dfr$posid, dfr$qEnd)
   }
    
    allAbund = round(res$theta)
    #allAbundVar = res$var.theta
    
    #sites.final = unlist(sites.final, use.names=FALSE)
    #sites.final$estAbund = allAbund[sites.final$Aliasposid]
    #sites.final$estAbund = allAbund[sites.final$Aliasposid]
    sites.final$estAbund = allAbund[sites.final$posid]
    #sites.final$varEstAbund = allAbundVar[sites.final$Aliasposid]
  }
  
  
  save(sites.final, file=paste0(workingDir, "/condensed/", condensedAlias, "/sites.final.RData"))
  
  
}

makeChuckData = function(aliases){
  as = sf = GRanges()
  pids = DNAStringSet()
  for(sample in aliases){
    load(paste0(sample, "/allSites.RData"))
    load(paste0(sample, "/sites.final.RData"))
    load(paste0(sample, "/primerIDData.RData"))
    
    as = append(as, allSites)
    sf = append(sf, sites.final)
    pids = append(pids, primerIDData[[1]])
  }
  primerIDs = pids[match(names(as), names(pids))]
  rm(pids)
  
  foo = DataFrame("rep"=Rle(as$rep), "start"=start(as), "end"=end(as), "strand"=strand(as), "chr"=seqnames(as))
  n = sapply(strsplit(as$clone, "/"), "[[", 4)
  n = strsplit(n, "-")
  n = sapply(n, function(x){paste0(x[1:(length(x)-1)], collapse="-")})
  foo$sample = Rle(n)
  mcols(primerIDs) = foo
  foo
}

  