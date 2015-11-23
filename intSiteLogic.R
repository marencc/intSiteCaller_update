## load hiReadsProcessor.R
libs <- c("plyr", "BiocParallel", "Biostrings", "GenomicAlignments" ,"hiAnnotator" ,"sonicLength", "GenomicRanges", "BiocGenerics", "ShortRead", "GenomicRanges", "igraph")
null <- suppressMessages(sapply(libs, library, character.only=TRUE))

codeDir <- get(load("codeDir.RData"))

stopifnot(file.exists(file.path(codeDir, "hiReadsProcessor.R")))
source(file.path(codeDir, "hiReadsProcessor.R"))
source(file.path(codeDir, "standardization_based_on_clustering.R"))


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


#' as tittled make PairwiseAlignmentsSingleSubject easily accessable
#' as needed by other functions
PairwiseAlignmentsSingleSubject2DF <- function(PA, shift=0) {
    stopifnot("PairwiseAlignmentsSingleSubject"  %in% class(PA))
    
    return(data.frame(
        width=width(pattern(PA)),
        score=score(PA),
        mismatch=width(pattern(PA))-score(PA),
        start=start(pattern(PA))+shift,
        end=end(pattern(PA))+shift
        ))
}


#' subset and substring
#' trim primer and ltrbit off of ltr side of read, R2 in protocol
#' both primer and ltrbit are required, otherwise disgard it
#' allow 2 mismatch for either primer or ltrbit
#' runSeq <- sapply(1:10000, function(i)
#'                  paste(sample(c("A","C","G","T"), 8, replace=TRUE),
#'                        collapse=""))
#' runSeq.p <- pairwiseAlignment(pattern=runSeq,
#'                               subject=primer,
#'                               substitutionMatrix=submat1,
#'                               gapOpening = 0,
#'                               gapExtension = 1,
#'                               type="overlap")
#' runSeq.p.df <- PairwiseAlignmentsSingleSubject2DF(runSeq.p)
#' table(runSeq.p.df$score)
#'   1    2    3    4    5    6    7
#' 211 3338 4330 1751  344   25    1
#' false positive rate 0.0025 and thus maxMisMatch=2,
#' (1/4)^(7-maxMisMatch)*choose(7-maxMisMatch) as expected
#' false positive rate combining both primer and ltr is 0.0025*0.0025=6.25E-6
#' @param reads.p DNAStringSet of reads, normally R2
#' @param primer character string of lenth 1, such as "GAAAATC"
#' @param ltrbit character string of lenth 1, such as "TCTAGCA"
#' @return DNAStringSet of reads with primer and ltr removed
#' 
trim_Ltr_side_reads <- function(reads.p, primer, ltrbit, maxMisMatch=2) {
    
    stopifnot(class(reads.p) %in% "DNAStringSet")
    stopifnot(!any(duplicated(names(reads.p))))
    stopifnot(length(primer)==1)
    stopifnot(length(ltrbit)==1)
    
    ## allows gap, and del/ins count as 1 mismatch
    submat1 <- nucleotideSubstitutionMatrix(match=1,
                                            mismatch=0,
                                            baseOnly=TRUE)
    
    ## p for primer
    ## search for primer from the beginning
    aln.p <- pairwiseAlignment(pattern=subseq(reads.p, 1, 1+nchar(primer)),
                               subject=primer,
                               substitutionMatrix=submat1,
                               gapOpening = 0,
                               gapExtension = 1,
                               type="overlap")
    aln.p.df <- PairwiseAlignmentsSingleSubject2DF(aln.p)
    
    ## l for ltrbit
    ## search for ltrbit fellowing primer
    aln.l <- pairwiseAlignment(pattern=subseq(reads.p, nchar(primer), nchar(primer)+nchar(ltrbit)+1),
                               subject=ltrbit,
                               substitutionMatrix=submat1,
                               gapOpening = 0,
                               gapExtension = 1,
                               type="overlap")
    aln.l.df <- PairwiseAlignmentsSingleSubject2DF(aln.l, shift=nchar(primer)-1)
    
    goodIdx <- (aln.p.df$score >= nchar(primer)-maxMisMatch &
                aln.l.df$score >= nchar(ltrbit)-maxMisMatch)
    
    reads.p <- subseq(reads.p[goodIdx], aln.l.df$end[goodIdx]+1)
    
    return(reads.p)
}
##trim_Ltr_side_reads(reads.p, primer, ltrbit)


#' subset and substring
#' trim primerID linker side of read, R1 in protocol
#' a primerIDlinker has N's in the middle
#' allow 3 mismatches for either part before and after Ns
#' see reasonning above
#' @param reads.l DNAStringSet of reads, normally R1
#' @param linker character string of lenth 1, such as
#'               "AGCAGGTCCGAAATTCTCGGNNNNNNNNNNNNCTCCGCTTAAGGGACT"
#' @param maxMisMatch=3
#' @return list of read.l and primerID
#' 
trim_primerIDlinker_side_reads <- function(reads.l, linker, maxMisMatch=3) {
    
    stopifnot(class(reads.l) %in% "DNAStringSet")
    stopifnot(!any(duplicated(names(reads.l))))
    stopifnot(length(linker)==1)
    
    pos.N <- unlist(gregexpr("N", linker))
    len.N <- length(pos.N)
    link1 <- substr(linker, 1, min(pos.N)-1)
    link2 <- substr(linker, max(pos.N)+1, nchar(linker))
    
    ## allows gap, and del/ins count as 1 mismatch
    submat1 <- nucleotideSubstitutionMatrix(match=1,
                                            mismatch=0,
                                            baseOnly=TRUE)
    
    ## search at the beginning for 1st part of linker
    aln.1 <- pairwiseAlignment(pattern=subseq(reads.l, 1, 2+nchar(link1)),
                               subject=link1,
                               substitutionMatrix=submat1,
                               gapOpening = 0,
                               gapExtension = 1,
                               type="overlap")
    aln.1.df <- PairwiseAlignmentsSingleSubject2DF(aln.1)
    
    ## search after 1st part of linker for the 2nd part of linker
    aln.2 <- pairwiseAlignment(pattern=subseq(reads.l, max(pos.N)-1, nchar(linker)+1),
                               subject=link2,
                               substitutionMatrix=submat1,
                               gapOpening = 0,
                               gapExtension = 1,
                               type="overlap")
    aln.2.df <- PairwiseAlignmentsSingleSubject2DF(aln.2, max(pos.N)-2)
    
    goodIdx <- (aln.1.df$score >= nchar(link1)-maxMisMatch &
                aln.2.df$score >= nchar(link2)-maxMisMatch)
    
    primerID <- subseq(reads.l[goodIdx],
                       aln.1.df$end[goodIdx]+1,
                       aln.2.df$start[goodIdx]-1)
    
    reads.l <- subseq(reads.l[goodIdx], aln.2.df$end[goodIdx]+1)
    
    stopifnot(all(names(primerID)==names(reads.l)))
    
    return(list("reads.l"=reads.l,
                "primerID"=primerID))
}
##trim_primerIDlinker_side_reads(reads.l, linker)


#' subseqing, trim off reads from where marker start to match
#' when human part of sequence is short, ltr side read will read in to 
#' linker, and linker side reads may read into ltrbit, primer, etc
#' allow 1 mismatch for linker common
#' @param reads DNAStringSet of reads
#' @param marker over reading marker
#' @return DNAStringSet of reads with linker sequences removed
#' 
trim_overreading <- function(reads, marker, maxMisMatch=3) {
    
    stopifnot(class(reads) %in% "DNAStringSet")
    stopifnot(!any(duplicated(names(reads))))
    stopifnot(length(marker)==1)
    
    
    submat1 <- nucleotideSubstitutionMatrix(match=1,
                                            mismatch=0,
                                            baseOnly=TRUE)
    
    ## allows gap, and del/ins count as 1 mismatch
    tmp <- pairwiseAlignment(pattern=reads,
                             subject=marker,
                             substitutionMatrix=submat1,
                             gapOpening = 0,
                             gapExtension = 1,
                             type="overlap")
    
    odf <- PairwiseAlignmentsSingleSubject2DF(tmp)
    
    odf$isgood <- FALSE
    ## overlap in the middle or at right
    odf$isgood <- with(odf, ifelse(mismatch<=maxMisMatch &
                                   start>1,
                                   TRUE, isgood))
    
    ## overlap at left
    odf$isgood <- with(odf, ifelse(mismatch<=maxMisMatch &
                                   start==1 &
                                   width>=nchar(marker)-1,
                                   TRUE, isgood))
    
    ## note with ovelrap alignmment, it only align with a minimum of 1/2 of the shorter one
    odf$cut <- with(odf, ifelse(isgood, odf$start-1, nchar(reads)))
    if( any(odf$cut < nchar(reads)) ) {
        odf$cut <- nchar(reads)-nchar(marker)/2
        odf$cut <- with(odf, ifelse(isgood, odf$start-1, cut))
    }
    
    reads <- subseq(reads, 1, odf$cut)
    

    ## debugging and machine learning code
    ## for linker_common, 16 bases, maxMisMatch=3 is the best
    ## mwidth <- as.integer(stringr::str_match(names(reads), "(\\d+):(\\d+)$")[,2])
    ## truecut <- ifelse( mwidth<nchar(reads), mwidth+1, nchar(reads) )
    ## ovldf <- data.frame(size=nchar(linker_common),
    ##                    width=width(pattern(tmp)),
    ##                    score=score(tmp),
    ##                    mismatch=width(pattern(tmp))-score(tmp),
    ##                    start=start(pattern(tmp)),
    ##                    end=end(pattern(tmp)),
    ##                    rlen=nchar(reads),
    ##                    mwidth=mwidth,
    ##                    cut=ifelse( mwidth<nchar(reads), mwidth+1, nchar(reads) ))
    ## ovldf$good <- ( abs(ovldf$cut-ovldf$start)<=2 )
    ## library(rpart)
    ## fit <- rpart(good ~ size + width + score + mismatch + start + end + rlen,
    ##             data=ovldf,
    ##             method="class")
}
##trim_overreading(reads.p, linker_common)
##trim_overreading(reads.l, largeLTRFrag)




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
  message("Entering ", workingDir)
  
  stats.bore <- data.frame(sample=alias)
  message("\nTrim reads with low quality bases")  
  
  ##filenames <- list(read1, read2)
  ##
  ##reads <- lapply(filenames, function(x) {
  ##  sapply(x, readFastq)
  ##})
  reads <- lapply(list(read1, read2), sapply, readFastq)
  
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
  rm(r)
  gc()
  
  stats.bore$Reads.p.afterTrim <- length(reads[[2]])
  stats.bore$Reads.l.afterTrim <- length(reads[[1]])
  print(t(stats.bore), quote=FALSE)
  
  message("\nTrim primer and ltrbit")
  
  #'.p' suffix signifies the 'primer' side of the amplicon (i.e. read2)
  #'.l' suffix indicates the 'liner' side of the amplicon (i.e. read1)
  
  reads.p <- trim_Ltr_side_reads(reads[[2]], primer, ltrbit)
  stats.bore$primed <- length(reads.p)
  stats.bore$LTRed <- length(reads.p)
  
  readslprimer <- trim_primerIDlinker_side_reads(reads[[1]], linker)
  reads.l <- readslprimer$reads.l
  primerIDs <-readslprimer$readslprimer$primerID
  
  stats.bore$linkered <- length(reads.l)
  save(primerIDs, file="primerIDData.RData")
  
  print(t(stats.bore), quote=FALSE) 
  
  ## check if reads were sequenced all the way by checking for opposite adaptor
  message("\nTrim reads.p over reading into linker")
  reads.p <- trim_overreading(reads.p, linker_common, 3)
  
  message("\nTrim reads.l over reading into ltr")
  ##reads.l <- trim_overreading(reads.l, largeLTRFrag)
  ## here we use first 20 because whole largeLTRFrag is too long
  ## with mismatch=3, the 20 bases can not be found in human genome
  reads.l <- trim_overreading(reads.l, substr(largeLTRFrag, 1, 20), 3)
  
  message("\nFilter on minimum length of ", mingDNA)
  reads.p <- subset(reads.p, width(reads.p) > mingDNA)
  reads.l <- subset(reads.l, width(reads.l) > mingDNA)
  
  stats.bore$reads.p_afterTrim <- length(reads.p)
  stats.bore$reads.l_afterTrim <- length(reads.l)
  
  
  message("\nRemove reads align to vector") 
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
  
  names(reads.p.u) <- seq_along(reads.p.u)
  names(reads.l.u) <- seq_along(reads.l.u)
  
  keys <- data.frame("R2"=match(reads.p, reads.p.u),
                     "R1"=match(reads.l, reads.l.u),
                     "names"=toload)
  
  save(keys, file="keys.RData")
  
  if(length(toload) > 0){
      ## devide reads by chunks of 30000
      chunks.p <- split(seq_along(reads.p.u), ceiling(seq_along(reads.p.u)/30000))
      for(i in c(1:length(chunks.p))){
          writeXStringSet(reads.p.u[chunks.p[[i]]],
                          file=paste0("R2-", i, ".fa"),
                          append=FALSE)
      }
      
      chunks.l <- split(seq_along(reads.l.u), ceiling(seq_along(reads.l.u)/30000))
      for(i in c(1:length(chunks.l))){    
          writeXStringSet(reads.l.u[chunks.l[[i]]],
                          file=paste0("R1-", i, ".fa"),
                          append=FALSE)
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
    mcols(algns.gr) <- algns[,c("matches", "misMatches", "qStart", "qEnd", "qSize", "tBaseInsert", "from")]
    rm(algns)
    algns.gr
  }
  
  load("keys.RData")
  
  hits.R2 <- processBLATData(read.psl(system("ls R2*.fa.psl.gz", intern=T),
                                      bestScoring=F, removeFile=F),
                             "R2")
  save(hits.R2, file="hits.R2.RData")
  
  hits.R1 <- processBLATData(read.psl(system("ls R1*.fa.psl.gz", intern=T),
                                      bestScoring=F, removeFile=F),
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
