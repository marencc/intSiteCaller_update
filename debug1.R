libs <- c("plyr", "dplyr", "BiocParallel", "Biostrings", "GenomicAlignments" ,"hiAnnotator" ,"sonicLength", "GenomicRanges", "BiocGenerics", "ShortRead", "GenomicRanges", "igraph")
null <- suppressMessages(sapply(libs, library, character.only=TRUE))

source("~/intSiteCaller/hiReadsProcessor.R")


processPsl <- function(algns, from){
    stopifnot(from == "R1" | from == "R2")
    
    algns$from <- from
    algns <- merge(algns, keys, by.x="qName", by.y=from)
    
    algns <- subset(algns, strand=="+" | strand=="-")
    
    needed <- c("tName", "tStart", "tEnd", "strand", "matches", "repMatches", "misMatches", "qStart", "qEnd", "qSize", "tBaseInsert", "from", "uPID")
    
    algns <- algns[, needed]
    
    algns <- dplyr::distinct(algns) %>%
        dplyr::filter( 100*(matches+repMatches)/qSize >minPercentIdentity &
                          qStart<=maxAlignStart ) %>%
            dplyr::select(uPID, tName, strand, tStart, tEnd, from)
    
    return(algns)
}

load("keys.RData")
keys$uPID <- paste(keys$R2, keys$R1, sep=":")

psl.R1 <- read.psl(list.files(".", "R1.*.fa.psl.gz"),
                   bestScoring=F, removeFile=F)
psl.R1 <- processPsl(psl.R1, "R1")
strand <- c("+", "-")
psl.R1$Cstrand <- strand[3-match(psl.R1$strand, strand)]
psl.R1$strand <- psl.R1$Cstrand


psl.R2 <- read.psl(list.files(".", "R2.*.fa.psl.gz"),
                   bestScoring=F,
                   removeFile=F)
psl.R2 <- processPsl(psl.R2, "R2")

psl.Pair <- merge(psl.R2, psl.R1,
                  by=c("uPID", "tName", "strand"),
                  suffixes = c(".R2",".R1"))

psl.Pair <- dplyr::mutate(psl.Pair,
                          position=ifelse(strand=="+", tStart.R2, tEnd.R2),
                          breakpoint=ifelse(strand=="+", tEnd.R1, tStart.R1),
                          chr=tName) %>%
    dplyr::filter( (strand=="+" & breakpoint-position>0 & breakpoint-position<2500) |
                      (strand=="-" & breakpoint-position<0 & breakpoint-position>-2500) ) %>%
    dplyr::select(uPID, chr, strand, position, breakpoint) %>%
        dplyr::arrange(uPID, chr, position, breakpoint) %>%
            dplyr::distinct()

keysCount <- dplyr::group_by(keys, uPID) %>% dplyr::summarize(count=n())

psl.Pair <- merge(psl.Pair, keysCount, by="uPID")



trim_primerIDlinker_side_reads <- function(reads.l, linker, maxMisMatch=3) {
    
    ## note how maxMisMatch was set explained in trim_ltr_side_reads
    
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
    
    return(list("reads.l"=reads.l, "primerID"=primerID))
}


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
    
    ##runSeq <- sapply(1:10000, function(i)
    ##                 paste(sample(c("A","C","G","T"), 8, replace=TRUE),
    ##                       collapse=""))
    ##runSeq.p <- pairwiseAlignment(pattern=runSeq,
    ##                              subject=primer,
    ##                              substitutionMatrix=submat1,
    ##                              gapOpening = 0,
    ##                              gapExtension = 1,
    ##                              type="overlap")
    ##runSeq.p.df <- PairwiseAlignmentsSingleSubject2DF(runSeq.p)
    ##table(runSeq.p.df$score)
    ##  1    2    3    4    5    6    7
    ##211 3338 4330 1751  344   25    1
    ##false positive rate 0.0025 and thus maxMisMatch=2, (1/4)^maxMisMatch*choose(7-maxMisMatch)
    ##false positive rate combining both primer and ltr is 0.0025*0.0025=6.25E-6
    
    goodIdx <- (aln.p.df$score >= nchar(primer)-maxMisMatch &
                aln.l.df$score >= nchar(ltrbit)-maxMisMatch)
    
    reads.p <- subseq(reads.p[goodIdx], aln.l.df$end[goodIdx]+1)
    
    return(reads.p)
}


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


trim_overreading <- function(reads, marker, maxMisMatch=3) {
    reads=reads.p
    marker=linker_common
    maxMisMatch=3
    reads=reads.l
    marker=substr(largeLTRFrag,1,20)
    maxMisMatch=3
    
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


