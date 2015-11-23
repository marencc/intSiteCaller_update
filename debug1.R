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


