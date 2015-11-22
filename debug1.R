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


