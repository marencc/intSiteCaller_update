PairwiseAlignmentsSingleSubject2DF <- function(PA) {
    stopifnot("PairwiseAlignmentsSingleSubject"  %in% class(PA))
    
    return(data.frame(
        width=width(pattern(PA)),
        score=score(PA),
        mismatch=width(pattern(PA))-score(PA),
        start=start(pattern(PA)),
        end=end(pattern(PA))
        ))
}


trim_overreading <- function(reads, marker, maxMisMatch=3) {
    reads=reads.p
    marker=linker_common
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
    ## aligned to middle or right
    odf$isgood <- with(odf, ifelse(mismatch<=maxMisMatch &
                                   start>1,
                                   TRUE, isgood))
    
    ## aligned to left
    odf$isgood <- with(odf, ifelse(mismatch<=maxMisMatch &
                                   start==1 &
                                   width>=nchar(marker)-1,
                                   TRUE, isgood))
    
    
    mwidth <- as.integer(stringr::str_match(names(reads), "(\\d+):(\\d+)$")[,2])
    
    ovldf <- data.frame(size=nchar(linker_common),
                        width=width(pattern(tmp)),
                        score=score(tmp),
                        mismatch=width(pattern(tmp))-score(tmp),
                        start=start(pattern(tmp)),
                        end=end(pattern(tmp)),
                        rlen=nchar(reads),
                        mwidth=mwidth,
                        cut=ifelse( mwidth<nchar(reads), mwidth+1, nchar(reads) ))
    
    ovldf$good <- ( abs(ovldf$cut-ovldf$start)<=2 )
    
    library(rpart)
    fit <- rpart(good ~ size + width + score + mismatch + start + end + rlen,
                 data=ovldf,
                 method="class")
}


