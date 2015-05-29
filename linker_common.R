#' Extract the second part of sequence of the linker sequence and
#' reverse complement; in case a linker has no N in the middle,
#' get the last 15 bases and reverse complement.
#' @param linkerSequence, vector of linker sequences
#' @return linkerCommon, vector of the sequence after N
linker_common <- function(linkerSequence) {
    stopifnot(require("Biostrings"))
    stopifnot(is.character(linkerSequence))
    splited <- strsplit(linkerSequence, 'N+')
    
    splited.len <- sapply(splited, length)
    stopifnot( all(splited.len == 2 | splited.len == 1))
    
    afterN_RCed <- lapply(splited, function(x) {
        ## linker with continuous N in the middle
        if ( length(x)==2 ) return(as.character(reverseComplement(DNAString(x[2])))) 
        ## linker with no N in the middle
        if ( length(x)==1 ) {
            if( nchar(x)<=15 ) stop("Linker should always longer than 15 bases")
            last15 <- substr(x, nchar(x)-15+1, nchar(x))
            return(as.character(reverseComplement(DNAString(last15)))) 
        }
        
        ## seperated N, this situation should not happen
        stop("Unknow situation:", paste(x, collapse=" "))
    })
    afterN_RCed <- unlist(afterN_RCed)
    
    stopifnot(length(afterN_RCed)==length(linkerSequence))
    stopifnot(nchar(afterN_RCed)>=15)
    
    return(afterN_RCed)
}





