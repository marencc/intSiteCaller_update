#' Extract the second part of sequence of the linker sequence and
#' reverse complement; in case a linker has no N in the middle,
#' get the last min_len bases and reverse complement.
#' @param linkerSequence, vector of linker sequences
#' @param min_len length of linker for the case of no N
#'    (with N minumum linker accepted as valid)
#' @return linkerCommon, vector of the sequence after N
linker_common <- function(linkerSequence, min_len=15) {
    stopifnot(require("Biostrings")) # should we test for libraries presence?
    stopifnot(is.character(linkerSequence)) # should we check types of args?
    splited <- strsplit(linkerSequence, 'N+')
    
    splited.len <- sapply(splited, length)
    stopifnot( all(splited.len == 2 | splited.len == 1))
    
    common_seq <- sapply(splited, function(x) {
        ## linker with continuous N in the middle
        if ( length(x)==2 ) 
            return(.get_rev_seq(x[2]))
        ## linker with no N in the middle
        if ( length(x)==1 ) {
            if( nchar(x) <= min_len ) stop(
                paste("Linker should always longer than", min_len, " bases"))
            last_seq <- substr(x, nchar(x)-min_len+1, nchar(x))
            return(.get_rev_seq(last_seq))
        }
    })
    stopifnot(nchar(common_seq)>=min_len) # small linkers are not used in prep
    return(common_seq)
}

.get_rev_seq <- function(seq) {
    as.character(reverseComplement(DNAString(seq)))
}
