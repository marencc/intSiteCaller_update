#' extract the second part of sequence of the linker sequence
#' @param linkerSequence, vector of linker sequences
#' @return linkerCommon, vector of the sequence after N
linker_common <- function(linkerSequence) {
    stopifnot(require("Biostrings"))
    stopifnot(is.character(linkerSequence))
    splited <- strsplit(linkerSequence, 'N+')
    stopifnot(sapply(splited, length) == 2)
    afterN <- sapply(splited, "[[", 2)
    stopifnot(length(afterN)==length(linkerSequence))
    
    afterN_RCed <- unlist(lapply(afterN, function(x) as.character(reverseComplement(DNAString(x)))))
    stopifnot(length(afterN_RCed)==length(linkerSequence))
    
    return(afterN_RCed)
}
