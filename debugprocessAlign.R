#    This source code file is a component of the larger INSPIIRED genomic analysis software package.
#    Copyright (C) 2016 Frederic Bushman
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

codeDir <- get(load("codeDir.RData"))
source(file.path(codeDir, "intSiteLogic.R"))

sampleID <- 5
sampleID <- 81
message(sampleID)

completeMetadata <- get(load("completeMetadata.RData"))[sampleID,]
print(t(completeMetadata), quote=FALSE)  

workingDir=completeMetadata$alias
minPercentIdentity=completeMetadata$minPctIdent
maxAlignStart=completeMetadata$maxAlignStart
maxLength=2500
refGenome="hg18"
refGenome="mm9"

## processAlignments <- function(workingDir, minPercentIdentity, maxAlignStart, maxLength, refGenome)

codeDir <- get(load("codeDir.RData"))
source(file.path(codeDir, "programFlow.R"))#for get_reference_genome function

setwd(workingDir)
message("Entering ", workingDir)

load("keys.RData")
keys$uPID <- paste(keys$R2, keys$R1, sep=":")
keys <- dplyr::group_by(keys, uPID) %>% dplyr::mutate(count=n()) 
ukeys <- (dplyr::group_by(keys, uPID) %>%
          dplyr::mutate(count=n(), uid=1:n()) %>%
          dplyr::top_n(n=1, uid) )


#' read psl gz files, assuming psl gz files don't have column header
#' @param pslFile character vector of file name(s)
#' @param toNull  character vector of column names to get rid of
#' @return data.frame, data.table of the psl table
#' @example 
readpsl <- function(pslFile, toNull=NULL) {
    cols <- c("matches", "misMatches", "repMatches", "nCount", "qNumInsert",
              "qBaseInsert", "tNumInsert", "tBaseInsert", "strand", "qName",
              "qSize", "qStart", "qEnd", "tName", "tSize", "tStart", "tEnd",
              "blockCount", "blockSizes", "qStarts", "tStarts")
    cols.class <- c(rep("numeric",8), rep("character",2), rep("numeric",3),
                    "character", rep("numeric",4), rep("character",3))
    
    psl <- lapply(pslFile, function(f) {
        message("Reading ",f)
        data.table::fread( paste("zcat", f), sep="\t" )
    })
    psl <- data.table::rbindlist(psl)
    colnames(psl) <- cols
    
    if(length(toNull)>0) psl[, toNull] <- NULL
    
    return(as.data.frame(psl))
}
##toNull <- c("blockCount", "blockSizes", "qStarts", "tStarts",
##            "nCount", "qNumInsert", "qBaseInsert", "tNumInsert",
##            "tBaseInsert")
##psl <- readpsl(pslFile, toNull=toNull)


#' Attach uPID, unique read pair id, unique in terms of base pairs;
#' Attach count, number of exact PCR duplicates;
#' Percent if identity filter;
#' Select only useful columns.
#' For efficiency, we will only load one of the PCR duplicates.
#' @param psl psl data frame
#' @param from R1 or R2
#' @param keys mapping between qNames in psl to real qNames
#' @param minPercentIdentity threshold
#' @note, there could be base errors in PCR duplicates,
#' if two reads align to the same location, they are usually considered
#' as PCR duplicates as well. This is a common practice.
processPsl <- function(psl, from, keys, minPercentIdentity=95){
    
    neededCols <- c("qName", "tName", "strand", "tStart", "tEnd", "qStart")
    stopifnot(neededCols %in% colnames(psl))
    stopifnot(from == "R1" | from == "R2")
    stopifnot(c("R2", "R1", "names", "uPID", "count") %in% colnames(keys))
    stopifnot(minPercentIdentity<=100 & minPercentIdentity>70)
    
    psl <- (dplyr::mutate(psl,
                          from=from,
                          POI=as.integer(100*(matches+repMatches)/qSize)) %>%
            dplyr::filter(strand=="+" | strand=="-") %>%
            dplyr::filter(POI>=minPercentIdentity) %>%
            dplyr::inner_join(keys, by = c("qName"=from)) %>%
            dplyr::select(uPID, count, tName, strand, tStart, tEnd,
                          from, POI, qStart) ) 
    
    return(psl)
}

toNull <- c("blockCount", "blockSizes", "qStarts", "tStarts",
            "nCount", "qNumInsert", "qBaseInsert", "tNumInsert",
            "tBaseInsert")

psl.R1 <- readpsl(list.files(".", "R1.*.fa.psl.gz"), toNull=toNull)
psl.R1 <- processPsl(psl.R1, from="R1", keys=ukeys, minPercentIdentity)
## complement R1 strand
strand <- c("+", "-")
psl.R1$Cstrand <- strand[3-match(psl.R1$strand, strand)]
psl.R1$strand <- psl.R1$Cstrand


psl.R2 <- readpsl(list.files(".", "R2.*.fa.psl.gz"), toNull=toNull)
psl.R2 <- processPsl(psl.R2, from="R2", keys=ukeys, minPercentIdentity)
psl.R2 <- dplyr::filter(psl.R2, qStart<=maxAlignStart)

## merge by unique pair id, chomosome, and strand
## note, strands of read1 have been complemented
psl.Pair <- dplyr::inner_join(psl.R2, psl.R1,
                              by = c("uPID"="uPID",
                                  "tName"="tName",
                                  "strand"="strand"))

colnames(psl.Pair) <- sub(".x$", ".R2", colnames(psl.Pair))
colnames(psl.Pair) <- sub(".y$", ".R1", colnames(psl.Pair))

rm(psl.R1, psl.R2)
gc()


psl.Pair <- dplyr::mutate(psl.Pair,
                          position=ifelse(strand=="+", tStart.R2, tEnd.R2),
                          breakpoint=ifelse(strand=="+", tEnd.R1, tStart.R1))

psl.Pair <- dplyr::filter(psl.Pair, abs(breakpoint-position)<maxLength)

psl.Pair <- dplyr::select(psl.Pair,
                          uPID,
                          chr=tName, strand, position, breakpoint,
                          count=count.R2)

psl.Pair <- (dplyr::group_by(psl.Pair, uPID) %>%
             dplyr::mutate(hits=n()) %>%
             dplyr::ungroup() )

psl.Pair.multi <- (dplyr::filter(psl.Pair, hits>1) %>%
                   dplyr::mutate(uPIDi=as.integer(as.factor(uPID))))

psl.Pair.multi.gr <- makeGRangesFromDataFrame(psl.Pair.multi,
                                              seqnames.field="chr",
                                              start.field="position",
                                              end.field="position",
                                              strand.field="strand",
                                              keep.extra.columns=TRUE)

psl.Pair.multi.gr.red <- reduce(psl.Pair.multi.gr,
                                min.gapwidth=5,
                                with.revmap=TRUE,
                                ignore.strand=FALSE)

pair.revmap <- lapply(psl.Pair.multi.gr.red$revmap,
                      function(idx) psl.Pair.multi$uPIDi[idx])

cid <- rep(NA, max(psl.Pair.multi$uPIDi))
for(i in 1:length(pair.revmap)) {
    ##message(i)
    posidx <- pair.revmap[[i]]
    if( all(is.na(cid[posidx])) ) {
        cid[posidx] <- i
    } else {
        precid <- unique(cid[posidx])
        precid <- precid[!is.na(precid)]
        mincid <- min(precid)
        cid[cid %in% precid] <- mincid
        cid[posidx] <- mincid
    }
}



