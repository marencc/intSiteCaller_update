#### debugging vector trimming parameters ####
## It trimmed 10% of perfect human reads ##

library(ggplot2)
library(Biostrings)
library(stringr)
library(gridExtra)

load("../completeMetadata.RData")

hits.p <- get(load("hits.v.p.RData"))
hits.l <- get(load("hits.v.l.RData"))

Vector <- readDNAStringSet("../vector_sim.fa")

primer <- paste0(completeMetadata$primer, completeMetadata$ltrBit)

primer.loci <- str_locate_all(Vector, primer)

aln.p <- qplot(hits.p$tStart, xlim=c(0,10000)) +
    geom_vline(xintercept =primer.loci[[1]][,"start"])

aln.l <- qplot(hits.l$tStart, xlim=c(0,10000)) +
    geom_vline(xintercept =primer.loci[[1]][,"start"])

## alignments of trimmed reads on vector
grid.arrange(aln.p, aln.l)


pOI <- c(hits.p$matches/hits.p$qSize, hits.l$matches/hits.l$qSize)
rLen <- c(hits.p$qSize, hits.l$qSize)
qplot(pOI)
qplot(rLen, pOI, xlab="Read length", ylab="Percent of identity")



hits.v.p <- get(load("hits.v.p.RData"))
hits.v.l <- get(load("hits.v.l.RData"))
load("reads.p.RData")
load("reads.l.RData")

#' find reads originating from vector
#' @param vectorSeq vector sequence fasta file
#' @param reads.p DNAStringSet, reads on primer side
#' @param reads.l DNAStringSet, reads on linker side
#' @return character, qNames for the vector reads
findVectorReads <- function(vectorSeq, reads.p, reads.l) {
    
    Vector <- readDNAStringSet(file.path("..", vectorSeq))
    
    blatParameters <- c(minIdentity=88, minScore=30, stepSize=3, 
                        tileSize=8, repMatch=112312, dots=1000, 
                        q="dna", t="dna", out="psl")
    
    
    hits.v.p <- try(read.psl(blatSeqs(query=reads.p, subject=Vector,     
                                      blatParameters=blatParameters, parallel=F),
                             bestScoring=F) )
    if( class(hits.v.p) == "try-error" ) hits.v.p <- data.frame()
    
    hits.v.l <- try(read.psl(blatSeqs(query=reads.l, subject=Vector, 
                                      blatParameters=blatParameters, parallel=F),
                             bestScoring=F) )
    if( class(hits.v.l) == "try-error" ) hits.v.l <- data.frame()
    
    hits.v.p <- dplyr::filter(hits.v.p, matches>0.85*qSize & strand=="+" & qStart<5) 
    hits.v.l <- dplyr::filter(hits.v.l, matches>0.85*qSize & strand=="-") 
    
    hits.v <- try(merge(hits.v.p[, c("qName", "tStart")],
                        hits.v.l[, c("qName", "tStart")],
                        by="qName")
                 ,silent = TRUE)
    if( class(hits.v) == "try-error" ) hits.v <- data.frame()
    
    hits.v <- dplyr::filter(hits.v, tStart.y>=tStart.x & tStart.y<tStart.x+2000)
    
    message(nrow(hits.v), " vector sequences found")
    return(hits.v$qName)
}
## vqName <- findVectorReads(vectorSeq, reads.p, reads.l)

vqName <- findVectorReads(vectorSeq, reads.p, reads.l)
reads.p <- reads.p[!names(reads.p) %in% vqName]
reads.l <- reads.l[!names(reads.l) %in% vqName]


try(merge(hits.v.p[, c("qName", "tStart")],
          hits.v.l[, c("qName", "tStart")],
          by="qName1")
    ,silent = TRUE)



findAndRemoveVector.eric <- function(reads, Vector, blatParameters, minLength=10){
    
    
    suppressWarnings(file.remove("hits.v.RData"))
    save(hits.v, file="hits.v.RData")    
    #collapse instances where a single read has multiple vector alignments
    hits.v <- reduce(GRanges(seqnames=hits.v$qName, IRanges(hits.v$qStart,
                                                            hits.v$qEnd)),
                     min.gapwidth=1200)
    names(hits.v) <- as.character(seqnames(hits.v))
    
    hits.v <- hits.v[start(hits.v)<=5 & width(hits.v)>minLength]
    
    reads[!names(reads) %in% names(hits.v)]
    
  }
  
  save(reads.p, file="reads.p.RData")
  tryCatch(reads.p <- findAndRemoveVector.eric(reads.p, Vector,
                                          blatParameters=blatParameters),
           error=function(e){print(paste0("Caught ERROR in intSiteLogic::findAndRemoveVector ",
               e$message))})
  suppressWarnings(file.rename("hits.v.RData", "hits.v.p.RData"))
  
  save(reads.l, file="reads.l.RData")
  tryCatch(reads.l <- findAndRemoveVector.eric(reads.l, Vector,
                                          blatParameters=blatParameters),
           error=function(e){print(paste0("Caught ERROR in intSiteLogic::findAndRemoveVector ",
               e$message))})
  suppressWarnings(file.rename("hits.v.RData", "hits.v.l.RData"))

