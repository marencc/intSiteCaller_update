## start in /home/yinghua/run20151024/

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


processPsl <- function(algns, from, keys){
    
    stopifnot(from == "R1" | from == "R2")
    stopifnot(c("R2", "R1", "names", "uPID", "count") %in% colnames(keys))
    
    algns$from <- from
    algns <- merge(algns, keys, by.x="qName", by.y=from)
    
    algns <- (dplyr::mutate(algns,  POI=100*(matches+repMatches)/qSize) %>%
              dplyr::filter(strand=="+" | strand=="-") %>%
              dplyr::select(uPID, count, tName, strand, tStart, tEnd, from, POI, qStart))
    
    return(algns)
}

psl.R1 <- read.psl(list.files(".", "R1.*.fa.psl.gz"),
                   bestScoring=F, removeFile=F)
psl.R1 <- processPsl(psl.R1, from="R1", keys=ukeys)
psl.R1 <- dplyr::filter(psl.R1, POI>minPercentIdentity)

## complement R1 strand
strand <- c("+", "-")
psl.R1$Cstrand <- strand[3-match(psl.R1$strand, strand)]
psl.R1$strand <- psl.R1$Cstrand



psl.R2 <- read.psl(list.files(".", "R2.*.fa.psl.gz"),
                   bestScoring=F, removeFile=F)
psl.R2 <- processPsl(psl.R2, from="R2", keys=ukeys)
psl.R2 <- dplyr::filter(psl.R2, POI>minPercentIdentity &
                        qStart<=maxAlignStart)

psl.Pair <- merge(psl.R2,
                  psl.R1,
                  by=c("uPID", "tName", "strand"),
                  suffixes = c(".R2",".R1"))

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
             dplyr::mutate(hits=n()))

psl.Pair.multi <- (dplyr::filter(psl.Pair, hits>1) %>%
                   dplyr::ungroup() %>% 
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

revmap <- psl.Pair.multi.gr.red$revmap
pair.revmap <- lapply(revmap, function(idx) psl.Pair.multi$uPIDi[idx])

cid <- rep(NA, max(psl.Pair.multi$uPIDi))
for(i in 1:length(pair.revmap)) {
    message(i)
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


