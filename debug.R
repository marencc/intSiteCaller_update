#### debugging vector trimming parameters ####
## It trimmed 10% of perfect human reads ##

library(dplyr)
library(ggplot2)
library(Biostrings)
library(stringr)
library(gridExtra)

load("completeMetadata.RData")

hitspFile <- list.files(".", pattern="hits.v.p.RData", recursive=TRUE)
hitslFile <- list.files(".", pattern="hits.v.l.RData", recursive=TRUE)

hits.p  <- do.call(rbind, lapply(hitspFile, function(f) get(load(f)) ))
hits.l  <- do.call(rbind, lapply(hitslFile, function(f) get(load(f)) ))

hits.p$POI <- hits.p$matches/hits.p$qSize
hits.l$POI <- hits.l$matches/hits.l$qSize

qplot(hits.p$POI)
qplot(hits.p$qStart, hits.p$POI, geom='bin2d') +
    scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red')

qplot(hits.l$POI)
qplot(hits.l$qStart, hits.l$POI, geom='bin2d') +
    scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red')


WAS <- read.csv("region.csv", header=TRUE)

hits.p$ID <- sub("^GTSP.*%", "", hits.p$qName)
hits.l$ID <- sub("^GTSP.*%", "", hits.l$qName)

hits.p.WAS <- hits.p[hits.p$ID %in% WAS$ID, ]
hits.l.WAS <- hits.l[hits.l$ID %in% WAS$ID, ]

needed <- c("matches", "misMatches", "strand", "qSize", "qStart",
            "tStart", "POI", "ID")

hits.p.WAS <- hits.p.WAS[, needed]
hits.l.WAS <- hits.l.WAS[, needed]

hits.p.WAS <- hits.p.WAS %>% group_by(ID) %>% mutate(value=POI) %>%
    top_n(1) %>% mutate(value=NULL)
hits.l.WAS <- hits.l.WAS %>% group_by(ID) %>% mutate(value=POI) %>%
    top_n(1) %>% mutate(value=NULL)


Vector <- readDNAStringSet("vector_WasLenti.fa")
primer <- paste0(completeMetadata$primer, completeMetadata$ltrBit)[1]
primer.loci <- str_locate_all(Vector, primer)

aln.p <- qplot(hits.p.WAS$tStart) + xlim(c(1, 10000)) +
    geom_vline(xintercept =primer.loci[[1]][,"start"])
aln.l <- qplot(hits.l.WAS$tStart) +  xlim(c(1, 10000)) +
    geom_vline(xintercept =primer.loci[[1]][,"start"])
grid.arrange(aln.p, aln.l)

grid.arrange(qplot(hits.p.WAS$POI), qplot(hits.l.WAS$POI))


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

Vector <- readDNAStringSet("../p746vector.fasta")

load("../completeMetadata.RData")
primerltr <- unique(paste0(completeMetadata$primer, completeMetadata$ltrBit))

stringr::str_locate_all(Vector, primerltr)


path="~/Nobles/vector/run20150929"
get_stat <- function(path) {
    sampleInfo <- read.table(file.path(path, "sampleInfo.tsv"), header=TRUE)
    
    stats.file <- list.files(path, pattern="^stats.RData$", recursive=TRUE, full.names=TRUE)
    
    stats <- plyr:::rbind.fill( lapply(stats.file, function(x) get(load(x))) )
    
    stats$sample <- as.character(stats$sample)
    rownames(stats) <- NULL
    
    stats <- merge(data.frame(sample=as.character(sampleInfo$alias)),
                   plyr:::rbind.fill(lapply(stats.file, function(x) get(load(x)))),
                   by="sample",
                   all.x=TRUE)
stats$sample <- as.character(stats$sample)
    stats$gtsp <- sub("-\\d+$", "", stats$sample)
    stats$Replicate <- sub("GTSP\\d+-", "", stats$sample)
    rownames(stats) <- NULL
    
    needed <- c("sample", "reads.p_afterTrim", "reads.l_afterTrim", "reads.p_afterVTrim", "reads.l_afterVTrim", "curated", "totalSonicLengths", "numUniqueSites", "totalEvents")
    return( stats[, needed] )
}

stat.new <- get_stat("~/Nobles/vector/run20150929")
stat.old <- get_stat("~/Nobles/run20150929")


stat.new <- get_stat("~/Frances/vector/run20150903")
stat.old <- get_stat("~/Frances/run20150903")


stat.new <- get_stat("~/Frances/vector/run20141104")
stat.old <- get_stat("~/Frances/run20141104")


stat.comp <- merge(stat.old, stat.new, by="sample")

toShow <- c("sample", "curated.x", "curated.y", "totalSonicLengths.x", "totalSonicLengths.y", "totalEvents.x", "totalEvents.y")

write.table(stat.comp[, toShow], row.names=F, quote=FALSE, sep="\t")
