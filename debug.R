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

