# do alignments
library("rtracklayer")
library("GenomicRanges")

blatPort  = as.integer(system("echo $LSB_JOBINDEX", intern=T))

toAlign = system("ls */*.fa", intern=T)

alignFile = toAlign[blatPort-get(load("bushmanBlatStartPort.RData"))+1] #blatPort is actually file offset by first blat port number and R is offset by 1

alias = strsplit(alignFile, "/")[[1]][1]

metadata = read.csv("processingParams.csv")

#we'll still hardcode the index directory path since it would be a pain to load
#the whole genome from an R package into a .2bit file for BLAT
indexPath = paste0("/home/aubreyba/genomeIndices/", metadata[metadata$alias==alias,"refGenome"], ".2bit")

system(paste0("/home/aubreyba/EAS/PMACS_scripts/BLATsamples.sh ", alignFile, " ", blatPort, " ", indexPath))