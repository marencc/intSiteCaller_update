# do alignments
library("rtracklayer")
library("GenomicRanges")

blatPort  = as.integer(system("echo $LSB_JOBINDEX", intern=T))

toAlign = system("ls */*.fa", intern=T)

alignFile = toAlign[blatPort-get(load("bushmanBlatStartPort.RData"))+1] #blatPort is actually file offset by first blat port number and R is offset by 1

indexPath = get(load("indexPath.RData"))

indexSeqInfo = seqinfo(TwoBitFile(indexPath))
isCircular(indexSeqInfo) = rep(F, length(indexSeqInfo))
genome(indexSeqInfo) = strsplit(strsplit(indexPath, ".2bit")[[1]], "/")[[1]][length(strsplit(indexPath,"/")[[1]])]

save(indexSeqInfo, file="indexSeqInfo.RData")

system(paste0("/home/aubreyba/EAS/PMACS_scripts/BLATsamples.sh ", alignFile, " ", blatPort, " ", indexPath))