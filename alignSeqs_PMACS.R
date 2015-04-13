# do alignments

blatPort  = as.integer(system("echo $LSB_JOBINDEX", intern=T))

toAlign = system("ls */*.fa", intern=T)

alignFile = toAlign[blatPort-get(load("bushmanBlatStartPort.RData"))+1] #blatPort is actually file offset by first blat port number and R is offset by 1

system(paste0("/home/aubreyba/EAS/PMACS_scripts/BLATsamples.sh ", alignFile, " ", blatPort))