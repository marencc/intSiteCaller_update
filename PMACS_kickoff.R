print("started analysis")

metadata = read.csv("~/EAS/intSiteValidation/sampleInfo.csv", stringsAsFactors=F)

metadata$read1 = paste0(getwd(), "/Data/demultiplexedReps/&", metadata$alias, "_S0_L001_R1_001.fastq.gz")
metadata$read2 = paste0(getwd(), "/Data/demultiplexedReps/&", metadata$alias, "_S0_L001_R2_001.fastq.gz")

metadata = metadata[,c("qualityThreshold", "badQualityBases", "qualitySlidingWindow", "primer", "ltrBit", "largeLTRFrag", "linkerSequence", "linkerCommon", "mingDNA", "read1", "read2", "alias", "vectorSeq", "minPctIdent", "maxAlignStart", "maxFragLength")]

names(metadata) = NULL

parameters = list()

for(i in c(1:nrow(metadata))){ #probably a nicer way to split a data frame into lists
  parameters = append(parameters, list(as.list(metadata[i,])))
}

save(parameters, file="parameters.RData")

suppressWarnings(dir.create("logs"))

bushmanJobID = "intSiteValidation" #allows simultaneous processing of datasets - make sure to use unique BLAT ports!
blatStartPort = 5560 #this can get a bit weird since spawning a bunch of blat threads could result in conflicts with other processes

save(bushmanJobID, file="~/EAS/intSiteValidation/bushmanJobID.RData")
save(blatStartPort, file="~/EAS/intSiteValidation/bushmanBlatStartPort.RData")

#demultiplex seqs
system(paste0('bsub -n1 -q plus -J "BushmanErrorCorrect_', bushmanJobID, '" -o logs/errorCorrectOutput.txt Rscript ~/EAS/PMACS_scripts/errorCorrectBC_PMACS.R'))
