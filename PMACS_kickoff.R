print("started analysis")

source("~/EAS/PMACS_scripts/programFlow.R") #to get the bsub function

metadata = read.csv(paste0(getwd(), "/sampleInfo.csv"), stringsAsFactors=F)
processingParams = read.csv(paste0(getwd(), "/processingParams.csv"), stringsAsFactors=F)

stopifnot(nrow(metadata) == nrow(processingParams))

metadata = merge(metadata, processingParams, "alias")

metadata$gender[with(metadata, gender==F)] = "F"
metadata$gender[with(metadata, gender=="m")] = "M"

metadata$read1 = paste0(getwd(), "/Data/demultiplexedReps/&", metadata$alias, "_S0_L001_R1_001.fastq.gz")
metadata$read2 = paste0(getwd(), "/Data/demultiplexedReps/&", metadata$alias, "_S0_L001_R2_001.fastq.gz")

stopifnot(all(c("qualityThreshold", "badQualityBases", "qualitySlidingWindow",
                "primer", "ltrBit", "largeLTRFrag", "linkerSequence", "linkerCommon",
                "mingDNA", "read1", "read2", "alias", "vectorSeq", "minPctIdent",
                "maxAlignStart", "maxFragLength", "gender") %in% names(metadata)))

parameters = list()

for(i in c(1:nrow(metadata))){ #probably a nicer way to split a data frame into lists
  parameters = append(parameters, list(as.list(metadata[i,])))
}

save(parameters, file="parameters.RData")

suppressWarnings(dir.create("logs"))

bushmanJobID = "intSiteValidation" #allows simultaneous processing of datasets - make sure to use unique BLAT ports!
blatStartPort = 5560 #this can get a bit weird since spawning a bunch of blat threads could result in conflicts with other processes

save(bushmanJobID, file=paste0(getwd(), "/bushmanJobID.RData"))
save(blatStartPort, file=paste0(getwd(), "/bushmanBlatStartPort.RData"))

cleanup = TRUE
save(cleanup, file=paste0(getwd(), "/cleanup.RData"))

#error-correct barcodes - kicks off subsequent steps
bsub(queue="plus",
     jobName=paste0("BushmanErrorCorrect_", bushmanJobID),
     logFile="logs/errorCorrectOutput.txt",
     command=paste0("Rscript -e \"source('~/EAS/PMACS_scripts/programFlow.R'); errorCorrectBC();\"")
)

