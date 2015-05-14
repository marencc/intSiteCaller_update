### SET RUN PARAMETERS HERE ###
bushmanJobID <- "intSiteValidation" #allows simultaneous processing of datasets - make sure to use unique BLAT ports!
blatStartPort <- 5560 #this can get a bit weird since spawning a bunch of blat threads could result in conflicts with other processes
codeDir <- "/home/aubreyba/EAS/intSiteCaller" #Aubrey's PMACS account - change for your application!
cleanup <- TRUE
### END RUN PARAMETERS ###

save(bushmanJobID, file=paste0(getwd(), "/bushmanJobID.RData"))
save(blatStartPort, file=paste0(getwd(), "/bushmanBlatStartPort.RData"))
save(codeDir, file=paste0(getwd(), "/codeDir.RData"))
save(cleanup, file=paste0(getwd(), "/cleanup.RData"))

sampleInfo <- read.csv("sampleInfo.csv", stringsAsFactors=F)
processingParams <- read.csv("processingParams.csv", stringsAsFactors=F)

#confirm that metadata is presented as we expect
stopifnot(nrow(sampleInfo) == nrow(processingParams))
stopifnot(!(is.null(sampleInfo$alias) | is.null(processingParams$alias)))
stopifnot(all(sampleInfo$alias %in% processingParams$alias))

completeMetadata <- merge(sampleInfo, processingParams, "alias")

completeMetadata$gender[with(completeMetadata, gender==F)] <- "F"
completeMetadata$gender[with(completeMetadata, gender=="m")] <- "M"

completeMetadata$read1 <- paste0(getwd(), "/Data/demultiplexedReps/", completeMetadata$alias, "_R1.fastq.gz")
completeMetadata$read2 <- paste0(getwd(), "/Data/demultiplexedReps/", completeMetadata$alias, "_R2.fastq.gz")

stopifnot(all(c("qualityThreshold", "badQualityBases", "qualitySlidingWindow",
                "primer", "ltrBit", "largeLTRFrag", "linkerSequence", "linkerCommon",
                "mingDNA", "read1", "read2", "alias", "vectorSeq", "minPctIdent",
                "maxAlignStart", "maxFragLength", "gender") %in% names(completeMetadata)))

save(completeMetadata, file="completeMetadata.RData")

suppressWarnings(dir.create("logs"))

source(paste0(codeDir, "/programFlow.R")) #to get the bsub function - doing this here so codeDir is available

#error-correct barcodes - kicks off subsequent steps
bsub(queue="plus",
     jobName=paste0("BushmanErrorCorrect_", bushmanJobID),
     logFile="logs/errorCorrectOutput.txt",
     command=paste0("Rscript -e \"source('", codeDir, "/programFlow.R'); errorCorrectBC();\"")
)
