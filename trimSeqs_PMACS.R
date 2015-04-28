source('~/EAS/PMACS_scripts/pairedIntSites.R')

sampleID = as.integer(system("echo $LSB_JOBINDEX", intern=T))
  
parameters = get(load("parameters.RData"))[[sampleID]]

alias = parameters[["alias"]]

suppressWarnings(dir.create(alias, recursive=TRUE))

status = tryCatch(eval(as.call(append(getTrimmedSeqs, unname(parameters[c("qualityThreshold", "badQualityBases",
                                                                          "qualitySlidingWindow", "primer", "ltrBit",
                                                                          "largeLTRFrag", "linkerSequence", "linkerCommon",
                                                                          "mingDNA", "read1", "read2", "alias", "vectorSeq")])))), 
                  error=function(e){print(paste0("Caught error: ", e$message))})

save(status, file="trimStatus.RData") #working directory is changed while executing getTrimmedSeqs
