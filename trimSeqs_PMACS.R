source('~/EAS/PMACS_scripts/pairedIntSites.R')

sampleID = as.integer(system("echo $LSB_JOBINDEX", intern=T))
  
parameters = get(load("parameters.RData"))[[sampleID]]

alias = parameters[[12]]

suppressWarnings(dir.create(alias, recursive=TRUE))

status = tryCatch(eval(as.call(append(getTrimmedSeqs, parameters[1:13]))), error=function(e){print(paste0("Caught error: ", e$message))})

save(status, file=paste0(alias, "/trimStatus.RData"))
