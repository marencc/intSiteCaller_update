source('~/EAS/PMACS_scripts/pairedIntSites.R')

sampleID  = as.integer(system("echo $LSB_JOBINDEX", intern=T))

parameters = get(load("parameters.RData"))[[sampleID]]

status = eval(as.call(append(list(processAlignments), unname(parameters[c("alias", "minPctIdent",
                                                                          "maxAlignStart", "maxFragLength")]))))

save(status, file="callStatus.RData") #working directory is changed while executing getTrimmedSeqs
