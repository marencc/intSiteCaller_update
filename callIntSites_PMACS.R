source('~/EAS/PMACS_scripts/pairedIntSites.R')

sampleID  = as.integer(system("echo $LSB_JOBINDEX", intern=T))

parameters = get(load("parameters.RData"))[[sampleID]]

alias = parameters[[14]]

alignMethod = "BLAT"

status = eval(as.call(append(list(processAlignments, alignMethod, T), parameters[c(14,16:18,8,13,19,20)])))

save(status, file=paste0(alias, "/callStatus.RData"))
