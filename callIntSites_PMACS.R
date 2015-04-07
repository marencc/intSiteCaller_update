source('~/EAS/PMACS_scripts/pairedIntSites.R')

sampleID  = as.integer(system("echo $LSB_JOBINDEX", intern=T))

parameters = get(load("parameters.RData"))[[sampleID]]

alias = parameters[[14]]

status = eval(as.call(append(list(processAlignments), parameters[c(12,14:16)])))
  
save(status, file=paste0(alias, "/callStatus.RData"))
