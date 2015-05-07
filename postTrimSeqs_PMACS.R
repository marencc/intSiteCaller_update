#sample index in callIntSites_PMACS.R is based on parameters.RData, so we
#should pass appropriate index values by using parameters.RData here
metadata = get(load("parameters.RData"))

bushmanJobID = get(load("bushmanJobID.RData"))
blatStartPort = get(load("bushmanBlatStartPort.RData"))

numAliases  = nrow(metadata)

numFastaFiles = length(system("ls */*.fa", intern=T))

#align seqs
system(paste0('bsub -n1 -q normal -w "done(BushmanPostTrimProcessing_', bushmanJobID, ')" -J "BushmanAlignSeqs_', bushmanJobID, '[', blatStartPort, '-', blatStartPort+numFastaFiles-1, ']" -o logs/alignOutput%I.txt Rscript ~/EAS/PMACS_scripts/alignSeqs_PMACS.R'))
system("sleep 10")
#call int sites (have to find out which ones worked)
successfulTrims = unlist(sapply(metadata, function(x){
  get(load(paste0(x[["alias"]], "/trimStatus.RData"))) == paste0(getwd(), "/", x[["alias"]])
}))

system(paste0('bsub -n1 -q max_mem30 -w "done(BushmanAlignSeqs_', bushmanJobID, ')" -J "BushmanCallIntSites_', bushmanJobID, '[', paste(which(successfulTrims), collapse=","), ']" -o logs/callSitesOutput%I.txt Rscript ~/EAS/PMACS_scripts/callIntSites_PMACS.R'))

system(paste0('bsub -n1 -q normal -w "done(BushmanCallIntSites_', bushmanJobID, ')" -J "BushmanCleanup_', bushmanJobID, '" -o logs/cleanupOutput.txt Rscript ~/EAS/PMACS_scripts/cleanup_PMACS.R'))
