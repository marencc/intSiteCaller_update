metadata = read.csv("allRepsBCs.csv")

bushmanJobID = get(load("bushmanJobID.RData"))
blatStartPort = get(load("bushmanBlatStartPort.RData"))

numAliases  = nrow(metadata)

numFastaFiles = length(system("ls */*.fa", intern=T))

#align seqs
system(paste0('bsub -n1 -q normal -w "done(BushmanPostTrimProcessing_', bushmanJobID, ')" -J "BushmanAlignSeqs_', bushmanJobID, '[', blatStartPort, '-', blatStartPort+numFastaFiles-1, ']" -o logs/alignOutput%I.txt Rscript ~/EAS/PMACS_scripts/alignSeqs_PMACS.R'))
system("sleep 10")
#call int sites (have to find out which ones worked)
successfulTrims = unlist(sapply(metadata$Alias, function(x){
  get(load(paste0(x, "/trimStatus.RData"))) == paste0(getwd(), "/", x)
}))

system(paste0('bsub -n1 -q max_mem30 -w "done(BushmanAlignSeqs_', bushmanJobID, ')" -J "BushmanCallIntSites_', bushmanJobID, '[', paste(which(successfulTrims), collapse=","), ']" -o logs/callSitesOutput%I.txt Rscript ~/EAS/PMACS_scripts/callIntSites_PMACS.R'))

if("Condensed.Alias" %in% names(metadata)){
  suppressWarnings(dir.create("condensed/"))

  metadata = metadata[successfulTrims,]
  
  condenseMetadata = lapply(split(metadata$Sample.Name, metadata$Condensed.Alias), as.character)
  condenseMetadata = condenseMetadata[sapply(condenseMetadata, length)>1]
  save(condenseMetadata, file="condenseMetadata.RData")
  
  system(paste0('bsub -n1 -q max_mem30 -w "done(BushmanCallIntSites_', bushmanJobID, ')" -J "BushmanCondense_', bushmanJobID, '[1-', length(condenseMetadata), ']" -o logs/condenseOutput%I.txt Rscript ~/EAS/PMACS_scripts/condense_PMACS.R'))
  
  #system(paste0('bsub -n1 -q max_mem30 -J "BushmanCondense_', bushmanJobID, '[1-', length(condenseMetadata), ']" -o logs/condenseOutput%I.txt Rscript ~/EAS/PMACS_scripts/condense_PMACS.R'))
  
  
}