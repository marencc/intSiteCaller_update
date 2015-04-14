cleanup = get(load("cleanup.RData"))

if(cleanup == "full" && file.exists("condensed")){ #do this first so that parameters.RData still exists
  parameters = get(load("parameters.RData"))
  aliases = sapply(parameters, "[[", 12)
  system(paste0("rm -r ", paste0(aliases, collapse=" ")), ignore.stderr=T)
}

if(cleanup != "none"){
  system("rm *.RData", ignore.stderr=T)
  system("rm Data/*.fasta", ignore.stderr=T)
  system("rm */hits.R*.RData", ignore.stderr=T)
  system("rm */p*.fa*", ignore.stderr=T)
  system("rm */keys.RData", ignore.stderr=T)
  system("rm */*Status.RData", ignore.stderr=T)
  system("rm -r logs", ignore.stderr=T)
  system("rm -r Data/demultiplexedReps", ignore.stderr=T)
}
