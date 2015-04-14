library("sonicLength")
library("GenomicRanges")
library("ShortRead")

sampleID  = as.integer(system("echo $LSB_JOBINDEX", intern=T))

aliases = get(load("condenseMetadata.RData"))[sampleID]

condensedAlias = names(aliases)

aliases = aliases[[1]]

as = GRanges()
sites.final = GRanges()
stats = NULL
pidData = list(DNAStringSet(), BStringSet())
multihits = GRangesList()
chimeras = GRangesList()
for(alias in aliases){
  rep = which(aliases == alias)
  if(file.exists(paste0(alias, "/allSites.RData")) & length(get(load(paste0(alias, "/allSites.RData")))) > 1){
    allSites = get(load(paste0(alias, "/allSites.RData")))
    multihitData = get(load(paste0(alias, "/multihitData.RData")))
    chimeraData = get(load(paste0(alias, "/chimeraData.RData")))
    allSites$rep = rep
    as = append(as, allSites)
    primerIDData = get(load(paste0(alias, "/primerIDData.RData")))
    pidData[[1]] = append(pidData[[1]], primerIDData[[1]])
    pidData[[2]] = append(pidData[[2]], primerIDData[[2]])
    #if(lengt)
    multihits = append(multihits, multihitData[[2]])
    chimeras = append(chimeras, chimeraData[[4]])
    stats = rbind(stats, get(load(paste0(alias, "/stats.RData"))))
  }
}
suppressWarnings(dir.create(paste0("condensed/", condensedAlias), recursive=FALSE))
allSites = as
save(allSites, file=paste0("condensed/", condensedAlias, "/allSites.RData"))

primerIDData = pidData
save(primerIDData, file=paste0("condensed/", condensedAlias, "/primerIDData.RData"))

multihitData = multihits
save(multihitData, file=paste0("condensed/", condensedAlias, "/multhihitData.RData"))

chimeraData = chimeras
save(chimeraData, file=paste0("condensed/", condensedAlias, "/chimeraData.RData"))

stats = colSums(Filter(is.numeric, stats))
save(stats, file=paste0("condensed/", condensedAlias, "/stats.RData"))

uniqueSites = allSites
allSitesSolostart = flank(uniqueSites, -1, start=T, both=F)
allSitesSolostart = flank(allSitesSolostart, 5, both=TRUE)      
sites.reduced = reduce(allSitesSolostart, min.gapwidth=5, with.revmap=T)
sites.reduced$counts = sapply(sites.reduced$revmap, length)  

allSites = uniqueSites[unlist(sites.reduced$revmap)]
allSites = split(allSites, Rle(values=c(1:length(sites.reduced)), lengths=sites.reduced$counts))
allSites = unlist(reduce(allSites, min.gapwidth=5))
mcols(allSites) = mcols(sites.reduced)
sites.final = allSites
allSites = uniqueSites

sites.final$intLoc = start(flank(sites.final, width=-1, start=TRUE, both=FALSE))
#sites.final$Aliasposid = paste0(sites.final$clone, "_", as.character(seqnames(sites.final)), as.character(strand(sites.final)), sites.final$intLoc)
sites.final$posid = paste0(as.character(seqnames(sites.final)), as.character(strand(sites.final)), sites.final$intLoc)
sites.final$estAbund = NA #tmp for now, will be overwritten if sonicAbund is calculated

save(sites.final, file=paste0("condensed/", condensedAlias, "/sites.final.RData"))
