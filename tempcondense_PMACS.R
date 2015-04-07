library("sonicLength")
library("GenomicRanges")
library("ShortRead")
library("parallel")

sampleID  = as.integer(system("echo $LSB_JOBINDEX", intern=T))

aliases = get(load("condenseMetadata.RData"))[sampleID]

condensedAlias = names(aliases)

aliases = aliases[[1]]

as = GRanges()
sites.final = GRanges()
stats = NULL
pidData = list(DNAStringSet(), BStringSet())
for(alias in aliases){
  rep = strsplit(alias, "-")[[1]][[length(strsplit(alias, "-")[[1]])]]
  if(file.exists(paste0(alias, "/allSites.RData")) & length(get(load(paste0(alias, "/allSites.RData")))) > 1){
    allSites = get(load(paste0(alias, "/allSites.RData")))
    allSites$rep = rep
    as = append(as, allSites)
    primerIDData = get(load(paste0(alias, "/primerIDData.RData")))
    pidData[[1]] = append(pidData[[1]], primerIDData[[1]])
    pidData[[2]] = append(pidData[[2]], primerIDData[[2]])
    stats = rbind(stats, get(load(paste0(alias, "/stats.RData"))))
  }
}
suppressWarnings(dir.create(paste0("condensed/", condensedAlias), recursive=FALSE))
allSites = as
save(allSites, file=paste0("condensed/", condensedAlias, "/allSites.RData"))
primerIDData = pidData
save(primerIDData, file=paste0("condensed/", condensedAlias, "/primerIDData.RData"))
stats = colSums(Filter(is.numeric, stats))
save(stats, file=paste0("condensed/", condensedAlias, "/stats.RData"))

load(paste0("condensed/", condensedAlias, "/allSites.RData"))


uniqueSites = allSites
allSitesSolostart = flank(uniqueSites, -1, start=T, both=F)
allSitesSolostart = flank(allSitesSolostart, 5, both=TRUE)      
sites.reduced = reduce(allSitesSolostart, min.gapwidth=5, with.mapping=T)
sites.reduced$counts = sapply(sites.reduced$mapping, length)  

allSites = uniqueSites[unlist(sites.reduced$mapping)]
allSites = split(allSites, Rle(values=c(1:length(sites.reduced)), lengths=sites.reduced$counts))
allSites = unlist(reduce(allSites, min.gapwidth=5))
mcols(allSites) = mcols(sites.reduced)
sites.final = allSites
allSites = uniqueSites

sites.final$intLoc = start(flank(sites.final, width=-1, start=TRUE, both=FALSE))
#sites.final$Aliasposid = paste0(sites.final$clone, "_", as.character(seqnames(sites.final)), as.character(strand(sites.final)), sites.final$intLoc)
sites.final$posid = paste0(as.character(seqnames(sites.final)), as.character(strand(sites.final)), sites.final$intLoc)
sites.final$estAbund = NA #tmp for now, will be overwritten if sonicAbund is calculated

# 
# 
# 
# load(paste0("condensed/", condensedAlias, "/sites.final.RData"))
# load(paste0("condensed/", condensedAlias, "/allSites.RData"))
# 
allSites$realStart[unlist(sites.final$mapping)] = rep(sites.final$intLoc, sites.final$counts)
# 
save(allSites, file=paste0("condensed/", condensedAlias, "/allSites.RData"))
# 
# print("done reading files")
# 
# sites.final$posid = paste0(as.character(seqnames(sites.final)), as.character(strand(sites.final)), start(flank(sites.final, -1, start=T, both=F)))
# 
# allSites = split(allSites[unlist(sites.final$mapping)], rep(c(1:length(sites.final)), sapply(sites.final$mapping, length)))
# 
# allSites.chunks = split(allSites, ceiling(seq_along(allSites)/2000))
# 
# print("starting the parallel job")
# 
# allSites = unlist(mclapply(allSites.chunks, function(y){lapply(y, function(x){
#   newreps = rep(0, max(x$rep))
#   newreps[unique(x$rep)] = c(1:length(unique(x$rep)))
#   x$rep = newreps[x$rep]
#   x
# })}, mc.cores=15), use.names=F)
# 
# chunks = split(allSites, ceiling(seq_along(allSites)/3000))
# 
# print("starting sonic abund")
# 
# allAbunds = mclapply(chunks, function(x){
#   
#   sonicSites = unlist(GRangesList(x))
#   
#   pos = paste0(as.character(seqnames(sonicSites)), as.character(strand(sonicSites)), sonicSites$realStart)
#   
#   if(any(x$rep>1)){ #has reps
#     dfr = data.frame("posid"=pos, "qEnd"=width(sonicSites), "rep"=sonicSites$rep)
#   }  else{
#     dfr = data.frame("posid"=pos, "qEnd"=width(sonicSites))
#   }
#   
#   #dfr$Aliasposid = paste0(dfr$Alias, "_", as.character(seqnames(sonicSites)), as.character(strand(sonicSites)), pos)
#   dfr$posid = paste0(as.character(seqnames(sonicSites)), as.character(strand(sonicSites)), sonicSites$realStart)
#   dfr = droplevels(unique(dfr))
#   
#   if(any(x$rep>1)){ #has reps
#     res = estAbund(dfr$posid, dfr$qEnd, dfr$rep)
#   }  else{
#     res = estAbund(dfr$posid, dfr$qEnd)
#   }
#   
#   allAbund = round(res$theta)
#   #allAbundVar = res$var.theta
# }, mc.cores=14)
# 
# names(allAbunds) = NULL
# 
# allAbunds = unlist(allAbunds)
# 
# sites.final$estAbund = allAbunds[sites.final$posid]

# 
# #names(foo) = sites.final$posid
# 
# hits = findOverlaps(flank(sites.final, 5, both=TRUE), flank(allSites, 5, both=TRUE))
# 
# #crash
# #realStarts = data.frame(realStart=sites.final[queryHits(hits)]$intLoc)
# 
# #better, not sure why
# realStarts = data.frame(realStart=sites.final$intLoc[queryHits(hits)])
# 
# #rownames(realStarts) = names(allSites[subjectHits(hits)])
# realStarts$names = names(allSites[subjectHits(hits)])
# 
# #allSites$realStart = realStarts$realStart[match(names(allSites), rownames(realStarts))]
# allSites$realStart = realStarts$realStart[match(names(allSites), realStarts$names)]
# 
# 
if(TRUE){
  
  sonicSites = allSites
  
  pos = paste0(as.character(seqnames(sonicSites)), as.character(strand(sonicSites)), sonicSites$realStart)
  
  if(any(allSites$rep>1)){ #has reps
    dfr = data.frame("posid"=pos, "qEnd"=width(sonicSites), "rep"=sonicSites$rep)
  }  else{
    dfr = data.frame("posid"=pos, "qEnd"=width(sonicSites))
  }
  
  #dfr$Aliasposid = paste0(dfr$Alias, "_", as.character(seqnames(sonicSites)), as.character(strand(sonicSites)), pos)
  dfr$posid = paste0(as.character(seqnames(sonicSites)), as.character(strand(sonicSites)), sonicSites$realStart)
  dfr = droplevels(unique(dfr))
  
  if(any(allSites$rep>1)){ #has reps
    res = estAbund(dfr$posid, dfr$qEnd, dfr$rep)
  }  else{
    res = estAbund(dfr$posid, dfr$qEnd)
  }
  
  allAbund = round(res$theta)
  #allAbundVar = res$var.theta

  sites.final$estAbund = allAbund[sites.final$posid]
  #sites.final$varEstAbund = allAbundVar[sites.final$Aliasposid]
}

sites.final$pctAbund = (100*sites.final$estAbund)/sum(sites.final$estAbund)

save(sites.final, file=paste0("condensed/", condensedAlias, "/sites.final.RData"))
