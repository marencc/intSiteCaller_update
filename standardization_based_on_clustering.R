standardizeSites <- function(unstandardizedSites){
  if( ! length(unstandardizedSites) > 0){
  return(unstandardizedSites)
  }
  #Get called start values for clustering  
  unstandardizedSites$Position <- ifelse(strand(unstandardizedSites) == "+", start(unstandardizedSites), end(unstandardizedSites))
  unstandardizedSites$Break <- ifelse(strand(unstandardizedSites) == "+", end(unstandardizedSites), start(unstandardizedSites))
  unstandardizedSites$Score <- 95
  unstandardizedSites$qEnd <- width(unstandardizedSites)
  
  #Positions clustered by 5L window and best position is chosen for cluster
  standardized <- clusterSites(
    psl.rd = unstandardizedSites,
    weight = rep(1, length(unstandardizedSites)) 
    )

  start(standardized) <- ifelse(strand(standardized) == "+", 
                                standardized$clusteredPosition, standardized$Break)
  end(standardized) <- ifelse(strand(standardized) == "-", 
                              standardized$clusteredPosition, standardized$Break)
  
  standardized$Position <- NULL
  standardized$Break <- NULL
  standardized$score <- NULL
  standardized$qEnd <- NULL
  standardized$clusteredPosition <- NULL
  standardized$clonecount <- NULL
  standardized$clusterTopHit <- NULL
  
  sort(standardized)
}  

  dereplicateSites <- function(sites){
    #Reduce sites which have the same starts, but loose range info
    #(no need to add a gapwidth as sites are standardized)
    sites.reduced <- flank(sites, -1, start=TRUE)
    sites.reduced <- unlist(reduce(sites.reduced, with.revmap=TRUE))
    sites.reduced$counts <- sapply(sites.reduced$revmap, length)
    
    #Order original sites by revmap  
    dereplicatedSites <- sites[unlist(sites.reduced$revmap)]
    
    #Skip this step and provide similar output if length(sites) = 0
    if(length(sites) > 0){
      dereplicatedSites <- split(dereplicatedSites, Rle(values = seq(length(sites.reduced)), lengths = sites.reduced$counts))
    }  
    
    #Dereplicate reads with same standardized starts and provide the longeset width
    dereplicatedSites <- unlist(reduce(dereplicatedSites))
    mcols(dereplicatedSites) <- mcols(sites.reduced)

    dereplicatedSites
  }
