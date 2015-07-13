  dereplicateSites_overlap <- function(uniqueReads){
    #do the dereplication, but loose the coordinates
    sites.reduced <- reduce(flank(uniqueReads, -5, both=TRUE), with.revmap=T)
    sites.reduced$counts <- sapply(sites.reduced$revmap, length)
    
    #order the unique sites as described by revmap
    dereplicatedSites <- uniqueReads[unlist(sites.reduced$revmap)]
    
    #if no sites are present, skip this step - keep doing the rest to provide a
    #similar output to a successful dereplication
    if(length(uniqueReads)>0){
      #split the unique sites as described by revmap (sites.reduced$counts came from revmap above)
      dereplicatedSites <- split(dereplicatedSites, Rle(values=seq(length(sites.reduced)), lengths=sites.reduced$counts))
    }
    
    #do the standardization - this will pick a single starting position and
    #choose the longest fragment as ending position
    dereplicatedSites <- unlist(reduce(dereplicatedSites, min.gapwidth=5))
    mcols(dereplicatedSites) <- mcols(sites.reduced)
    
    dereplicatedSites
  }
  
  standardizeSites_overlap <- function(unstandardizedSites){
    if(length(unstandardizedSites)>0){
      dereplicated <- dereplicateSites(unstandardizedSites)
      dereplicated$dereplicatedSiteID <- seq(length(dereplicated))
      
      #order the original object to match
      unstandardizedSites <- unstandardizedSites[unlist(dereplicated$revmap)]
      
      #graft over the seqnames, starts, ends, and metadata
      trueBreakpoints <- start(flank(unstandardizedSites, -1, start=F))
      ##standardizedStarts <- rep(start(dereplicated), dereplicated$counts)
      standardizedStarts <- rep(start(flank(dereplicated, -1, start=T)), dereplicated$counts)
      standardized <- GRanges(seqnames=seqnames(unstandardizedSites),
                              ranges=IRanges(start=pmin(standardizedStarts, trueBreakpoints),
                                             end=pmax(standardizedStarts, trueBreakpoints)),
                              strand=strand(unstandardizedSites),
                              seqinfo=seqinfo(unstandardizedSites))
      mcols(standardized) <- mcols(unstandardizedSites)
      
      standardized
    }
    else{
      unstandardizedSites
    }
  }
