#    This source code file is a component of the larger INSPIIRED genomic analysis software package.
#    Copyright (C) 2016 Frederic Bushman
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

library(dplyr)

cluster_multihits <- function(unclusteredMultihits, method, unique = TRUE){
  if(method == "gt.reduce.cluster"){
    flank.sites <- flank(unclusteredMultihits, -1, start = TRUE)
    red.sites <- reduce(flank.sites, min.gapwidth = 5L, with.revmap = TRUE)
    revmap <- red.sites$revmap
    edgelist <- unique(matrix(
      c(Rle(unclusteredMultihits$readPairKey[sapply(revmap, "[", 1)], 
            sapply(revmap, length)), 
        unclusteredMultihits$readPairKey[unlist(revmap)]), 
      ncol = 2))
    graph <- graph.edgelist(edgelist, directed = FALSE)
    clusters <- clusters(graph)
    unclusteredMultihits$cid <- clusters$membership[unclusteredMultihits$readPairKey]
    split(unclusteredMultihits, unclusteredMultihits$cid)
  }else if(method == "iter.reduce.cluster"){
    flank.sites <- flank(unclusteredMultihits, -1, start = TRUE)
    red.sites <- reduce(flank.sites, min.gapwidth = 5L, with.revmap=TRUE)
    pair.revmap <- lapply(red.sites$revmap, function(idx) unclusteredMultihits$readPairKey[idx])
    cid <- rep(NA, max(as.numeric(factor(unclusteredMultihits$readPairKey))))
    for(i in 1:length(pair.revmap)) {
      posidx <- pair.revmap[[i]]
      if( all(is.na(cid[posidx])) ) {
        cid[posidx] <- i
      } else {
        precid <- unique(cid[posidx])
        precid <- precid[!is.na(precid)]
        mincid <- min(precid)
        cid[cid %in% precid] <- mincid
        cid[posidx] <- mincid
      }
    }
    sites$cid <- as.numeric(factor(cid[unclusteredMultihits$readPairKey]))
    split(unclusteredMultihits, unclusteredMultihits$cid)
  }else if(method == "gt.findOverlaps.cluster"){
    if(unique){
      multihits.split <- unique(split(unclusteredMultihits, unclusteredMultihits$ID))
    }else{
      multihits.split <- split(unclusteredMultihits, unclusteredMultihits$ID)
    }
    multihits.split <- flank(multihits.split, -1, start=T)
    overlaps <- findOverlaps(multihits.split, multihits.split, maxgap=5)
    edgelist <- matrix(c(queryHits(overlaps), subjectHits(overlaps)), ncol=2)
    clusteredMultihitData <- clusters(graph.edgelist(edgelist, directed=F))
    clusteredMultihitNames <- split(names(multihits.split), clusteredMultihitData$membership)
    if(unique){
      clusteredMultihitPositions <- GRangesList(lapply(clusteredMultihitNames, function(x){
        unname(granges(unique(unlist(multihits.split[x]))))
      }))
    }else{
      clusteredMultihitPositions <- GRangesList(lapply(clusteredMultihitNames, function(x){
        unname(granges(unlist(multihits.split[x])))
      }))
    }
    clusteredMultihitPositions
  }
}
