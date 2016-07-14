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
  standardized$Score <- NULL
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
  sites.reduced <- unlist(reduce(sites.reduced, min.gapwidth = 0L, with.revmap=TRUE))
  sites.reduced$counts <- sapply(sites.reduced$revmap, length)
  
  #Order original sites by revmap  
  dereplicatedSites <- sites[unlist(sites.reduced$revmap)]
  
  #Skip this step and provide similar output if length(sites) = 0
  if(length(sites) > 0){
    dereplicatedSites <- split(dereplicatedSites, Rle(values = seq(length(sites.reduced)), lengths = sites.reduced$counts))
  }  
  
  #Dereplicate reads with same standardized starts and provide the longeset width
  dereplicatedSites <- unlist(reduce(dereplicatedSites, min.gapwidth = 0L))
  mcols(dereplicatedSites) <- mcols(sites.reduced)

  dereplicatedSites
}
