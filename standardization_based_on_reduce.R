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
