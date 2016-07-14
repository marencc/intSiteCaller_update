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

options(stringsAsFactors = FALSE)

cur_dir <- getwd()
args <- commandArgs(trailingOnly=TRUE)
if( length(args)>0 ) cur_dir <- args[1]

stats.file <- list.files(cur_dir, pattern="^stats.RData$", recursive=TRUE, full.names=TRUE)

tmp.statlist <- lapply(setNames(stats.file, stats.file), function(x) {
    a <- load(x)
    get(a)
})
stats <- plyr:::rbind.fill(tmp.statlist)
stats$sample <- as.character(stats$sample)
rownames(stats) <- NULL

sampleinfo <- read.table("sampleInfo.tsv", header=TRUE)

samplestats <- merge(stats, sampleinfo, by.x="sample", by.y="alias", all.y=TRUE)

samplestats$workdir <- cur_dir
    
write.table(samplestats, "", sep = "\t", row.names=FALSE, quote=FALSE)

