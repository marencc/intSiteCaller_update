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

library(digest)
options(stringsAsFactors = FALSE)

wideScreen <- function(howWide=as.numeric(strsplit(system('stty size', intern=T), ' ')[[1]])[2]) {
   options(width=as.integer(howWide))
}

rdate2md5 <- function(rdata_file) {
    stopifnot(file.exists(rdata_file))
    
    junk <- sapply(RData_objects, function(x){
        object <- load(x)
        md5 <- digest:::digest(get(object))
        rm(list=object)
        return(md5)
    } )
    df <- data.frame(RData=names(junk), digest=unname(junk))
    return(df)
}

old.digest <- "intSiteValidation.digest"

RData_objects <- sort(list.files(path=".", pattern="RData$", recursive=TRUE))

md5.new <- rdate2md5(RData_objects)

if( !file.exists(old.digest) ) {
    write.table(file=old.digest, md5.new, row.names=FALSE, quote=FALSE, sep="\t")
    message("No prvious digest found, digest were generated and saved in ", old.digest)
    q()
}

md5.old <- read.table(old.digest, header=TRUE)

md5.all <- merge(md5.old, md5.new, by="RData", all.x=TRUE, all.y=TRUE)
colnames(md5.all) <- c("RData", "digest.old", "digest.new")
md5.all <- subset(md5.all, !grepl("codeDir.RData|completeMetadata.RData|primerIDData.RData",
                                 RData,
                                 ignore.case=TRUE) )

md5.all$same <- with(md5.all, digest.old==digest.new)

if( any(!md5.all$same) ) {
    message("FAIL")
    message("The following RData's are different")
    wideScreen()
    print(subset(md5.all, !same))
} else {
    message("PASS")
}

