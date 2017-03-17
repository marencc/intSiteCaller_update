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

keptpattern <- c("*.txt$",
                 "*.csv$",
                 "*.tsv$",
                 "*.xls$",
                 "logs/",
                 "html$",
                 "pdf$",
                 "doc$",
                 "docx$",
                 "*.RData",
		 "*.yml")

allfiles <- list.files(recursive=TRUE)

filelist <- sapply(keptpattern, function(pattern) grep(pattern, allfiles, value=TRUE))

filekeep <- sort(unique(unlist(unname(filelist))))

filedel <- setdiff(allfiles, filekeep)

message("The following files will be deleted\n",
        paste(filedel, collapse="\n") )

cat("enter yes to continue: ")
yes <- readLines(con="stdin", 1)
if( yes!="yes" ) q()

file.remove(filedel)

