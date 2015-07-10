keptpattern <- c("*.txt$",
                 "*.csv$",
                 "*.tsv$",
                 "logs/",
                 "html$",
                 "*.RData")

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

