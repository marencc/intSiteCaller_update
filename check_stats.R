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

