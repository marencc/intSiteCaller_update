options(stringsAsFactors = FALSE)

stats.file <- list.files(".", pattern="stats.RData", recursive=TRUE, full.names=TRUE)

junk <- lapply(setNames(stats.file, stats.file), function(x) {
    a <- load(x)
    get(a)
})
stats <- do.call(rbind, junk)
stats$sample <- as.character(stats$sample)
rownames(stats) <- NULL

sampleinfo <- read.table("sampleInfo.tsv", header=TRUE)

samplestats <- merge(stats, sampleinfo, by.x="sample", by.y="alias", all.y=TRUE)

write.table(samplestats, "", sep = "\t", row.names=FALSE, quote=FALSE)
