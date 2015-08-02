libs <- (c("plyr", "ggplot2", "reshape2", "knitr"))
sapply(libs, require, character.only=TRUE)

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

stats$gtsp <- NULL

stats.mdf <-  melt(stats, id.vars="sample")
stats.mdf$gtsp <- sub("-\\d+$", "", stats.mdf$sample)
stats.mdf$replicate <- sub("GTSP\\d+-", "", stats.mdf$sample)

stats.mdf.samplelist <- split(stats.mdf, stats.mdf$gtsp)

theme_default <- theme_bw() + 
    theme(text = element_text(size=14),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.x = element_text(size=14, face="bold", angle = 45, hjust = 1),
          axis.title.x = element_text(size=14, face="bold"),
          axis.text.y = element_text(size=14, face="bold"),
          axis.title.y = element_text(size=14, face="bold"),
          ##legend.position=c(0.9,0.85),
          ##legend.key.size=1,
          legend.text=element_text(size=8),
          legend.position="top",
          legend.box = "horizontal",
          legend.title = element_text("Replicates"))


ggplot(plyr::rbind.fill(stats.mdf.samplelist[1:4]), aes(variable, value, fill=replicate)) +
    geom_bar(position=position_dodge(width = 0.8), stat="identity") + 
            theme_default +
                facet_wrap( ~ gtsp, ncol=1, scales = "free_y")


ggplot(plyr::rbind.fill(stats.mdf.samplelist[ plotList[[i]] ]), 
       aes(variable, value, fill=replicate)) +
    geom_bar(position=position_dodge(width = 0.8), stat="identity") + 
    theme_default +
    facet_wrap( ~ gtsp, ncol=1, scales = "free_y")


plotsPerPage <- 4
plotList <- split(1:length(stats.mdf.samplelist), 
                  (1:length(stats.mdf.samplelist)-1)%/%plotsPerPage)

for(i in seq(plotList)) {
    ##cat(plotList[[i]],"\n")
    
    ggplot(plyr::rbind.fill(stats.mdf.samplelist[ plotList[[i]] ]), 
           aes(variable, value, fill=replicate)) +
        geom_bar(position=position_dodge(width = 0.8), stat="identity") + 
        theme_default +
        facet_wrap( ~ gtsp, ncol=1, scales = "free_y")
    
}
