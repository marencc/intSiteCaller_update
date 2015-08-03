#### libraries ####
libs <- (c("plyr", "ggplot2", "scales", "reshape2", "RMySQL", "knitr", "markdown"))
sapply(libs, require, character.only=TRUE)

options(stringsAsFactors = FALSE)

#### directory to process, current or given ####
cur_dir <- getwd()
args <- commandArgs(trailingOnly=TRUE)
if( length(args)>0 ) cur_dir <- args[1]

#### code Dir, from Rscript or search ####
codeDir <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly=FALSE), value=T)))
if( length(codeDir)!=1 ) codeDir <- list.files(path="~", pattern="intSiteCaller$", recursive=TRUE, include.dirs=TRUE, full.names=TRUE)
stopifnot(file.exists(file.path(codeDir, "intSiteCaller.R")))
stopifnot(file.exists(file.path(codeDir, "stats.Rmd")))

message("Processing ", cur_dir)

stats.file <- list.files(cur_dir, pattern="^stats.RData$", recursive=TRUE, full.names=TRUE)

tmp.statlist <- lapply(setNames(stats.file, stats.file), function(x) {
    a <- load(x)
    get(a)
})
tmp.statlist <- lapply(setNames(stats.file, stats.file), function(x) get(load(x)))
    
stats <- plyr:::rbind.fill(tmp.statlist)
stats$sample <- as.character(stats$sample)
rownames(stats) <- NULL

stats$gtsp <- NULL

#### get sampleInfo if available ####
gtsps <- unique(sub("-\\d+$", "", stats$sample))
getPatientInfo <- function(gtsps=gtsps) {
    junk <- sapply(dbListConnections(MySQL()), dbDisconnect)
    dbConn <- dbConnect(MySQL(), group="intSitesDev237") 
    patientInfo <- data.frame(gtsp=gtsps, info="")
    if( dbGetQuery(dbConn, "SELECT 1")==1 ) {
        gtspsin <- paste(sprintf("'%s'", gtsps), collapse = ",")
        sql <- sprintf("SELECT * FROM specimen_management.gtsp WHERE specimenaccnum IN (%s)", gtspsin)
        sampleInfo <- suppressWarnings( dbGetQuery(dbConn, sql) )
        colnames(sampleInfo) <- tolower(colnames(sampleInfo))
        
        patientInfo <- dplyr::select(sampleInfo,
                                     gtsp=specimenaccnum,
                                     trial,
                                     patient,
                                     timepoint,
                                     celltype,
                                     vcn,
                                     sampleprepmethod)
                                     ##seqmethod)
        patientInfo$info <- with(patientInfo, paste(trial, patient, timepoint, celltype, vcn, sampleprepmethod))
        patientInfo <- merge(data.frame(gtsp=gtsps), 
                             data.frame(gtsp=patientInfo$gtsp, info=patientInfo$info), 
                             all.x=TRUE)
    }
    return(patientInfo)    
}
gtspInfo <- getPatientInfo(gtsps = unique(sub("-\\d+$", "", stats$sample)))

stats.mdf <-  melt(stats, id.vars="sample")
stats.mdf$gtsp <- sub("-\\d+$", "", stats.mdf$sample)
stats.mdf$Replicate <- sub("GTSP\\d+-", "", stats.mdf$sample)

stats.mdf <- merge(stats.mdf, gtspInfo, by="gtsp", all.x=TRUE)
stats.mdf$gtspinfo <- with(stats.mdf, paste(gtsp, info))

stats.mdf.listBygtsp <- split(stats.mdf, stats.mdf$gtsp)

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
          legend.box = "horizontal")
          

plotsPerPage <- 2
plotList <- split(1:length(stats.mdf.listBygtsp), 
                  (1:length(stats.mdf.listBygtsp)-1)%/%plotsPerPage)

#### begin generating markdown ####
makeReport <- function() {
    RmdFile <- file.path(codeDir, "stats.Rmd")
    mdFile <- paste0(basename(cur_dir), ".stat.md")
    htmlFile <- paste0(basename(cur_dir), ".stat.html")
    pdfFile <- paste0(basename(cur_dir), ".stat.pdf")
    
    fig.path <- sprintf("%s.%s.%s.knitr.fig", basename(cur_dir),
                        format(Sys.Date(), format="%Y%m%d"),
                        Sys.getpid())
    
    unlink(fig.path, force=TRUE, recursive=TRUE)
    
    options(knitr.table.format='html')
    knit(RmdFile, output=mdFile)
    markdownToHTML(mdFile, htmlFile)
    message("output file: ", htmlFile, "\n")
    if( system("which wkhtmltopdf", ignore.stdout=TRUE, ignore.stderr=TRUE)==0 ) {
        cmd <- sprintf("wkhtmltopdf -s Letter %s %s", htmlFile, pdfFile)
        system(cmd)
        message("output file: ", pdfFile, "\n")
    }

    unlink(fig.path, force=TRUE, recursive=TRUE)
    unlink(mdFile, force=TRUE)
}
makeReport()
q()

#### saved test code ####
i=1
p <- ggplot(plyr::rbind.fill(stats.mdf.listBygtsp[ plotList[[i]] ]), 
            aes(variable, value, fill=Replicate)) +
    geom_bar(position=position_dodge(width = 0.8), stat="identity") + 
    scale_y_log10() +
    geom_vline(xintercept = 1:(ncol(stats)-2)+0.5, linetype=4) +
    theme_default 
print(p)
