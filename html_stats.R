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

#### libraries ####
libs <- c("dplyr", "ggplot2", "reshape2", "scales", "RMySQL", "knitr", "markdown", "yaml", "DBI")
null <- suppressMessages(sapply(libs, library, character.only=TRUE))

options(stringsAsFactors = FALSE)
options(dplyr.width = Inf)

# load configuration file
config <- yaml.load_file("INSPIIRED.yml")

get_args <- function() {
    suppressMessages(library(argparse))
    
    codeDir <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly=FALSE), value=T)))
    if( length(codeDir)!=1 ) codeDir <- list.files(path="~", pattern="intSiteCaller$", recursive=TRUE, include.dirs=TRUE, full.names=TRUE)
    stopifnot(file.exists(file.path(codeDir, "intSiteCaller.R")))
    stopifnot(file.exists(file.path(codeDir, "stats.Rmd")))
    
    p <- ArgumentParser(formatter_class='argparse.RawTextHelpFormatter')
    p$add_argument("dataDir", nargs='?', default='.')
    p$add_argument("-c", "--codeDir", type="character", nargs=1,
                   default=codeDir,
                   help="Directory of code")
    args <- p$parse_args(commandArgs(trailingOnly=TRUE))
    
    args$dataDir <- normalizePath(args$dataDir, mustWork=TRUE)
    
    return(args)
}
args <- get_args()

#### load stats.RData, sampleInfo.tsv ####
message("\nProcessing ", args$dataDir)
sampleInfo <- read.table(file.path(args$dataDir, "sampleInfo.tsv"), header=TRUE)

stats.file <- list.files(args$dataDir, pattern="^stats.RData$", recursive=TRUE, full.names=TRUE)

stats <- plyr:::rbind.fill( lapply(stats.file, function(x) get(load(x))) )

# JKE
#message('JKE:')
#write.table(stats)


stats$sample <- as.character(stats$sample)
rownames(stats) <- NULL

stats <- merge(data.frame(sample=as.character(sampleInfo$alias)),
               plyr:::rbind.fill(lapply(stats.file, function(x) get(load(x)))),
               by="sample",
               all.x=TRUE)

#write.table(stats)


stats$sample <- as.character(stats$sample)
stats$gtsp <- sub("-\\d+$", "", stats$sample)
stats$Replicate <- sub("GTSP\\d+-", "", stats$sample)
rownames(stats) <- NULL
stats[is.na(stats)] <- 0 # otherwise NA + 1 = NA

#### get sampleInfo if available ####
##gtsps <- unique(sub("-\\d+$", "", stats$sample))
getPatientInfo <- function(gtsps=gtsps) {
    junk <- sapply(dbListConnections(MySQL()), dbDisconnect)

    ### dbConn <- dbConnect(MySQL(), group="intsites_miseq.read") 

    if (config$dataBase == "mysql"){
      dbConn <- dbConnect(MySQL(), group=config$mysqlSpecimenManagementGroup)
    }else{
      dbConn <- dbConnect(RSQLite::SQLite(), dbname=config$SQLiteDBpath) }



    patientInfo <- data.frame(gtsp=gtsps, info="")
    gtspsin <- paste(sprintf("'%s'", gtsps), collapse = ",")
    sql <- sprintf("SELECT * FROM specimen_management.gtsp WHERE specimenaccnum IN (%s)", gtspsin)
    sampleInfo <- suppressWarnings( dbGetQuery(dbConn, sql) )

    # JKE
    message('gtspsin: ', gtspsin)
    write.table(sampleInfo)

    colnames(sampleInfo) <- tolower(colnames(sampleInfo))
    
    patientInfo <- dplyr::select(sampleInfo,
                                 gtsp=specimenaccnum,
                                 trial,
                                 patient,
                                 timepoint,
                                 celltype,
                                 vcn,
                                 sampleprepmethod)
    
    patientInfo$info <- with(patientInfo, paste(trial, patient, timepoint, celltype, vcn, sampleprepmethod))
    patientInfo <- merge(data.frame(gtsp=gtsps), 
                         data.frame(gtsp=patientInfo$gtsp, info=patientInfo$info), 
                         all.x=TRUE)
    return(patientInfo)    
}
##gtspInfo <- getPatientInfo(gtsps = unique(sub("-\\d+$", "", stats$sample)))
gtspInfo <- try( getPatientInfo(gtsps=unique(stats$gtsp)) )
if( "try-error" %in% class(gtspInfo) ) {
    gtspInfo <- data.frame(gtsp=unique(stats$gtsp),
                           info="NA")
}
stats$gtsp <- sub("-\\d+$", "", stats$sample)
stats <- merge(stats, gtspInfo, by="gtsp", all.x=TRUE)

#### prepare melted data frame to be used by ggplot ####
plotCols <- ! colnames(stats) %in% c("gtsp", "Replicate", "info")
stats.mdf <- melt(subset(stats, select=plotCols), id.vars="sample")
stats.mdf$gtsp <- sub("-\\d+$", "", stats.mdf$sample)
stats.mdf$Replicate <- sub("GTSP\\d+-", "", stats.mdf$sample)
stats.mdf <- merge(stats.mdf, gtspInfo, by="gtsp", all.x=TRUE)
stats.mdf$gtspinfo <- with(stats.mdf, paste(gtsp, info))
stats.mdf$info <- NULL

## split by gtsp number
stats.mdf.listBygtsp <- split(stats.mdf, stats.mdf$gtsp)

#### set default print theme ####
theme_default <- function() {
    theme0 <- theme_bw() + 
        theme(text = element_text(size=14),
              axis.title.x=element_blank(),
              axis.title.y=element_blank(),
              axis.text.x = element_text(size=14, face="bold", angle = 45, hjust = 1),
              axis.title.x = element_text(size=14, face="bold"),
              axis.text.y = element_text(size=14, face="bold"),
              axis.title.y = element_text(size=14, face="bold"),
              legend.text=element_text(size=8),
              legend.position="top",
              legend.box = "horizontal")
    
    return(theme0)
}

#### begin generating markdown ####
makeReport <- function() {
    RmdFile <- file.path(args$codeDir, "stats.Rmd")
    mdFile <- paste0(basename(args$dataDir), ".stat.md")
    htmlFile <- paste0(basename(args$dataDir), ".stat.html")
    pdfFile <- paste0(basename(args$dataDir), ".stat.pdf")
    
    fig.path <- sprintf("%s.%s.%s.knitr.fig", basename(args$dataDir),
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
save.image(file="debug.RData")
q()

#### saved test code ####
i=1
plotList[[i]]

pl <- plotList.gtsp[[1]]

p <- ggplot(plyr::rbind.fill(stats.mdf.listBygtsp[ pl ]), 
            aes(variable, value, fill=Replicate)) +
    geom_bar(position=position_dodge(width = 0.8), stat="identity") + 
    scale_y_log10() +
    geom_vline(xintercept = 1:(nlevels(stats.mdf$variable)-1) +0.5, linetype=4) +
    theme_default() 
print(p)
