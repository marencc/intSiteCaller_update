#' 1. for a given  runid build directory for caller
#' 2. pull fastq.gz file from microb120:/media/THING1 for runid
#' 3. check if any .csv file is in the directory
#' 4. if found, output sampleInfo.tsv

#' at present script only works for microb120 and user yinghua
#' location of fastq files  and database is hard-coded
runid <- commandArgs(trailingOnly=TRUE)[1]

library("methods", quietly=TRUE)
library("RMySQL", quietly = TRUE) #also loads DBI
options(stringsAsFactors = FALSE)

group <- "intsites_miseq.read"
VECTOR_ROOT_PATH <- "/media/THING1/dryga/intSitePipe/vector"

codeDir <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly=FALSE), value=T)))
if( length(codeDir)!=1 ) codeDir <- list.files(path="~", pattern="intSiteCaller$", recursive=TRUE, include.dirs=TRUE, full.names=TRUE)
if( length(codeDir)!=1 ) stop("Cannot determine codeDir")
stopifnot(file.exists(file.path(codeDir, "default_processingParams.tsv")))
stopifnot(file.exists(file.path(codeDir, "intSiteCaller.R")))

#### get all samples already in the database ####
null <- sapply(dbListConnections(MySQL()), dbDisconnect)
dbConn <- dbConnect(MySQL(), group=group)
sql <- "select * from samples"
samples <- suppressWarnings( dbGetQuery(dbConn,sql) )

message("Run Directory: ", runid,"\n")

sshuser <- Sys.getenv("USER")

cmd <- sprintf("ssh %s@microb120.med.upenn.edu ls /media/THING1/Illumina/%s/Data/Intensities/BaseCalls/Undetermined_*.fastq.gz", sshuser, runid)

message("Checking files")
message(cmd)
fastq <- system(cmd, intern=TRUE)
if( length(fastq)<3 ) {
    cmd <- sprintf("ssh %s@microb120.med.upenn.edu ls /media/THING1/Illumina/%s/Undetermined_*.fastq.gz", sshuser, runid)
    message(cmd)
    fastq <- system(cmd, intern=TRUE)
}
fastq <- grep("M03249|M00142|M03543", fastq, value=TRUE)


if( ! length(fastq)==3 ) stop( paste(fastq, collapse="\n") )
stopifnot(length(fastq)==3)
stopifnot(any(grepl("R1", fastq)))
stopifnot(any(grepl("R2", fastq)))
stopifnot(any(grepl("I1", fastq)))
message("\nFastq files:\n", paste(fastq, collapse="\n"))

stopifnot(length(runid)==1)
message("RunID:\t", runid)
write(runid, file="miseqid.txt")
message("\n1. miseqid saved as miseqid.txt\n")

#### prepare processingParams.tsv ####
if( !file.exists("processingParams.tsv") ) {
    stopifnot( file.copy(file.path(codeDir, "default_processingParams.tsv"),
                         file.path(".", "processingParams.tsv") ) )
}
processingParams <- read.table("processingParams.tsv", header=TRUE)
message("\n2. processing parameters saved as processingParams.tsv\n")

#### prepare sampleInfo.tsv, check if alias names already in database ####
needed <- c("alias", "linkerSequence", "bcSeq",	"gender", "primer", "ltrBit", "largeLTRFrag", "vectorSeq")

csv.file <- list.files(".", pattern="csv$")
if( length(csv.file)!=1 ) stop("\n",
              paste(csv.file, collapes="\n"),
              "\nNone or 1+ csv files found, quit, proceed manually")
if( !grepl(runid, csv.file) ) stop("csv filename must contain runid")

csv.tab <- read.csv(csv.file)
if(!all(needed %in% colnames(csv.tab) )) stop(paste(needed, collapse=" "), " colums are needed in the csv file")

if( any(!colnames(csv.tab) %in% c(needed, colnames(processingParams)) ) ) {
    message("Extra colums in csv file:\n", paste(
        setdiff(colnames(csv.tab), c(needed, colnames(processingParams))), collapse="\t") )
}


## set csv priority if common fields exist
if( any(colnames(csv.tab) %in% colnames(processingParams)) ) {
    comField <- intersect(colnames(csv.tab), colnames(processingParams))
    needed <- unique(c(needed, comField))
    processingParams <- processingParams[, -which(colnames(processingParams) %in% comField)]
}

tsv.tab <- subset(csv.tab, select=needed)

metadata <- merge(tsv.tab, processingParams)
metadata$miseqid <- runid

if( any(duplicated(tsv.tab$alias)) ) stop(
    "\nThe following alias are duplicated:\n",
    paste(tsv.tab$alias[duplicated(tsv.tab$alias)], collapse="\n"),
    "\n")

conflict <- merge(metadata, samples, by.x="alias", by.y="sampleName")
if( nrow(conflict) > 0 ) {
    conflict <-conflict[, c("alias", "refGenome.x", "refGenome.y",
                            "miseqid.x", "miseqid.y")]
    message("Some alias names exist in the database, consider modify the ", csv.file, " file\n")
    write.table(conflict, sep="\t", quote=FALSE)
    message()
    if( any(conflict$miseqid.x!=conflict$miseqid.y) ) message(
                "\nThe following replicates must be incremented, ",
                "GTSP numbers must be different for different runs.\n",
                paste(conflict$alias[conflict$miseqid.x!=conflict$miseqid.y], collapse="\n") )
    
    if( any(conflict$refGenome.x==conflict$refGenome.y) ) message(
                "\nThe following replicates will be reprocessed\n",
                paste(conflict$alias[conflict$refGenome.x==conflict$refGenome.y], collapse="\n") )
    
    if( any(conflict$refGenome.x!=conflict$refGenome.y) ) message(
                "\nThe following replicates will be reprocessed on different freeze\n",
                paste(conflict$alias[conflict$refGenome.x!=conflict$refGenome.y], collapse="\n") )
    
}

write.table(tsv.tab, file="sampleInfo.tsv", sep="\t", row.names=FALSE, quote=FALSE)
write.table(processingParams, file="processingParams.tsv", sep="\t", row.names=FALSE, quote=FALSE)
message("\n3. sampleInfo saved as tsv to sampleInfo.tsv\n")

### get Data/fastq files ####
dir.create("Data", showWarnings=FALSE)
cmd <- sprintf("scp %s@microb120.med.upenn.edu:%s Data", sshuser, fastq)
message("\nDownloading fastq.gz files")
gotfastq <- sapply(cmd, system)
if( !all(gotfastq==0) ) message("Connection to microb120 was not successful, get files manually")
message("\n4. fastq saved as in Data\n")

#### get vector file ####
vectorfiles <- unique(tsv.tab$vectorSeq)
cmd <- sprintf(paste0("scp %s@microb120.med.upenn.edu:", VECTOR_ROOT_PATH, "/%s ."), sshuser, vectorfiles)
gotvector <- sapply(cmd, system)
if( !all(gotvector==0) ) stop("Connection to svn was not successful or vector files not found")

if( !all(file.exists(as.character(vectorfiles))) ) stop("fetching ", paste(vectorfiles, collapse=" "), "was not successful, proceed manually")

message("\n5. vector saved in", paste(vectorfiles, collapse="\n"))

#### done ####
if( system("which tree > /dev/null 2>&1", ignore.stderr=TRUE)==0 ) system("tree")
message("\nNow ready to execute\n Rscript ",
        file.path(codeDir, "intSiteCaller.R"),
        " -j run", runid)


