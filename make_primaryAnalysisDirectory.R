#' 1. decide the runid from the current directory, e.g,
#'    within run20150609, the runid is 150609
#' 2. pull fastq.gz file from microb120:/media/THING1 for runid
#' 3. check if any .csv file is in the directory
#' 4. if found, output sampleInfo.tsv

options(stringsAsFactors = FALSE)


#### prepare Data directory ####
runid <- stringr::str_match(getwd(), "[1-2][4-9][0-1][0-9][0-3][0-9]$")
if( is.na(runid) | nchar(runid)!=6 ) stop(
                            "Can not determine runid, please run this script in a run directory such as run20150609" )
message("RunID: ", runid,"\n")

cat("enter microb120 user name, hit RETURN if same as on pmacs: ")
sshuser <- readLines(con="stdin", 1)
if( sshuser=="" ) sshuser <- Sys.getenv("USER")

cmd <- sprintf("ssh %s@microb120.med.upenn.edu ls /media/THING1/Illumina/%s_M*/Data/Intensities/BaseCalls/Undetermined_*.fastq.gz", sshuser, runid)

message("Cheching files")
message(cmd)
fastq <- system(cmd, intern=TRUE)
stopifnot(length(fastq)==3)
stopifnot(any(grepl("R1", fastq)))
stopifnot(any(grepl("R2", fastq)))
stopifnot(any(grepl("I1", fastq)))

dir.create("Data", showWarnings=FALSE)
cmd <- sprintf("scp %s@microb120.med.upenn.edu:%s Data", sshuser, fastq)
message("\n1. Downloading fastq.gz files")
gotfastq <- sapply(cmd, system)
if( !all(gotfastq==0) ) message("Connection to microb120 was not successful, get files manually")


#### prepare sampleInfo.tsv ####
needed <- c("alias", "linkerSequence", "bcSeq",	"gender", "primer", "ltrBit", "largeLTRFrag", "vectorSeq")

csv.file <- list.files(".", pattern="csv$")
if( length(csv.file)!=1 ) stop(
              "None or 1+ csv files found, quit, proceed manually")
if( !grepl(runid, csv.file) ) stop("csv filename must contain runid such as eneTherapy-20150505-sampleInfo.csv")

csv.tab <- read.csv(csv.file)
if(!all(needed %in% colnames(csv.tab) )) stop(paste(needed, collapse=" "), " colums are needed in the csv file")

tsv.tab <- subset(csv.tab, select=needed)
write.table(tsv.tab, file="sampleInfo.tsv", sep="\t", row.names=FALSE, quote=FALSE)
message("\n2. sampleInfo saved as tsv to sampleInfo.tsv\n")

#### get vector file ####
vectorfiles <- unique(tsv.tab$vectorSeq)

cat("enter svn user name, hit RETURN if same as on pmacs: ")
svnuser <- readLines(con="stdin", 1)
if( svnuser=="" ) svnuser <- Sys.getenv("USER")

cmd <- sprintf("wget --user=%s --ask-password https://microb98.med.upenn.edu/svn/yinghua/niravpipeline/vector/%s", svnuser, vectorfiles)
gotvector <- sapply(cmd, system)

if( !all(gotvector==0) ) stop("Connection to svn was not successful or vector files not found")

if( !all(file.exists(as.character(vectorfiles))) ) stop("fetching ", paste(vectorfiles, collapse=" "), "was not successful, proceed manually")

message("\n3. Got vector files, please procees and execute\n Rscript path/to/intSiteCaller.R")
  
