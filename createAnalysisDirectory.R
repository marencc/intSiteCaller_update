library('methods',  quietly=T)
library('RMySQL',   quietly=T)
library('optparse', quietly=T)
library('yaml')

# parse command line arguments
option_list = list( make_option( c("-s", "--sampleInfoFile"), action="store", default=NULL, type='character', help="path to the sample infor file"))
opt = parse_args(OptionParser(option_list=option_list))


# Check provided parameters
if (is.null(opt$s)) stop("A tab delimited sample information file needs to be provided via the -s or --sampleInfoFile flags")
if (!file.exists(opt$s)) stop("The sample information file does not exist.")


# Check for configuration file
if (!file.exists('./intSiteCallerConfig.yaml')) stop("The setup script could not find the expected intSiteCallerConfig.yaml configuration file.\n")


# load configuration file
config <- yaml.load_file("intSiteCallerConfig.yaml")


# determine the path to the directory that this script resides in
codeDir <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly=F), value=T)))


# Ensure that the user provided a sequencing run id in the config file
if (nchar(config$sampleRunId) < 1)  stop("A run id was not provided in the configuration file")


# Track runid via a data file
write(config$sampleRunId, file="miseqid.txt")


# Get all samples already in the database
null <- sapply(dbListConnections(MySQL()), dbDisconnect)
dbConn <- dbConnect(MySQL(), group=config$databaseConnectionGroup)
samples <- suppressWarnings( dbGetQuery(dbConn, "select * from samples") )


# Retrieve fastq files
dir.create("Data", showWarnings=F)

for (f in names(config[['SequencingFiles']]))
{
   cmdReturn <- vector()
  
   if (file.exists(config[['SequencingFiles']][[f]]))
   {
      cmdReturn <- sapply(paste0('cp ', config[['SequencingFiles']][[f]],  ' ./Data/Undetermined_S0_L001_', f, '_001.fastq.gz'), system)
   }
   else
   {
      cmdReturn <- sapply(paste0('scp ', config$remoteUser, ':', config[['SequencingFiles']][[f]],  ' ./Data/Undetermined_S0_L001_', f, '_001.fastq.gz'), system)
   }

   if (cmdReturn != 0) stop(paste0('The ', f, ' sequencing file could not be copied'))
   if (file.info(paste0('./Data/Undetermined_S0_L001_', f, '_001.fastq.gz'))$size == 0) stop('The ', f, ' sequencing file has zero size')
}


# Create processing paramaeter file
processingParamsHeader <- c('qualityThreshold', 'badQualityBases', 'qualitySlidingWindow', 'mingDNA', 'minPctIdent', 'maxAlignStart', 'maxFragLength', 'refGenome')

processingParamsValues <-  vector()
for (v in config$ProcessingParameters)
{
   processingParamsValues <- c(processingParamsValues, v)
}

write(processingParamsHeader, file='processingParams.tsv', sep ="\t", ncolumns=length(processingParamsHeader), append=F)
write(processingParamsValues, file='processingParams.tsv', sep ="\t", ncolumns=length(processingParamsHeader), append=T)
processingParams <- read.table("processingParams.tsv", header=TRUE)


# Required columns in the sample information file
needed <- c("alias", "linkerSequence", "bcSeq", "gender", "primer", "ltrBit", "largeLTRFrag", "vectorSeq")


# Fail safe, the sample run in needs to be in the csv file name
if( !grepl(config$sampleRunId, opt$s) ) stop("As a fail safe, the run id defined in the config file must be part of the sample information file")


# Read the sample information into a dataframe
csv.tab <- read.table(opt$s, sep="\t", header=TRUE)


# Make sure that the sample information file has tall the required columns 
if(!all(needed %in% colnames(csv.tab) )) stop(paste(needed, collapse=" "), " colums are needed in the csv file")


# Make sure that there are no extra columns outside of the required columns and columns found in the processing parameters
if( any(!colnames(csv.tab) %in% c(needed, colnames(processingParams)) ) ) {
    stop("Extra columns in csv file:\n", paste(setdiff(colnames(csv.tab), c(needed, colnames(processingParams))), collapse="\t") ) }


# Create a table from the provided sample information file using only the required columns defined in the 'needed'
tsv.tab <- subset(csv.tab, select=needed)


# Combine the processing param file information with the data from the sample information file.
# Add the run id to list. 
metadata <- merge(tsv.tab, processingParams)
metadata$miseqid <- config$sampleRunId


# Make sure that each sample alias is found only once in the sample information file
if( any(duplicated(tsv.tab$alias)) ) stop("The following alias are duplicated:\n",
    paste(tsv.tab$alias[duplicated(tsv.tab$alias)], collapse="\n"))


# Create a new table by joining together the metadata table and the samples table from the intSiteDB
# Join the tables on the sample name which is called 'alias' in the sample information table and called 'sampleName' in the intSiteDB.
conflict <- merge(metadata, samples, by.x="alias", by.y="sampleName")

# If a new table is created, then the sample names in the sample information file are already stored in the intSiteDB.

if( nrow(conflict) > 0 ) 
{
    # Create a sub-table with select columns
    conflict <-conflict[, c("alias", "refGenome.x", "refGenome.y", "miseqid.x", "miseqid.y")]

    # Present the table to the user
    message("Some alias names exist in the database, consider modify the ", opt$s, " file\n")
    write.table(conflict, sep="\t", quote=FALSE)
      
    # If any of the table rows have different run ids then the sample ids need to be incremented, ie. GTSP0567-1 -> GTSP0567-2
    if( any(conflict$miseqid.x != conflict$miseqid.y) ) message(
                "\nThe following replicates must be incremented, ",
                "GTSP numbers must be different for different runs.\n",
                paste(conflict$alias[conflict$miseqid.x!=conflict$miseqid.y], collapse="\n") )
    
    # If any of the table rows have the same run ids and the same  genome references then tell users that this appears to be a reprocessing request 
    if( any(conflict$refGenome.x==conflict$refGenome.y) ) message(
                "\nThe following replicates will be reprocessed\n",
                paste(conflict$alias[conflict$refGenome.x==conflict$refGenome.y], collapse="\n") )
    
    # If any of the rows have the same run ids but different genome references, tell the users that these may be more appropriate for a different DB.
    if( any(conflict$refGenome.x!=conflict$refGenome.y) ) message(
                "\nThe following replicates will be reprocessed on different freeze\n",
                paste(conflict$alias[conflict$refGenome.x!=conflict$refGenome.y], collapse="\n") )
}

write.table(tsv.tab, file="sampleInfo.tsv", sep="\t", row.names=FALSE, quote=FALSE)
write.table(processingParams, file="processingParams.tsv", sep="\t", row.names=FALSE, quote=FALSE)


# Retrieve vector files defined in sample information file
vectorfiles <- unique(tsv.tab$vectorSeq)

for (v in vectorfiles)
{
   cmdReturn <- vector()

   if (file.exists(paste0(config$vectorDataPath, '/', v)))
   {
      cmdReturn <- sapply(paste0('cp ', config$vectorDataPath, '/', v,  ' .'), system)
   }  
   else
   {
      cmdReturn <- sapply(paste0('scp ', config$remoteUser, ':', config$vectorDataPath, '/', v,  ' .'), system)
   }

   if (cmdReturn != 0) stop(paste0('The ', v, ' vector file could not be retrieved.'))
   if (file.info(v)$size == 0) stop('The ', v, ' vector file has zero size')
}

if( system("which tree > /dev/null 2>&1", ignore.stderr=TRUE)==0 ) system("tree")
message("\nNow ready to execute\n Rscript ", file.path(codeDir, "intSiteCaller.R"), " -j <distinctId>")
