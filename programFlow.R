bsub = function(cpus = 1, queue = "normal", wait = NULL, jobName = NULL, logFile = NULL, command = NULL){ 
  cmd = paste0("bsub ",
               "-n", as.character(cpus),
               " -q ", queue)
  
  if(!is.null(wait)){
    cmd = paste0(cmd, " -w \"", wait, "\"")
  }
  
  if(!is.null(jobName)){
    cmd = paste0(cmd, " -J \"", jobName, "\"")
  }
  
  if(!is.null(logFile)){
    cmd = paste0(cmd, " -o ", logFile)
  }
  
  cmd = paste0(cmd, " ", command) #no default, should crash if no command provided
  system(cmd)
}

alignSeqs = function(){
  # do alignments
  library("rtracklayer")
  library("GenomicRanges")
  
  blatPort  = as.integer(system("echo $LSB_JOBINDEX", intern=T))
  
  toAlign = system("ls */*.fa", intern=T)
  
  alignFile = toAlign[blatPort-get(load("bushmanBlatStartPort.RData"))+1] #blatPort is actually file offset by first blat port number and R is offset by 1
  
  alias = strsplit(alignFile, "/")[[1]][1]
  
  metadata = read.csv("processingParams.csv")
  
  #we'll still hardcode the index directory path since it would be a pain to load
  #the whole genome from an R package into a .2bit file for BLAT
  indexPath = paste0("/home/aubreyba/genomeIndices/", metadata[metadata$alias==alias,"refGenome"], ".2bit")
  
  system(paste0("/home/aubreyba/EAS/PMACS_scripts/BLATsamples.sh ", alignFile, " ", blatPort, " ", indexPath))
  
}

callIntSites = function(){
  source('~/EAS/PMACS_scripts/pairedIntSites.R')
  
  sampleID  = as.integer(system("echo $LSB_JOBINDEX", intern=T))
  
  parameters = get(load("parameters.RData"))[[sampleID]]
  
  status = eval(as.call(append(list(processAlignments), unname(parameters[c("alias", "minPctIdent",
                                                                            "maxAlignStart", "maxFragLength",
                                                                            "refGenome")]))))
  
  save(status, file="callStatus.RData") #working directory is changed while executing getTrimmedSeqs
}

cleanup = function(){
  cleanup = get(load("cleanup.RData"))
  
  if(cleanup){
    system("rm *.RData", ignore.stderr=T)
    system("rm Data/*.fasta", ignore.stderr=T)
    system("rm */hits.R*.RData", ignore.stderr=T)
    system("rm */R*.fa*", ignore.stderr=T)
    system("rm */keys.RData", ignore.stderr=T)
    system("rm */*Status.RData", ignore.stderr=T)
    system("rm -r logs", ignore.stderr=T)
    system("rm -r Data/demultiplexedReps", ignore.stderr=T)
  }
}

demultiplex = function(){
  #Demultiplexing is currently a single-core process - perhaps it could be made more efficient by having each 
  #error-correct worker do its own mini-demultiplex with the barcodes that it error corrected, then write their
  #own shorter fastq files which can be cat'd together after everything is done
  
  library("ShortRead")
  
  I1 = readFasta(list.files("Data", pattern="correctedI1-.", full.names=T)) #readFasta("Data/correctedI1.fasta")
  
  metadata = read.csv("sampleInfo.csv") #only need bcSeq, so sampleInfo.csv is ok here
  
  I1 = I1[as.vector(sread(I1)) %in% metadata$bcSeq]
  samples = metadata[match(as.character(sread(I1)), metadata$bcSeq), "alias"]
  
  #only necessary if using native data - can parse out description w/ python
  I1Names =  sapply(strsplit(as.character(ShortRead::id(I1)), " "), "[[", 1)#for some reason we can't dynamically set name/id on ShortRead!
  
  rm(I1)
  
  suppressWarnings(dir.create("Data/demultiplexedReps"))
  
  #R1 
  R1 = readFastq("Data/Undetermined_S0_L001_R1_001.fastq.gz")
  R1Names = sapply(strsplit(as.character(ShortRead::id(R1)), " "), "[[", 1)#for some reason we can't dynamically set name/id on ShortRead!
  names(R1Names) = NULL
  
  R1 = R1[match(I1Names, R1Names)]
  R1 = split(R1, samples)
  for (i in 1:length(R1)){writeFastq(R1[[i]], paste0("Data/demultiplexedReps/&", names(R1[i]), "_S0_L001_R1_001.fastq.gz"), mode="w")}
  
  rm(R1, R1Names)
  
  #R2
  R2 = readFastq("Data/Undetermined_S0_L001_R2_001.fastq.gz")
  R2Names = sapply(strsplit(as.character(ShortRead::id(R2)), " "), "[[", 1) #for some reason we can't dynamically set name/id on ShortRead!
  names(R2Names) = NULL
  
  R2 = R2[match(I1Names, R2Names)]
  R2 = split(R2, samples)
  
  for (i in 1:length(R2)){writeFastq(R2[[i]], paste0("Data/demultiplexedReps/&", names(R2[i]), "_S0_L001_R2_001.fastq.gz"), mode="w")}
  
  rm(R2, R2Names, I1Names, samples) 
}

errorCorrectBC = function(){
  library("ShortRead")
  
  metadata = read.csv("sampleInfo.csv") #only number of samples, so one file is ok here
  
  bushmanJobID = get(load("bushmanJobID.RData"))
  
  I1 = readFastq("Data/Undetermined_S0_L001_I1_001.fastq.gz")
  I1 = trimTailw(I1, 2, "0", 12)
  I1 = I1[width(I1)==max(width(I1))]
  
  I1 = split(I1, ceiling(seq_along(I1)/500000))
  
  for(chunk in names(I1)){
    writeFasta(I1[[chunk]], file=paste0("Data/trimmedI1-", chunk, ".fasta"))
  }
  
  bsub(jobName=paste0("BushmanErrorCorrectWorker_", bushmanJobID, "[1-", length(I1),"]"),
       logFile="logs/errorCorrectWorkerOutput%I.txt",
       command="python ~/EAS/PMACS_scripts/processGolay.py"
  )
  
  bsub(queue="max_mem64",
       wait=paste0("done(BushmanErrorCorrectWorker_", bushmanJobID, ")"),
       jobName=paste0("BushmanDemultiplex_", bushmanJobID),
       logFile="logs/demultiplexOutput.txt",
       command=paste0("Rscript -e \"source('~/EAS/PMACS_scripts/programFlow.R'); demultiplex();\"")
  )
  
  #trim seqs
  bsub(queue="plus",
       wait=paste0("done(BushmanDemultiplex_", bushmanJobID, ")"),
       jobName=paste0("BushmanTrimSeqs_", bushmanJobID, "[1-", nrow(metadata), "]"),
       logFile="logs/trimOutput%I.txt",
       command=paste0("Rscript -e \"source('~/EAS/PMACS_scripts/programFlow.R'); trimSeqs();\"")
  )
  
  system("sleep 5") #allow time for jobs to enter queue
  
  #post-trim processing, also kicks off alignment and int site calling jobs
  bsub(wait=paste0("done(BushmanTrimSeqs_", bushmanJobID, ")"),
       jobName=paste0("BushmanPostTrimProcessing_", bushmanJobID),
       logFile="logs/postTrimOutput.txt",
       command=paste0("Rscript -e \"source('~/EAS/PMACS_scripts/programFlow.R'); postTrimSeqs();\"")    
  )
}

postTrimSeqs = function(){
  #sample index in callIntSites_PMACS.R is based on parameters.RData, so we
  #should pass appropriate index values by using parameters.RData here
  metadata = get(load("parameters.RData"))
  
  bushmanJobID = get(load("bushmanJobID.RData"))
  blatStartPort = get(load("bushmanBlatStartPort.RData"))
  
  numAliases  = nrow(metadata)
  
  numFastaFiles = length(system("ls */*.fa", intern=T))
  
  #align seqs
  bsub(wait=paste0("done(BushmanPostTrimProcessing_", bushmanJobID, ")"),
       jobName=paste0("BushmanAlignSeqs_", bushmanJobID, "[", blatStartPort, "-", blatStartPort+numFastaFiles-1, "]"),
       logFile="logs/alignOutput%I.txt",
       command=paste0("Rscript -e \"source('~/EAS/PMACS_scripts/programFlow.R'); alignSeqs();\"")
  )
  
  system("sleep 5")
  #call int sites (have to find out which ones worked)
  successfulTrims = unlist(sapply(metadata, function(x){
    get(load(paste0(x[["alias"]], "/trimStatus.RData"))) == paste0(getwd(), "/", x[["alias"]])
  }))
  
  bsub(queue="max_mem30",
       wait=paste0("done(BushmanAlignSeqs_", bushmanJobID, ")"),
       jobName=paste0("BushmanCallIntSites_", bushmanJobID, "[", paste(which(successfulTrims), collapse=","), "]"),
       logFile="logs/callSitesOutput%I.txt",
       command=paste0("Rscript -e \"source('~/EAS/PMACS_scripts/programFlow.R'); callIntSites();\"")
  )
  
  bsub(wait=paste0("done(BushmanCallIntSites_", bushmanJobID, ")"),
       jobName=paste0("BushmanCleanup_", bushmanJobID),
       logFile="logs/cleanupOutput.txt",
       command=paste0("Rscript -e \"source('~/EAS/PMACS_scripts/programFlow.R'); cleanup();\"")
  )
}

trimSeqs = function(){
  source('~/EAS/PMACS_scripts/pairedIntSites.R')
  
  sampleID = as.integer(system("echo $LSB_JOBINDEX", intern=T))
  
  parameters = get(load("parameters.RData"))[[sampleID]]
  
  alias = parameters[["alias"]]
  
  suppressWarnings(dir.create(alias, recursive=TRUE))
  
  status = tryCatch(eval(as.call(append(getTrimmedSeqs, unname(parameters[c("qualityThreshold", "badQualityBases",
                                                                            "qualitySlidingWindow", "primer", "ltrBit",
                                                                            "largeLTRFrag", "linkerSequence", "linkerCommon",
                                                                            "mingDNA", "read1", "read2", "alias", "vectorSeq")])))), 
                    error=function(e){print(paste0("Caught error: ", e$message))})
  
  save(status, file="trimStatus.RData") #working directory is changed while executing getTrimmedSeqs
}