primerIDAlignSeqs = function(subjectSeqs=NULL, patternSeq=NULL, 
                             qualityThreshold1=0.75, qualityThreshold2=0.50, 
                             doAnchored=FALSE, doRC=TRUE, returnUnmatched=FALSE, 
                             returnRejected=FALSE, showStats=FALSE, ...) {
  
  .checkArgs_SEQed()
  
  qualityThreshold1 <- as.numeric(qualityThreshold1)
  qualityThreshold2 <- as.numeric(qualityThreshold2)
  
  ## make sure there are Ns in the patternSeq for considering primerIDs
  if(length(unlist(gregexpr("N",patternSeq)))<4) {
    stop("There should be minimum of atleast 4 Ns in patternSeq to be ",
         "considered as a primerID sequence.")
  }
  
  ## Get the right orientation of the supplied patternSeq to peform proper search 
  ## at later two step search phase. 
  if(doRC) {
    patternSeq <- tryCatch(doRCtest(subjectSeqs, patternSeq),
                           error=function(e) patternSeq)
  }
  
  primerIDpos <- unlist(gregexpr("N", patternSeq))
  
  ## Perform primerID extraction by breaking the pattern into two parts...
  ## for sanity sakes due to homopolymers ##
  pattern1 <- as.character(subseq(DNAString(patternSeq), 1, primerIDpos[1]-1))
  pattern2 <- as.character(subseq(DNAString(patternSeq), 
                                  primerIDpos[length(primerIDpos)]+1))
  
  pattern1.hits <- pairwiseAlignSeqs(subjectSeqs, pattern1, "middle", 
                                     qualityThreshold=qualityThreshold1, 
                                     doRC=FALSE, returnUnmatched=TRUE, ...)
  pattern2.hits <- pairwiseAlignSeqs(subjectSeqs, pattern2, "middle", 
                                     qualityThreshold=qualityThreshold2,
                                     doRC=FALSE, ...)
  
  ## Set aside reads which has no match to the pattern1...
  ## most likely mispriming if primerID is on 5' or read was too loong if on 3'
  ## Hits returned from pairwiseAlignSeqs will be filtered for low scored hits...
  ## so no need to check for those from subjectSeqs
  unmatched <- pattern1.hits[["Absent"]]
  pattern1.hits <- pattern1.hits[["hits"]]
  
  ## remove reads which only have a match to one of the patterns...crossover most likely!
  rejected1 <- setdiff(names(pattern1.hits), names(pattern2.hits))
  rejected2 <- setdiff(names(pattern2.hits), names(pattern1.hits))
  rejected <- c(pattern1.hits[names(pattern1.hits) %in% rejected1], 
                pattern2.hits[names(pattern2.hits) %in% rejected2])
  rm("rejected1","rejected2")
  if(showStats) { 
    message("Removed ", length(rejected),
            " read(s) for only matching one of pattern1 or pattern2") 
  }
  
  ## use only reads which match to both sides of the patterns.
  good.rows <- intersect(names(pattern1.hits), names(pattern2.hits))
  pattern1.hits <- pattern1.hits[names(pattern1.hits) %in% good.rows]
  pattern2.hits <- pattern2.hits[names(pattern2.hits) %in% good.rows]
  
  stopifnot(identical(names(pattern1.hits), names(pattern2.hits)))
  stopifnot(identical(names(pattern1.hits), good.rows))
  
  ## make sure there is no overlap of ranges between pattern1.hits & 
  ## pattern2.hits if there is, then no primerID was found...remove it
  badAss <- end(pattern1.hits) >= start(pattern2.hits)
  if(any(badAss)) {
    rejected <- c(rejected, pattern1.hits[badAss], pattern2.hits[badAss])
    pattern1.hits <- pattern1.hits[!badAss]
    pattern2.hits <- pattern2.hits[!badAss]
    good.rows <- good.rows[!badAss]
    bore <- table(badAss)["TRUE"]
    message("Removed ", ifelse(is.na(bore),0,bore),
            " read(s) for not having primerID present between pattern 1 & 2")
  }
  
  hits <- IRanges(start=start(pattern1.hits), 
                  end=end(pattern2.hits), 
                  names=good.rows)        
  
  primerIDs <- IRanges(start=end(pattern1.hits)+1, 
                       end=start(pattern2.hits)-1, 
                       names=good.rows)        
  
  if(length(hits)==0) {
    stop("No hits found that matched both sides of patternSeq with ",
         "primerID in the middle.")
  }    
  
  ## do anchored search for only sequences that matched both sides of patternSeq
  if(doAnchored) {
    message("Found ", length(primerIDs), 
            " total primerIDs before anchored filter.")
    
    ## get anchors of bases flanking Ns
    anchorBase.s <- substr(patternSeq, primerIDpos[1]-1, primerIDpos[1]-1)
    anchorBase.e <- substr(patternSeq, primerIDpos[length(primerIDpos)]+1, 
                           primerIDpos[length(primerIDpos)]+1)
    
    anchorBase.s.i <- trimSeqs(subjectSeqs, 
                               resize(pattern1.hits,width=1,fix="end"))
    anchorBase.e.i <- trimSeqs(subjectSeqs, 
                               resize(pattern2.hits,width=1,fix="start"))
    rows <- anchorBase.s==as.character(anchorBase.s.i) & 
      anchorBase.e==as.character(anchorBase.e.i)
    
    unAnchored <- hits[!rows]
    primerIDs <- primerIDs[rows]
    hits <- hits[rows]
    message("Found ", length(primerIDs), 
            " total primerIDs after anchored filter.")
    rm("rows","anchorBase.s.i","anchorBase.e.i")
  }
  rm("good.rows","pattern1.hits","pattern2.hits")
  cleanit <- gc()
  
  ## also remove any primerIDs that were too short or big than intended.
  nSize <- 2
  badAss <- !width(primerIDs) %in% 
    (length(primerIDpos)-nSize):(length(primerIDpos)+nSize)
  if(any(badAss)) {
    rejectedprimerIDs <- hits[badAss]
    hits <- hits[!badAss]
    primerIDs <- primerIDs[!badAss]
    message("Removed ", table(badAss)["TRUE"],
            " read(s) for not having right primerID length")
  }
  
  hits <- IRangesList("hits"=hits, "primerIDs"=primerIDs)
  
  if(exists("unAnchored")) {
    if(length(unAnchored)>0) { 
      hits <- append(hits, IRangesList("unAnchoredprimerIDs"=unAnchored)) 
    }
  }
  
  if(returnUnmatched & length(unmatched)>0) {
    hits <- append(hits, IRangesList("Absent"=unmatched))
  }
  
  if(returnRejected) {
    if(length(rejected)>0) { 
      hits <- append(hits, IRangesList("Rejected"=rejected)) 
    }
    if(exists("rejectedprimerIDs")) { 
      if(length(rejectedprimerIDs)>0) { 
        hits <- append(hits, IRangesList("RejectedprimerIDs"=rejectedprimerIDs)) 
      } 
    }
  }    
  return(hits)
}


trimSeqs = function(dnaSet, coords, side="middle", offBy=0) {
  stopifnot(class(dnaSet) %in% c("DNAStringSet", "DNAString"))
  stopifnot(class(coords)=="IRanges")
  
  if(length(dnaSet)==0 | length(coords)==0) {
    stop("dnaSet/coords is empty. Please supply reads/coords to be trimmed.")
  }
  
  # check if both dnaSet and coords has 'names' attribute, 
  # if yes then check if they have matching names, else check lengths. 
  if(is.null(names(dnaSet)) | is.null(names(coords))) {
    stopifnot(length(dnaSet)==length(coords))
  } else {
    rows <- match(names(coords), names(dnaSet))
    if(any(is.na(rows))) {
      stop("Some of the names in coords are not present in dnaSet")
    }
    if(!is.ordered(rows)) {
      dnaSet <- dnaSet[rows]
      if(!identical(names(dnaSet), names(coords))) {
        stop("Names are not identical between dnaSet and coords parameter")
      }
    }
  }        
  
  # temp helper function to show messages #
  .showMessage <- function(x) {
    message("Following sequences were removed from trimming since their ",
            "coordinates+offBy were out of sequence length: ", 
            paste(x,collapse=", "))
  }
  
  # trim by side and check if any of the coords are off the sequence length in dnaSet
  if(tolower(side)=="left") {
    test <- end(coords)+offBy > width(dnaSet) | end(coords)+offBy < 1
    if(any(test)) {
      .showMessage(names(dnaSet)[test])
    }
    subseq(dnaSet[!test], start=end(coords[!test])+offBy)
  } else if (tolower(side)=="right") {
    test <- start(coords)-offBy > width(dnaSet) | end(coords)+offBy < 1
    if(any(test)) {
      .showMessage(names(dnaSet)[test])
    }
    subseq(dnaSet[!test], end=start(coords[!test])-offBy)
  } else {
    test <- start(coords)+offBy > width(dnaSet) | 
      end(coords)-offBy > width(dnaSet) | start(coords)+offBy < 1
    if(any(test)) {
      .showMessage(names(dnaSet)[test])
    }
    subseq(dnaSet[!test], 
           start=start(coords[!test])+offBy, 
           end=end(coords[!test])-offBy)
  }    
}


pairwiseAlignSeqs = function(subjectSeqs=NULL, patternSeq=NULL, side="left", 
                             qualityThreshold=1, showStats=FALSE, bufferBases=5, 
                             doRC=TRUE, returnUnmatched=FALSE, 
                             returnLowScored=FALSE, parallel=FALSE, ...) {
  dp <- NULL
  
  .checkArgs_SEQed()
  
  if(parallel) {
    subjectSeqs2 <- chunkize(subjectSeqs)
    hits <- bplapply(subjectSeqs2, function(x) 
      pairwiseAlignSeqs(x, patternSeq, side, qualityThreshold, showStats=FALSE, 
                        bufferBases, doRC, returnUnmatched,  returnLowScored, 
                        parallel=FALSE, ...), BPPARAM=dp)    
    hits <- do.call(c, hits)
    if(is(hits,"CompressedIRangesList")) {
      attrs <- unique(names(hits))
      hits <- sapply(attrs, 
                     function(x) unlist(hits[names(hits)==x],use.names=FALSE))
      IRangesList(hits)
    } else {
      hits
    }
  } else {
    qualityThreshold <- as.numeric(qualityThreshold)
    
    ## only get the relevant side of subject sequence with extra bufferBases to 
    ## account for indels/mismatches & save memory while searching and avoid 
    ## searching elsewhere in the sequence
    if(tolower(side)=="left") {
      badSeqs <- DNAStringSet()
      culprits <- width(subjectSeqs) < (nchar(patternSeq)+bufferBases)
      if(any(culprits)) {
        badSeqs <- subjectSeqs[culprits]
        message(length(badSeqs),
                " sequences were removed from aligning since they were",
                " shorter than pattern getting aligned: ",
                (nchar(patternSeq)+bufferBases),"bp")            
        subjectSeqs <- subjectSeqs[!culprits]            
      }
      subjectSeqs2 <- subseq(subjectSeqs, start=1, 
                             end=(nchar(patternSeq)+bufferBases))
      overFromLeft <- rep(0,length(subjectSeqs))
    } else if (tolower(side)=="right") { 
      overFromLeft <- width(subjectSeqs)-(nchar(patternSeq)+bufferBases)
      overFromLeft[overFromLeft<1] <- 1
      subjectSeqs2 <- subseq(subjectSeqs, start=overFromLeft)
    } else {
      subjectSeqs2 <- subjectSeqs
      overFromLeft <- rep(0, length(subjectSeqs))
    }
    
    ## search both ways to test which side yields more hits!
    if(doRC) {
      patternSeq <- tryCatch(doRCtest(subjectSeqs2, patternSeq, 
                                      qualityThreshold),
                             error=function(e) patternSeq)
    }
    
    ## type=overlap is best for primer trimming...see Biostrings Alignment vignette
    if(any(names(match.call()) %in% c("type","gapOpening","gapExtension"))) {
      hits <- pairwiseAlignment(subjectSeqs2, patternSeq, ...)        
    } else {
      hits <- pairwiseAlignment(subjectSeqs2, patternSeq, type="overlap", 
                                gapOpening=-1, gapExtension=-1, ...)
    }
    stopifnot(length(hits)==length(subjectSeqs2))
    
    scores <- round(score(hits))
    highscored <- scores >= round(nchar(patternSeq)*qualityThreshold)*2
    
    if(!any(highscored)) {
      stop("No hits found which passed the qualityThreshold")
    }
    
    ## basically a small subset of highscored
    unmatched <- nchar(hits) <= round(nchar(patternSeq)*.1) 
    
    # no point in showing stats if all sequences are a potential match! #
    if(showStats & qualityThreshold!=0) {
      bore <- as.numeric(table(highscored)['FALSE'])
      message("Total of ",ifelse(is.na(bore),0,bore),
              " did not have the defined pattern sequence (", patternSeq,
              ") that passed qualityThreshold on the ", side, " side")
    }
    
    ## extract starts-stops of the entire pattern hit ##
    starts <- start(pattern(hits))
    ends <- end(pattern(hits))
    namesq <- names(subjectSeqs)
    hits <- IRanges(start=starts+overFromLeft-ifelse(side=="right",2,0), 
                    end=ends+overFromLeft-ifelse(side=="right",2,0), 
                    names=namesq)
    rm("scores","subjectSeqs2","subjectSeqs","starts","ends","namesq")
    
    ## no need to test if there were any multiple hits since pairwiseAlignment 
    ## will only output one optimal alignment...see the man page.
    if(!returnLowScored & !returnUnmatched) {      
      hits <- hits[highscored]
    } else {
      hitstoreturn <- IRangesList("hits"=hits[highscored])
      if(returnLowScored & length(hits[!highscored])>0) {
        hitstoreturn <- append(hitstoreturn, 
                               IRangesList("Rejected"=hits[!highscored]))
      }
      
      if(returnUnmatched & length(hits[unmatched])>0) {
        hitstoreturn <- append(hitstoreturn, 
                               IRangesList("Absent"=hits[unmatched]))
      }
      hits <- hitstoreturn
      rm(hitstoreturn)
    }             
    cleanit <- gc()
    return(hits)
  }
}


function(query=NULL, subject=NULL, standaloneBlat=TRUE, port=5560, 
         host="localhost", parallel=TRUE, numServers=1L,
         gzipResults=TRUE,
         blatParameters=c(minIdentity=90, minScore=10, stepSize=5, 
                          tileSize=10, repMatch=112312, dots=50, 
                          maxDnaHits=10, q="dna", t="dna", 
                          out="psl")) {
  
  if(length(system("which blat",intern = TRUE))==0) { 
    stop("Command blat not found!")
  }
  
  ## get all BLAT options from the system for comparison to blatParameters later
  suppressWarnings(blatOpts <- unique(sub("\\s+-(\\S+)=\\S+.+", "\\1", 
                                          grep("\\s+-.+=", 
                                               system("blat",intern=TRUE), 
                                               value=TRUE))))
  
  suppressWarnings(gfClientOpts <- unique(sub("\\s+-(\\S+)=\\S+.+", "\\1", 
                                              grep("\\s+-.+=", 
                                                   system("gfClient", 
                                                          intern=TRUE), 
                                                   value=TRUE))))
  
  gfServerOpts <- c("tileSize","stepSize","minMatch","maxGap", "trans","log",
                    "seqLog","syslog","logFacility","mask","repMatch",
                    "maxDnaHits", "maxTransHits","maxNtSize","maxAsSize",
                    "canStop")                                               
  
  if(!standaloneBlat) { 
    message("Using gfClient protocol to perform BLAT.")
    if(is.null(port)) {
      stop("The port paramter is empty. ",
           "Please define the port used to start gfServer with")
    }
  }
  
  ## check the subject parameter
  if(is.null(subject) | length(subject)==0) {
    stop("The subject parameter is empty. ", 
         "Please supply subject sequences or a path to 2bit or nib files to ",
         "serve as reference/target")
  } else {
    subjectFile <- NULL
    if(is.atomic(subject)) {
      if (any(grepl("\\.2bit$|\\.nib$", subject, ignore.case=TRUE))) {
        if(standaloneBlat) { 
          stop("Standalone BLAT cannot be used when subject is an indexed ",
               "nib or 2bit file.") 
        }
        indexFileDir <- dirname(subject)
        subjectFile <- list.files(path=indexFileDir, 
                                  pattern=basename(subject), full.names=TRUE)
        if(length(subjectFile)==0) { 
          stop("The file(s) supplied in subject parameter doesn't exist.") }
      } else {
        ## change object type if necessary for troubleshooting purpose in later 
        ## steps
        subject <- DNAStringSet(subject)
      }
    }
    
    if(is.null(subjectFile)) {
      ## subjectFile is still null so it means that subject is a DNAStringSet
      if(is.null(names(subject))) { ## add names of subject if not present
        names(subject) <- paste("subject", 1:length(subject))
      }
      
      ## write out the subject sequences into a fasta file
      filename.seq <- paste("subjectFile.fa",runif(1),"tempyS",sep=".")      
      writeXStringSet(subject, filepath=filename.seq, format="fasta")                                  
      subjectFile <- filename.seq
    }
  }
  
  ## check the query parameter
  if(is.null(query) | length(query)==0) {
    stop("The query parameter is empty. Please supply reads to be aligned")
  } else {
    queryFiles <- NULL
    if(is.atomic(query)) {
      if (any(grepl("\\.fna$|\\.fa$|\\.fastq$|\\.fasta$|\\*", query, 
                    ignore.case=TRUE))) {
        ## detect whether query paramter is a regex or list of files
        if(any(grepl("\\*|\\$|\\+|\\^",query))) {
          queryFiles <- list.files(path=dirname(query), pattern=basename(query), 
                                   full.names=TRUE)            
        } else {
          queryFiles <- query
        }
        
        if(parallel) {
          ## split the fasta files into smaller chunks for parallel BLATing
          queryFiles <- unlist(sapply(queryFiles,
                                      function(f) 
                                        splitSeqsToFiles(f, bpworkers(),
                                                         "tempyQ")), 
                               use.names=FALSE)                    
        }
      } else {
        ## change object type if necessary for troubleshooting purpose in later steps
        query <- DNAStringSet(query)
      }
    }
    
    if(is.null(queryFiles)) {
      ## queryFiles is still null so it means that query is a DNAStringSet           
      if(is.null(names(query))) {  ## fix names of query if not present
        names(query) <- paste("read", 1:length(query),sep="-")
      }  
      
      ## write out the query sequences into fasta files
      if(parallel) {
        queryFiles <- splitSeqsToFiles(query, bpworkers(), "tempyQ")
      } else {
        queryFiles <- paste("queryFile.fa",runif(1),"tempyQ",sep=".")          
        writeXStringSet(query, filepath=queryFiles, format="fasta")                
      }
    }
  }
  
  ## perform the Blatting of queryFiles vs subjectFile/indexFiles using 
  ## gfClient/standalone BLAT  
  
  ## do some formatting ##
  queryFiles <- as.character(queryFiles)
  subjectFile <- as.character(subjectFile)
  
  dp <- if(parallel){ bpparam() } else { SerialParam() }
  
  ## BLAT it ##
  if(standaloneBlat) {        
    blatOpts <- blatParameters[names(blatParameters) %in% blatOpts]
    stopifnot(length(subjectFile)==1)
    filenames <- bplapply(queryFiles, function(x) {
      filename.out <- paste(x, blatOpts["out"], sep=".")
      cmd <- paste("blat", paste(paste0("-",names(blatOpts)), blatOpts, 
                                 collapse=" ", sep="="), "-noHead", 
                   subjectFile, x, filename.out)
      message(cmd)
      system(cmd)
      
      ## no need to save splitted files!
      if(grepl("\\.tempyQ$",x)) { system(sprintf("rm %s",x)) } 
      
      if(gzipResults) { 
        system(paste("gzip", filename.out))
        filename.out <- paste(filename.out, "gz", sep=".") 
      }
      filename.out
    }, BPPARAM=dp)
    
    if(grepl("\\.tempyS$",subjectFile)) { system(sprintf("rm %s",subjectFile)) }
    
  } else {
    # start the gfServer if not started already! #
    killFlag <- FALSE
    port <- port + 0:(numServers-1)
    for(n in 1:numServers) {
      searchCMD <- sprintf("gfServer status %s %s", host, port[n])
      if(system(searchCMD,ignore.stderr=TRUE)!=0) {
        message(sprintf("Starting gfServer # %s.", n))
        startgfServer(seqDir=subjectFile, host=host, port=port[n], 
                      gfServerOpts=blatParameters[names(blatParameters) 
                                                  %in% gfServerOpts])
        killFlag <- TRUE
      }
    }
    
    gfClientOpts <- blatParameters[names(blatParameters) %in% gfClientOpts]
    stopifnot(length(subjectFile)>0)
    filenames <- bpmapply(function(port, x) {
      filename.out <- paste(x, gfClientOpts["out"], sep=".")
      cmd <- paste("gfClient",
                   paste(paste0("-", names(gfClientOpts)), gfClientOpts, 
                         collapse=" ", sep="="), "-nohead", host, port, "/", x, 
                   filename.out)
      message(cmd)
      system(cmd)
      
      ## no need to save splitted files!
      if(grepl("\\.tempyQ$", x)) { file.remove(x) } 
      
      if(gzipResults) { 
        system(paste("gzip", filename.out))
        filename.out <- paste(filename.out, "gz", sep=".") 
      }
      filename.out      
    }, rep(port, length=length(queryFiles)), queryFiles, 
    SIMPLIFY=FALSE, USE.NAMES=FALSE, BPPARAM=dp)  
    
    ## only kill if the gfServer was started from within this function ## 
    if(killFlag) {
      # stop to conserve memory #
      message("Kill gfServer.")        
      sapply(port, function(x) stopgfServer(port=x))
    }
  }
  return(unlist(filenames))
}



.checkArgs_SEQed = .checkArgs_SEQed <- function() {
  
  checks <- expression( 
    if("subjectSeqs" %in% names(formals())) { 
      if(is.null(subjectSeqs) | length(subjectSeqs)==0) {
        stop("subjectSeqs paramter is empty. Please supply reads to be aligned")
      }      
      
      ## give names if not there for troubleshooting purpose in later steps
      removeSubjectNamesAfter <- FALSE
      if(is.null(names(subjectSeqs))) {
        removeSubjectNamesAfter <- TRUE
        names(subjectSeqs) <- paste("read", 1:length(subjectSeqs))
      } 
    },
    
    if("patternSeq" %in% names(formals())) { 
      if(is.null(patternSeq) | length(patternSeq)==0) {
        stop("patternSeq paramter is empty. Please supply reads to be aligned")
      } else if (length(patternSeq)>1) {
        stop("More than 1 patternSeq defined. Please only supply one pattern.")
      }
    },
    
    if("reads" %in% names(formals())) { 
      if(is.null(reads) | length(reads)==0) {
        stop("reads paramter is empty. Please supply reads to be aligned")
      }
    },
    
    if("Vector" %in% names(formals())) { 
      if(is.null(Vector) | length(Vector)==0) {
        stop("Vector paramter is empty. Please supply Vector to be aligned")
      }
      
      if(!grepl("DNAString",class(Vector))) {
        stop("Vector paramter is not of DNAString class")
      }
    },
    
    if("parallel" %in% names(formals())) { 
      dp <- if(parallel) { bpparam() } else { SerialParam() }
    },
    
    if("sampleInfo" %in% names(formals())) { 
      stopifnot(is(sampleInfo,"SimpleList"))
    },
    
    if("feature" %in% names(formals())) {
      if(is.null(feature)) {
        stop("Please define a feature to extract.")
      }
    }      
  )
  
  eval.parent(checks)
}

read.psl = function(pslFile=NULL, bestScoring=TRUE, asGRanges=FALSE, 
                    removeFile=TRUE, parallel=FALSE) {
  qName <- dp <- NULL
  files <- pslFile
  .checkArgsSetDefaults_ALIGNed()
  
  ## setup psl columns + classes
  cols <- pslCols()
  
  hits <- bplapply(files, function(x) {
    message(x)
    ## add extra fields incase pslx format ##
    ncol <- max(count.fields(x, sep = "\t"))
    if(ncol > length(cols)) {
      for(f in 1:(ncol-length(cols))) {
        cols[paste0("V",f)] <- "character"
      }
    }
    hits.temp <- read.delim(x, header=FALSE, col.names=names(cols), 
                            stringsAsFactors=FALSE, colClasses=cols)    
    if(bestScoring) {  
      ## do round one of bestScore here to reduce file size          
      hits.temp$score <- with(hits.temp, 
                              matches-misMatches-qBaseInsert-tBaseInsert)
      isBest <- with(hits.temp, ave(score, qName, FUN=function(x) x==max(x)))
      hits.temp <- hits.temp[as.logical(isBest),]
      rm("isBest")
    }
    hits.temp    
  }, BPPARAM=dp)  
  hits <- unique(rbind.fill(hits))
  
  if(nrow(hits)==0) {
    if(removeFile) { file.remove(pslFile) }
    stop("No hits found")
  }
  
  ## do round two of bestScore incase any got missed in round one
  if(bestScoring) {
    message("\t cherry picking!")
    hits$score <- with(hits, matches-misMatches-qBaseInsert-tBaseInsert)    
    isBest <- with(hits, ave(score, qName, FUN=function(x) x==max(x)))
    hits <- hits[as.logical(isBest),]
    rm("isBest")
  }
  
  if(asGRanges) {
    hits <- pslToRangedObject(hits, useTargetAsRef=TRUE)
  }
  
  message("Ordering by qName")
  if(is(hits,"GRanges")) {
    hits <- hits[sort(hits$qName, index.return=T)$ix] #hacked together because version dependancies suck
  } else {
    hits <- arrange(hits, qName)
  }  
  
  if(removeFile) { file.remove(pslFile) }
  
  return(hits)
}


blatSeqs = function(query=NULL, subject=NULL, standaloneBlat=TRUE, port=5560, 
                    host="localhost", parallel=TRUE, numServers=1L,
                    gzipResults=TRUE,
                    blatParameters=c(minIdentity=90, minScore=10, stepSize=5, 
                                     tileSize=10, repMatch=112312, dots=50, 
                                     maxDnaHits=10, q="dna", t="dna", 
                                     out="psl")) {
  
  if(length(system("which blat",intern = TRUE))==0) { 
    stop("Command blat not found!")
  }
  
  ## get all BLAT options from the system for comparison to blatParameters later
  suppressWarnings(blatOpts <- unique(sub("\\s+-(\\S+)=\\S+.+", "\\1", 
                                          grep("\\s+-.+=", 
                                               system("blat",intern=TRUE), 
                                               value=TRUE))))
  
  suppressWarnings(gfClientOpts <- unique(sub("\\s+-(\\S+)=\\S+.+", "\\1", 
                                              grep("\\s+-.+=", 
                                                   system("gfClient", 
                                                          intern=TRUE), 
                                                   value=TRUE))))
  
  gfServerOpts <- c("tileSize","stepSize","minMatch","maxGap", "trans","log",
                    "seqLog","syslog","logFacility","mask","repMatch",
                    "maxDnaHits", "maxTransHits","maxNtSize","maxAsSize",
                    "canStop")                                               
  
  if(!standaloneBlat) { 
    message("Using gfClient protocol to perform BLAT.")
    if(is.null(port)) {
      stop("The port paramter is empty. ",
           "Please define the port used to start gfServer with")
    }
  }
  
  ## check the subject parameter
  if(is.null(subject) | length(subject)==0) {
    stop("The subject parameter is empty. ", 
         "Please supply subject sequences or a path to 2bit or nib files to ",
         "serve as reference/target")
  } else {
    subjectFile <- NULL
    if(is.atomic(subject)) {
      if (any(grepl("\\.2bit$|\\.nib$", subject, ignore.case=TRUE))) {
        if(standaloneBlat) { 
          stop("Standalone BLAT cannot be used when subject is an indexed ",
               "nib or 2bit file.") 
        }
        indexFileDir <- dirname(subject)
        subjectFile <- list.files(path=indexFileDir, 
                                  pattern=basename(subject), full.names=TRUE)
        if(length(subjectFile)==0) { 
          stop("The file(s) supplied in subject parameter doesn't exist.") }
      } else {
        ## change object type if necessary for troubleshooting purpose in later 
        ## steps
        subject <- DNAStringSet(subject)
      }
    }
    
    if(is.null(subjectFile)) {
      ## subjectFile is still null so it means that subject is a DNAStringSet
      if(is.null(names(subject))) { ## add names of subject if not present
        names(subject) <- paste("subject", 1:length(subject))
      }
      
      ## write out the subject sequences into a fasta file
      filename.seq <- paste("subjectFile.fa",runif(1),"tempyS",sep=".")      
      writeXStringSet(subject, filepath=filename.seq, format="fasta")                                  
      subjectFile <- filename.seq
    }
  }
  
  ## check the query parameter
  if(is.null(query) | length(query)==0) {
    stop("The query parameter is empty. Please supply reads to be aligned")
  } else {
    queryFiles <- NULL
    if(is.atomic(query)) {
      if (any(grepl("\\.fna$|\\.fa$|\\.fastq$|\\.fasta$|\\*", query, 
                    ignore.case=TRUE))) {
        ## detect whether query paramter is a regex or list of files
        if(any(grepl("\\*|\\$|\\+|\\^",query))) {
          queryFiles <- list.files(path=dirname(query), pattern=basename(query), 
                                   full.names=TRUE)            
        } else {
          queryFiles <- query
        }
        
        if(parallel) {
          ## split the fasta files into smaller chunks for parallel BLATing
          queryFiles <- unlist(sapply(queryFiles,
                                      function(f) 
                                        splitSeqsToFiles(f, bpworkers(),
                                                         "tempyQ")), 
                               use.names=FALSE)                    
        }
      } else {
        ## change object type if necessary for troubleshooting purpose in later steps
        query <- DNAStringSet(query)
      }
    }
    
    if(is.null(queryFiles)) {
      ## queryFiles is still null so it means that query is a DNAStringSet           
      if(is.null(names(query))) {  ## fix names of query if not present
        names(query) <- paste("read", 1:length(query),sep="-")
      }  
      
      ## write out the query sequences into fasta files
      if(parallel) {
        queryFiles <- splitSeqsToFiles(query, bpworkers(), "tempyQ")
      } else {
        queryFiles <- paste("queryFile.fa",runif(1),"tempyQ",sep=".")          
        writeXStringSet(query, filepath=queryFiles, format="fasta")                
      }
    }
  }
  
  ## perform the Blatting of queryFiles vs subjectFile/indexFiles using 
  ## gfClient/standalone BLAT  
  
  ## do some formatting ##
  queryFiles <- as.character(queryFiles)
  subjectFile <- as.character(subjectFile)
  
  dp <- if(parallel){ bpparam() } else { SerialParam() }
  
  ## BLAT it ##
  if(standaloneBlat) {        
    blatOpts <- blatParameters[names(blatParameters) %in% blatOpts]
    stopifnot(length(subjectFile)==1)
    filenames <- bplapply(queryFiles, function(x) {
      filename.out <- paste(x, blatOpts["out"], sep=".")
      cmd <- paste("blat", paste(paste0("-",names(blatOpts)), blatOpts, 
                                 collapse=" ", sep="="), "-noHead", 
                   subjectFile, x, filename.out)
      message(cmd)
      system(cmd)
      
      ## no need to save splitted files!
      if(grepl("\\.tempyQ$",x)) { system(sprintf("rm %s",x)) } 
      
      if(gzipResults) { 
        system(paste("gzip", filename.out))
        filename.out <- paste(filename.out, "gz", sep=".") 
      }
      filename.out
    }, BPPARAM=dp)
    
    if(grepl("\\.tempyS$",subjectFile)) { system(sprintf("rm %s",subjectFile)) }
    
  } else {
    # start the gfServer if not started already! #
    killFlag <- FALSE
    port <- port + 0:(numServers-1)
    for(n in 1:numServers) {
      searchCMD <- sprintf("gfServer status %s %s", host, port[n])
      if(system(searchCMD,ignore.stderr=TRUE)!=0) {
        message(sprintf("Starting gfServer # %s.", n))
        startgfServer(seqDir=subjectFile, host=host, port=port[n], 
                      gfServerOpts=blatParameters[names(blatParameters) 
                                                  %in% gfServerOpts])
        killFlag <- TRUE
      }
    }
    
    gfClientOpts <- blatParameters[names(blatParameters) %in% gfClientOpts]
    stopifnot(length(subjectFile)>0)
    filenames <- bpmapply(function(port, x) {
      filename.out <- paste(x, gfClientOpts["out"], sep=".")
      cmd <- paste("gfClient",
                   paste(paste0("-", names(gfClientOpts)), gfClientOpts, 
                         collapse=" ", sep="="), "-nohead", host, port, "/", x, 
                   filename.out)
      message(cmd)
      system(cmd)
      
      ## no need to save splitted files!
      if(grepl("\\.tempyQ$", x)) { file.remove(x) } 
      
      if(gzipResults) { 
        system(paste("gzip", filename.out))
        filename.out <- paste(filename.out, "gz", sep=".") 
      }
      filename.out      
    }, rep(port, length=length(queryFiles)), queryFiles, 
    SIMPLIFY=FALSE, USE.NAMES=FALSE, BPPARAM=dp)  
    
    ## only kill if the gfServer was started from within this function ## 
    if(killFlag) {
      # stop to conserve memory #
      message("Kill gfServer.")        
      sapply(port, function(x) stopgfServer(port=x))
    }
  }
  return(unlist(filenames))
}

splitSeqsToFiles = function(x, totalFiles=4, suffix="tempy", 
                            filename="queryFile.fa") {
  if(is.atomic(x)) {
    message("Splitting file ",x)
    totalSeqs <- length(fasta.info(x, use.names=FALSE))
    chunks <- round(totalSeqs/totalFiles)
    ## incase totalSeqs is lower than number of files to be created!
    chunks <- ifelse(chunks>0, chunks, totalSeqs) 
    
    starts <- seq(0, totalSeqs, by=chunks) ## create chunks of starts    
    for(skippy in starts[starts!=totalSeqs]) {
      filename.out <- paste(x, skippy, runif(1), suffix, sep=".")
      ## no need to read the entire file...save memory by reading in N lines
      query.tmp <- readBStringSet(x,nrec=chunks, skip=skippy) 
      writeXStringSet(query.tmp, filepath=filename.out, format="fasta")            
    }
    return(list.files(path=dirname(x), 
                      pattern=paste0(basename(x),".*", suffix, "$"), 
                      full.names=TRUE))
  } else if (class(x)=="DNAStringSet") {
    message("Splitting Reads.")
    totalSeqs <- length(x)
    chunks <- round(totalSeqs/totalFiles)
    starts <- seq(1, totalSeqs, by=chunks)
    stops <- unique(c(seq(chunks, totalSeqs, by=chunks), totalSeqs))
    stopifnot(length(starts)==length(stops))        
    for(skippy in 1:length(starts)) {
      filename.out <- paste(filename, skippy, runif(1), suffix, sep=".")            
      writeXStringSet(x[starts[skippy]:stops[skippy]], filepath=filename.out,
                      format="fasta")            
    }            
    return(list.files(path=".", 
                      pattern=paste(filename,".*",suffix,"$",sep=""), 
                      full.names=TRUE))
  } else {
    stop("Dont know what is supplied in parameter x.")
  }
}

.checkArgsSetDefaults_ALIGNed <- function() {
  
  checks <- expression(
    if("pslFile" %in% names(formals()) | "files" %in% names(formals())) { 
      if(is.null(files) | length(files)==0) {
        stop("files parameter empty. Please supply a filename to be read.")
      }
      
      if(any(grepl("\\*|\\$|\\+|\\^",files))) {
        ## vector of filenames
        files <- list.files(path=dirname(files), pattern=basename(files), 
                            full.names=TRUE)      
      }
      
      if(length(files)==0) { 
        stop("No file(s) found with given paramter in files:", files) 
      }
    },
    
    if("psl.rd" %in% names(formals())) {
      if(!is.null(psl.rd)) {
        if(length(psl.rd)==0 | !is(psl.rd,"GRanges")) {
          stop("psl.rd paramter is empty or not a GRanges object")
        }
      }
    },
    
    if("parallel" %in% names(formals())) { 
      dp <- if(parallel) { bpparam() } else { SerialParam() }
    }
  )
  
  eval.parent(checks)
}

pslCols <- function(withClass=TRUE) {
  cols <- c("matches", "misMatches", "repMatches", "nCount", "qNumInsert", 
            "qBaseInsert", "tNumInsert", "tBaseInsert", "strand", "qName", 
            "qSize", "qStart", "qEnd", "tName", "tSize", "tStart", "tEnd", 
            "blockCount", "blockSizes", "qStarts", "tStarts")
  
  cols.class <- c(rep("numeric",8), rep("character",2), rep("numeric",3),
                  "character", rep("numeric",4), rep("character",3))
  
  if(withClass) {
    structure(cols.class, names=cols)
  } else {
    cols
  }
}

pslToRangedObject <- function(x, useTargetAsRef=TRUE, isblast8=FALSE) {
  if(useTargetAsRef) {
    metadataCols <- c(setdiff(names(x), c("tName","tStart","tEnd","strand")),
                      ifelse(isblast8, NA, "tStarts"))
    out <- GRanges(seqnames=x$tName, IRanges(start=x$tStart, end=x$tEnd),
                   strand=x$strand)     
  } else {
    metadataCols <- c(setdiff(names(x), c("qName","qStart","qEnd","strand")),
                      ifelse(isblast8, NA, "qStarts"))
    out <- GRanges(seqnames=x$qName, IRanges(start=x$qStart, end=x$qEnd),
                   strand=x$strand)
  }
  
  for(f in na.omit(metadataCols)) {
    mcols(out)[[f]] <- x[,f]
  } 
  
  out
}