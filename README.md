# intSiteCaller

***


## Introduction
This code is designed to take fastq files that are produced by the MiSeq and return integration sites, multihits, and chimeras.  RData and fasta files are used for data persistance in place of a database, thus allowing massive parallelization and the LSF job submission system on the [PMACS HPC](http://www.med.upenn.edu/hpc/hardware-physical-environment.html)


## Inputs and Usage
                                    
Analysis is started by having the user create the following directory structure:

```
primaryAnalysisDirectory
├── Data
│   ├── Undetermined_S0_L001_I1_001.fastq.gz
│   ├── Undetermined_S0_L001_R1_001.fastq.gz
│   └── Undetermined_S0_L001_R2_001.fastq.gz
├── PMACS_kickoff.R
├── processingParams.csv
└── sampleInfo.csv
``` 
#### Primary Analysis Directory

* `Data/Undetermined_S0_L001_*_001.fastq.gz` are the fastq files returned by the MiSeq (R1, R2, and I1)

* `PMACS_kickoff.R` contains additional metadata parameters that must be set by the user.  These are located in the first four lines of the file.
  * `bushmanJobID` is 
	* `blatStartPort`
	* `codeDir` is the directory where the main code logic files (`programFlow.R`, `intSiteLogic.R`, and `BLATsamples.sh`) are located
	* `cleanup`	is `TRUE` if temporary files are to be removed from the primary analysis directory, `FALSE` otherwise
    
* `processingParams.csv` contains 'dryside' sample metadata:
	* `alias` is the human-readable sample description
	* `qualityThreshold, badQualityBases, qualitySlidingWindow`
		* After error-correcting and demultiplexing, intSiteCaller trims raw MiSeq reads based on Illumina-issued quality scores.  `badQualityBases` is the number of bases below `qualityThreshold`, using [standard Illumina ASCII Q-score encoding](http://support.illumina.com/content/dam/illumina-support/documents/myillumina/a557afc4-bf0e-4dad-9e59-9c740dd1e751/casava_userguide_15011196d.pdf) (p.41-2), that can be observed in a window of `qualitySlidingWindow` before the read is trimmed.
	* `primer` is the primer sequence as seen in MiSeq read 2
	* `ltrBit` is the LTR sequence as seen in MiSeq read 2
	* `largeLTRFrag` is 43nt of the LTR sequence as seen from MiSeq read **1**
	* `linkerCommon` is 15nt of the linker sequence as seen from MiSeq read **2**
	* `mingDNA` is the minimum length of genomic DNA (in nt) that is allowed to be passed onto the alignment step
	* `vectorSeq` is a filepath (either absolute or relative to the *primary analysis directory*) to the vector sequence in fasta format -- it is encouraged to place the vector sequence directly in the primary analysis directory, although that is not a requirement
	* `minPctIdent` is the minimum percent identity that a query sequence has to a putative alignment on the target sequence in order to be considered 'valid'
	* `maxAlignStart` is the maximum number of nucleotides of the query sequence that do *not* match the putative target sequence before a matched nucleotide is seen
	* `maxFragLength` is the maximum length of a properly-paired alignment (in nt)
	* `refGenome` is the reference genome to be used for alignment - this is passed in as a standard text string (ex. 'hg18', 'hg19', 'mm8', 'mm9')

* `sampleInfo.csv` contains 'wetside' sample metadata:
	* `alias` is the human-readable sample description
	* `linkerSequence` is the linker sequence as seen in MiSeq read 1.  N's indicate the presence of a primerID
	* `bcSeq` is the barcode used during sample preparation 
	* `gender` is either 'M' of 'F' for male/female, respectively


After creating the directory structure, the following command is issued from within the primary analysis directory:

`bsub -n1 -q normal -J "BushmanKickoff_analysis" -o logs/kickoffOutput.txt Rscript PMACS_kickoff.R`

The rest of the processing is fully automated and shouldn't take more than 4 hours to process 1.5e7 raw reads.

Temporary files and logs are currently stored by default.  Theses files can be automatically removed post-run by adjusting the `cleanup` variable in `PMACS_kickoff.R.`
                                


## Outputs

#### Integration sites
This code returns integration sites in two formats.  `allSites.RData` is a `GRanges` object that contains a single record for each Illumina read.  `sites.final.RData` is a `GRanges` object of dereplicated integration sites along with a tally of how many reads were seen for each site (irregardless of sonic breakpoint).  The `revmap` column in `sites.final.RData` links unique reads from `allSites` to dereplicated sites in `sites.final`.


#### Multihits
Multihits are stored in `multihitData.RData` which is a `GRangesList`.  The first item in this list is a `GRanges` object where each record represents a properly-paired alignment.  Individual multihit reads can be identified by analysing the `ID` column, which cooresponds to the unique Illumina read identifier.  The second item in `multihitData` is a `GRanges` object of dereplicated multihits, which lists each unique genomic integration site as a unique record.  The `revmap` column pairs records from `multihitData[[1]]` to `multihitData[[2]]`.  The third item is a `GRanges` object of multihit clusters.  This is still in development.

#### Chimeras
Chimeras are stored in `chimeraData.RData` which is a list that contains some basic chimera frequency statistics and a `GRangesList` object.  Each `GRanges` object contains two records, one for the read1 alignment and another for the read2 alignment

#### PrimerIDs
PrimerIDs (if present in the linker sequence) are stored in `primerIDData.RData`.  This file is a base R `list` containing a `DNAStringSet` and a `BStringSet` containing the sequences and quality scores of the primerIDs.

#### Stats
Processing statistics are returned in the `stats.RData` file.  This file contains a single `data.frame`.  Detail of the specific columns provided in this dataframe will be added later and can inferred from `intSiteLogic.R`.



## Dependencies

This code is highly dependent on Bioconductor packages for processing DNA data and collapsing/expanding alignments.

The following R packages and their subsesequent dependencies are required for proper operation of `intSiteCaller`:
* `ShortRead`
* `hiReadsProcessor`
* `GenomicRanges`
* `rtracklayer`
* `BSgenome`
* Any `BSgenome.*.UCSC.*` package cooresponding to reference genomes specified in `processingParams.csv`

Specific versioning analysis has not yet been performed.  The `sessionInfo()` from an environment successfully running this code is included for reference below until specific versioning requirements can be identified.

```
other attached packages:
 [1] ShortRead_1.24.0        GenomicAlignments_1.2.1 Rsamtools_1.18.2       
 [4] BiocParallel_1.0.0      hiReadsProcessor_0.1.5  xlsx_0.5.7             
 [7] xlsxjars_0.6.1          rJava_0.9-6             hiAnnotator_1.0        
[10] BSgenome_1.34.0         Biostrings_2.34.1       XVector_0.6.0          
[13] plyr_1.8.1              RMySQL_0.9-3            DBI_0.3.1              
[16] rtracklayer_1.26.2      iterators_1.0.7         foreach_1.4.2          
[19] GenomicRanges_1.18.3    GenomeInfoDb_1.2.3      IRanges_2.0.1          
[22] S4Vectors_0.4.0         BiocGenerics_0.12.1     sonicLength_1.4.4   
```


## Code Structure

- Primary read trimming and integration site calling logic is contained in `intSiteLogic.R`.
- Branching and condensing of PMACS jobs is performed in `programFlow.R`
- Barcode error correcting logic is performed in `errorCorrectIndices/golay.py` as wrapped by `errorCorrectIndices/processGolay.py` and alignment parameters are contained in `BLATsamples.sh`.
- All core code files are checked into the repository.  Pointing `codeDir` in `PMACS_kickoff.R` to the cloned repository directory is sufficient to run `intSiteCaller`.
- Flowcharts will be added to graphically describe the flow of the overall program as well as the flow/logic of individual functions


## Tests

A sample dataset is included for verification of integration site calling accuracy.  The `testCases` directory contains a subdirectory, `intSiteValidation`, and a compressed folder, `intSiteValidationOUTPUT.tar.gz`.  To analyze the test data, move `intSiteValidation/*` to an LSF envrionment, make any necessary changes to `PMACS_kickoff.R`, and execute the code as described in the 'Usage' section.  Upon successful execution, the directory should resemble `intSiteValidationOUTPUT.tar.gz`.  Note that this subset of data contains samples with some borderline cases.  For example, clone7 samples should all fail, and many of the clone1-clone4 samples should return no multihits or chimeras.  The current implementation of the code handles these gracefully.