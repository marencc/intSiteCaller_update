# intSiteCaller

***


## Introduction
This code is designed to take fastq files that are produced by the MiSeq and return integration sites, multihits, and chimeras.  RData and fasta files are used for data persistance in place of a database, thus allowing massive parallelization and the LSF job submission system on the [PMACS HPC](http://www.med.upenn.edu/hpc/hardware-physical-environment.html)


## Inputs
                                    
Analysis is started by having the user create the following directory structure:

```
primaryAnalysisDirectory
├── Data
│   ├── Undetermined_S0_L001_I1_001.fastq.gz
│   ├── Undetermined_S0_L001_R1_001.fastq.gz
│   └── Undetermined_S0_L001_R2_001.fastq.gz
├── processingParams.csv
└── sampleInfo.csv
``` 
#### Primary Analysis Directory

* `Data/Undetermined_S0_L001_*_001.fastq.gz` are the fastq files returned by the MiSeq (R1, R2, and I1)
    
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

## Usage

After creating the above directory structure, the following command is issued:

```Rscript intSiteCaller.R```

The rest of the processing is fully automated and shouldn't take more than 4 hours to process 1.5e7 raw reads.

`intSiteCaller.R` can handle the following options
* `-j`, `--jobID` - Unique name by which to identify this intance of intSiteCaller [default: intSiteCallerJob]
* `-c`, `--codeDir` - Directory where intSiteCaller code is stored, can be relative or absolute [default: .]
* `-p`, `--primaryAnalysisDir` - Location of primary analysis directory, can be relative or absolute [default: .]
* `-C`, `--cleanup` - Remove temporary files upon successful execution of intSiteCaller
* `-h`, `--help` - Show the help message and exit



## Outputs

#### Integration sites
This code returns integration sites in two formats.  `allSites.RData` is a `GRanges` object that contains a single record for each Illumina read.  `sites.final.RData` is a `GRanges` object of dereplicated integration sites along with a tally of how many reads were seen for each site (irregardless of sonic breakpoint).  The `revmap` column in `sites.final.RData` links unique reads from `allSites` to dereplicated sites in `sites.final`.


#### Multihits
Multihits are stored in `multihitData.RData` which is a `list`.  The first item in this list is a `GRanges` object where each record represents a properly-paired alignment.  Individual multihit reads can be identified by analysing the `ID` column, which cooresponds to the unique Illumina read identifier.  The second item in `multihitData` is a `list` of `list`s.  Each object in the primary list represents a multihit cluster.  The two sub-objects are:
* `potentialSites` is a `GRanges` object where each record is a soloStart representation of each potential integration site of a given multihit cluster.
* `widths` is a `data.frame` where the `Var1` column is the width of the multihit cluster and the `Freq` column is the number of reads where that width was seen for the given multihit cluster.  This is analagous to PCR breakpoints and counts for unique sites.

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
* `argparse`
* Any `BSgenome.*.UCSC.*` package cooresponding to reference genomes specified in `processingParams.csv`

Specific versioning analysis has not yet been performed.

Additionally, BLAT code requires the availability of the `blat` and `python` command.  `blat` is available on PMACS at `/opt/software/blatSrc/v35/bin/x86_64`, however is not included in the default `PATH`.  This directory will need to be added to the user's default `PATH` in order to successfullly run BLAT alignments. 

`intSiteCaller` confirms the presence of all dependancies and will throw an error if a dependancy is not met.

## Code Structure

- Primary read trimming and integration site calling logic is contained in `intSiteLogic.R`.
- Branching and condensing of PMACS jobs is performed in `programFlow.R`
- Barcode error correcting logic is performed in `errorCorrectIndices/golay.py` as wrapped by `errorCorrectIndices/processGolay.py`.
- All code files are checked into the repository.
- Flowcharts will be added to graphically describe the flow of the overall program as well as the flow/logic of individual functions


## Tests

A sample dataset is included for verification of integration site calling accuracy.  The `testCases` directory contains a subdirectory, `intSiteValidation`, and a compressed folder, `intSiteValidationOUTPUT.tar.gz`.  To analyze the test data, move `intSiteValidation/*` to an LSF envrionment and execute the code as described in the 'Usage' section.  Upon successful execution, the directory should resemble `intSiteValidationOUTPUT.tar.gz`.  Note that this subset of data contains samples with some borderline cases.  For example, clone7 samples should all fail, and many of the clone1-clone4 samples should return no multihits or chimeras.  The current implementation of the code handles these gracefully.