library(testthat)
source("../read_sample_files.R")
source("../linker_common.R")

sample_file <- "../testCases/intSiteValidation/sampleInfo.tsv"
processing_file <- "../testCases/intSiteValidation/processingParams.tsv"
completeMetadata <- read_sample_processing_files(sample_file, processing_file)

sample_columns <- c("alias", "linkerSequence", "bcSeq", "gender",
    "primer", "ltrBit", "largeLTRFrag", "linkerCommon", "vectorSeq")

processing_columns <- c( "qualityThreshold", "badQualityBases", "qualitySlidingWindow",
    "mingDNA", "minPctIdent", "maxAlignStart", "maxFragLength", "refGenome")

all_columns <- c(sample_columns, processing_columns)

test_that("Can merge sample info and proc params", {
    expect_equal(nrow(completeMetadata), 41)
    expect_equal(ncol(completeMetadata), length(all_columns))
})

test_that("all columns are present", {
    expect_true(all(c(
        "qualityThreshold", "badQualityBases", "qualitySlidingWindow",
        "primer", "ltrBit", "largeLTRFrag", "linkerSequence", "linkerCommon",
        "mingDNA", "alias", "vectorSeq", "minPctIdent",
        "maxAlignStart", "maxFragLength", "gender") %in% names(completeMetadata))
    ) 
})

test_that("number of samples stays the same after merging", {
    expect_equal(nrow(completeMetadata), nrow(read.delim(sample_file)))
})

test_that("all processing parameters are the same", {
    sapply(processing_columns, function(x)
        expect_equal(length(unique(completeMetadata[x])), 1)
    )
})

test_that("Should fail if there are no files", {
    expect_error(read_sample_processing_files("this file does not exist", processing_file))
    expect_error(read_sample_processing_files(sample_file, "this file does not exist neither"))
})

