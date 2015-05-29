library(testthat)
source("../linker_common.R")

linker_input <- c(
    "CGCCAACGTAGCGTTGCACNNNNNNNNNNNNCTCCGCTTAAGGGACT",
    "GACGCGTAAGGTCGGGTAGNNNNNNNNNNNNCTCCGCTTAAGGGACT" )
    

expected_linker_common_output <- c("AGTCCCTTAAGCGGAG", "AGTCCCTTAAGCGGAG")

test_that("Reverse complement of the sequence after N", {
    expect_equal(expected_linker_common_output, linker_common(linker_input) )
})

test_that("Reverse complement of the sequence after N, vector of length 1", {
    expect_equal(expected_linker_common_output[1], linker_common(linker_input[1]) )
})


linker_noN_input <- c(
    "CGCCAACGTAGCGTTGCACCTCCGCTTAAGGGACT",
    "GACGCGTAAGGTCGGGTAGCTCCGCTTAAGGGACT" )

test_that("Fail if no N in the middle", {
    expect_error(linker_common(linker_noN_input))
})

linker_moreN_input <- c(
    "CGCCAACGTAGCGNNNTTGCACCTCNCGCTTAAGGGACT",
    "GACGCGTAAGGTCGGGTAGCTCCGCTTAAGGGACT" )

test_that("Fail if more N in the middle", {
    expect_error(linker_common(linker_noN_input))
})

