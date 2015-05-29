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


linker_input <- c(
    "CGCCAACGTAGCGTTGCACCTCCGCTTAAGGGACT",
    "GACGCGTAAGGTCGGGTAGCTCCGCTTAAGGGACT" )

expected_linker_common_output <- c("AGTCCCTTAAGCGGA", "AGTCCCTTAAGCGGA")

test_that("Linker with no N", {
    expect_equal(expected_linker_common_output, linker_common(linker_input) )
})


linker_input <- c(
    "CGCCAACGTAGCGTTGCACNNNNNNNNNNNNCTCCGCTTAAGGGACT",
    "GACGCGTAAGGTCGGGTAGNNNNNNNNNNNNCTCCGCTTAAGGGACT",
    "CGCCAACGTAGCGTTGCACCTCCGCTTAAGGGACT",
    "GACGCGTAAGGTCGGGTAGCTCCGCTTAAGGGACT" )

expected_linker_common_output <- c(
    "AGTCCCTTAAGCGGAG",
    "AGTCCCTTAAGCGGAG",
    "AGTCCCTTAAGCGGA",
    "AGTCCCTTAAGCGGA")

test_that("Mixed kinds of linkers", {
    expect_equal(expected_linker_common_output, linker_common(linker_input) )
})


linker_input <- c(
    "CGCCAACGTAGCGNNNTTGCACCTCNCGCTTAAGGGACT",
    "GACGCGTAAGGTCGGGTAGCTCCGCTTAAGGGACT" )

test_that("Should fail if seperated N's in the middle", {
    expect_error(linker_common(linker_input))
})

