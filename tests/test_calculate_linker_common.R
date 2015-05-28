library(testthat)
source("../linker_common.R")

linker_input <- c(
    "CGCCAACGTAGCGTTGCACNNNNNNNNNNNNCTCCGCTTAAGGGACT",
    "GACGCGTAAGGTCGGGTAGNNNNNNNNNNNNCTCCGCTTAAGGGACT"
    )

expected_linker_common_output <- c("AGTCCCTTAAGCGGAG", "AGTCCCTTAAGCGGAG")

test_that("Second sequence after N", {
    expect_equal(expected_linker_common_output, linker_common(linker_input) )
})


linker_noN_input <- c(
    "CGCCAACGTAGCGTTGCACCTCCGCTTAAGGGACT",
    "GACGCGTAAGGTCGGGTAGCTCCGCTTAAGGGACT"
    )

test_that("Fail if no N in the middle", {
    expect_error(linker_common(linker_noN_input))
})

