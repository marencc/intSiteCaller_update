#    This source code file is a component of the larger INSPIIRED genomic analysis software package.
#    Copyright (C) 2016 Frederic Bushman
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

#' read sample information and processing information
#' and form one large dataframe.
#' @param processing_param_file the same parameters for all samples
#' @param sample_info_file information about sample and techincal DNA used.
#' @return df with cols:
#'  alias, primer, ltrBit, largeLTRFrag, linkerSequence, 
#'  linkerCommon, bcSeq, vectorSeq, gender
#'  qualityThreshold, badQualityBases, qualitySlidingWindow,
#'  maxAlignStart, maxFragLength, minPctIdent, mingDNA
read_sample_processing_files <- function(sample_info_path, processing_info_path) {
    .check_input(sample_info_path, processing_info_path)
    sampleInfo <- read.delim(sample_info_path, stringsAsFactors=F)
    sampleInfo$linkerCommon <- linker_common(sampleInfo$linkerSequence)
    processingParams <- read.delim(processing_info_path, stringsAsFactors=F)
    complete_meta <- cbind(sampleInfo, processingParams)
    .check_column_names(complete_meta)
    .check_sex_is_correct(complete_meta$gender)
    .check_sampleName_is_correct(complete_meta$alias)
    complete_meta
}

.check_input <- function(sample_info_path, processing_info_path) {
    stopifnot(file.exists(sample_info_path) & file.exists(processing_info_path))
}

.columns_sample_proc <- c( "qualityThreshold", "badQualityBases", "qualitySlidingWindow",
    "primer", "ltrBit", "largeLTRFrag", "linkerSequence", "linkerCommon",
    "mingDNA", "alias", "vectorSeq", "minPctIdent",
    "maxAlignStart", "maxFragLength", "gender", "bcSeq", "refGenome"
)

.check_column_names <- function(samples) {
    stopifnot(all(.columns_sample_proc %in% names(samples)))
}

.valid_sex <- c('m', 'f')

.check_sex_is_correct <- function(gender) {
    stopifnot(all(gender %in% .valid_sex))
}

.check_sampleName_is_correct <- function(samples) {
    stopifnot( ! (is.null(samples)))
    stopifnot(length(samples) == length(unique(samples)))
}
