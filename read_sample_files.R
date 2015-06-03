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
    cbind(sampleInfo, processingParams)
}

.check_input <- function(sample_info_path, processing_info_path) {
    stopifnot(file.exists(sample_info_path) & file.exists(processing_info_path))
}
