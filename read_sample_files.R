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
    return data.frame(alias=c('a', 'b'), gender=c('F', 'M'))
}
