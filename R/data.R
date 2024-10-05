#' @title Human transcription factor list
#'
#' @description A character vector of gene Ensembl IDs which were identified as
#' human TFs by Lambert SA et al (PMID: 29425488). Candidate proteins were
#' manually examined by a panel of experts based on available data. Proteins
#' with experimentally demonstrated DNA binding specificity were considered TFs.
#' Other proteins, such as co-factors and RNA binding proteins, were classified
#' as non-TFs. **Citation:** Lambert SA, Jolma A, Campitelli LF, et al. The
#' Human Transcription Factors. Cell. 2018 Feb 8;172(4):650-665. doi:
#' 10.1016/j.cell.2018.01.029. Erratum in: Cell. 2018 Oct 4;175(2):598-599.
#' PMID: 29425488.
#'
#' @usage data("humanTranscriptionFactorList", package = "TENET")
#'
#' @docType data
#'
#' @format A character vector containing 1,639 Ensembl IDs of known human TFs.
#'
#' @examples
#' data("humanTranscriptionFactorList", package = "TENET")
#'
#' @source \url{http://humantfs.ccbr.utoronto.ca/download.php}
"humanTranscriptionFactorList"
