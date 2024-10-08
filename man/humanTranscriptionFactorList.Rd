% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{humanTranscriptionFactorList}
\alias{humanTranscriptionFactorList}
\title{Human transcription factor list}
\format{
A character vector containing 1,639 Ensembl IDs of known human TFs.
}
\source{
\url{http://humantfs.ccbr.utoronto.ca/download.php}
}
\usage{
data("humanTranscriptionFactorList", package = "TENET")
}
\description{
A character vector of gene Ensembl IDs which were identified as
human TFs by Lambert SA et al (PMID: 29425488). Candidate proteins were
manually examined by a panel of experts based on available data. Proteins
with experimentally demonstrated DNA binding specificity were considered TFs.
Other proteins, such as co-factors and RNA binding proteins, were classified
as non-TFs. \strong{Citation:} Lambert SA, Jolma A, Campitelli LF, et al. The
Human Transcription Factors. Cell. 2018 Feb 8;172(4):650-665. doi:
10.1016/j.cell.2018.01.029. Erratum in: Cell. 2018 Oct 4;175(2):598-599.
PMID: 29425488.
}
\examples{
data("humanTranscriptionFactorList", package = "TENET")

}
\keyword{datasets}
