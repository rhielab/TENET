% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/step7TopGenesUserPeakOverlap.R
\name{step7TopGenesUserPeakOverlap}
\alias{step7TopGenesUserPeakOverlap}
\title{Identify if RE DNA methylation sites linked to top genes/transcription
factors are located within a specific distance of specified genomic regions}
\usage{
step7TopGenesUserPeakOverlap(
  TENETMultiAssayExperiment,
  peakData,
  geneAnnotationDataset = NA,
  DNAMethylationArray = NA,
  hypermethGplusAnalysis = TRUE,
  hypomethGplusAnalysis = TRUE,
  topGeneNumber = 10,
  distanceFromREDNAMethylationSites = 100,
  coreCount = 1
)
}
\arguments{
\item{TENETMultiAssayExperiment}{Specify a MultiAssayExperiment object
containing expression and methylation SummarizedExperiment objects, such as
one created by the TCGADownloader function. This MultiAssayExperiment object
should also contain the results from the \code{step5OptimizeLinks} and
\code{step6DNAMethylationSitesPerGeneTabulation} functions in its metadata.}

\item{peakData}{Specify a data frame, matrix, or GRanges object with
genomic regions (peaks) of interest, organized in a bed-like manner (see
\url{https://genome.ucsc.edu/FAQ/FAQformat.html#format1}), or a path to a
directory containing .bed, .narrowPeak, .broadPeak, and/or .gappedPeak files
with peaks of interest. The files may optionally be compressed
(.gz/.bz2/.xz).}

\item{geneAnnotationDataset}{Specify a gene annotation dataset which is
used to identify names for genes by their Ensembl IDs. The argument must be
either a GRanges object (such as one imported via \code{rtracklayer::import}) or a
path to a GFF3 or GTF file. Both GENCODE and Ensembl annotations are
supported. Other annotation datasets may work, but have not been tested.
See the "Input data" section of the vignette for information on the required
dataset format.
Specify NA to use the names for genes listed in the "geneName" column of the
elementMetadata of the rowRanges of the "expression" SummarizedExperiment
object within the TENETMultiAssayExperiment object. Defaults to NA.}

\item{DNAMethylationArray}{Specify the name of a DNA methylation probe
array supported by the sesameData package (see
\code{?sesameData::sesameData_getManifestGRanges}). If an array is specified,
RE DNA methylation sites and their locations in that array's manifest are
cross-referenced with RE DNA methylation site IDs included in the rownames
of the methylation dataset provided in the "methylation"
SummarizedExperiment object within the TENETMultiAssayExperiment object, and
only those overlapping will be considered for analysis. If set to NA, all RE
DNA methylation sites with locations listed in the rowRanges of the
"methylation" SummarizedExperiment object are used. Defaults to NA.}

\item{hypermethGplusAnalysis}{Set to TRUE to create data frames with the
peak overlap information for the unique hypermethylated RE DNA methylation
sites linked to the top genes/TFs by most hypermethylated RE DNA methylation
sites with G+ links. Defaults to TRUE.}

\item{hypomethGplusAnalysis}{Set to TRUE to create data frames with the peak
overlap information for the unique hypomethylated RE DNA methylation sites
linked to the top genes/TFs by most hypomethylated RE DNA methylation sites
with G+ links. Defaults to TRUE.}

\item{topGeneNumber}{Specify the number of top genes/TFs, based on the most
linked RE DNA methylation sites of a given analysis type, for which to
generate data showing overlap with user peak files for the RE DNA
methylation sites linked to those genes. Defaults to 10.}

\item{distanceFromREDNAMethylationSites}{Specify either 0 or a positive
integer in base pairs to be the distance from the linked RE DNA methylation
sites within which an RE DNA methylation site will be considered to overlap
a peak. Defaults to 100.}

\item{coreCount}{Argument passed as the mc.cores argument to mclapply. See
\code{?parallel::mclapply} for more details. Defaults to 1.}
}
\value{
Returns the MultiAssayExperiment object given as the
TENETMultiAssayExperiment argument with an additional list of information
named 'step7TopGenesUserPeakOverlap' in its metadata with the output of this
function. This list is subdivided into hypermethGplus or hypomethGplus
results as selected by the user, which are further subdivided into data
frames with data for the unique RE DNA methylation sites linked to the top
overall genes, and for top TF genes only. Each of these data frames contain
a row for each of the unique RE DNA methylation sites linked to the top
genes/TFs for the specified analysis types. These note the locations of the
RE DNA methylation sites and the specified search windows, whether the RE
DNA methylation site is linked to each of the top genes/TFs, and whether
each RE DNA methylation site was found within the specified distance to any
peak in each of the user's peak files.
}
\description{
This function takes the top genes/transcription factors (TFs) by number of
linked RE DNA methylation sites identified by the
\code{step6DNAMethylationSitesPerGeneTabulation} function, up to the number
specified by the user, and identifies if the RE DNA methylation sites linked
to those genes/TFs from the hyper- or hypomethylated G+ analysis quadrants
are found in the vicinity of genomic regions (peaks) of interest, supplied
by the user in the form of .bed, .narrowPeak, .broadPeak, and/or gappedPeak
files in a specified directory or given as a data frame or GRanges object.
}
\examples{
\dontshow{if (interactive()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
## This example uses the example MultiAssayExperiment provided in the
## TENET.ExperimentHub package to do peak overlapping for all unique RE DNA
## methylation sites linked to the top 10 genes, by number of linked hyper-
## or hypomethylated RE DNA methylation sites, using a GRanges object
## containing the genomic coordinates of peaks of interest. Gene names, and
## the locations of RE DNA methylation sites, will be retrieved from the
## rowRanges of the 'expression' and 'methylation' SummarizedExperiment
## objects in the example MultiAssayExperiment. The analysis will be
## performed using one CPU core.

## Load the example TENET MultiAssayExperiment object
## from the TENET.ExperimentHub package
exampleTENETMultiAssayExperiment <-
    TENET.ExperimentHub::exampleTENETMultiAssayExperiment()

## Load the example peak GRanges object from the TENET.ExperimentHub package
exampleTENETPeakRegions <- TENET.ExperimentHub::exampleTENETPeakRegions()

## Use the example datasets to do the peak overlapping
returnValue <- step7TopGenesUserPeakOverlap(
    TENETMultiAssayExperiment = exampleTENETMultiAssayExperiment,
    peakData = exampleTENETPeakRegions
)

## This example uses the example MultiAssayExperiment provided in the
## TENET.ExperimentHub package to do peak overlapping for all unique RE DNA
## methylation sites linked to the top 5 genes only, by number of linked
## hypomethylated RE DNA methylation sites only, using bed-like files
## containing the genomic coordinates of peaks of interest in the user's R
## working directory. Gene names will be retrieved from the rowRanges of the
## 'expression' SummarizedExperiment object in the example
## MultiAssayExperiment, and RE DNA methylation sites and their locations
## retrieved from the HM450 array via the sesameData package. A window of
## 500 base pairs will be used to identify if the RE DNA methylation sites
## lie within the vicinity of peaks. The analysis will be performed using 8
## CPU cores.

## Load the example TENET MultiAssayExperiment object from the
## TENET.ExperimentHub package
exampleTENETMultiAssayExperiment <-
    TENET.ExperimentHub::exampleTENETMultiAssayExperiment()

## Use the example dataset to do the peak overlapping
returnValue <- step7TopGenesUserPeakOverlap(
    TENETMultiAssayExperiment = exampleTENETMultiAssayExperiment,
    peakData = ".",
    DNAMethylationArray = "HM450",
    hypermethGplusAnalysis = FALSE,
    topGeneNumber = 5,
    distanceFromREDNAMethylationSites = 500,
    coreCount = 8
)
\dontshow{\}) # examplesIf}
}
