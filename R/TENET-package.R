## Note: The TENET.ExperimentHub importFrom exists to silence an R CMD check
## warning. TENET.ExperimentHub is listed in Imports because it is required
## by function examples.

## Note: The sesame importFrom exists to silence an R CMD check warning. sesame
## is listed in Imports because it is required indirectly by TCGAbiolinks to
## download methylation data.

#' @keywords internal
#' @title TENET (Tracing regulatory Element Networks using Epigenetic Traits)
#' @description TENET identifies key transcription factors and regulatory
#' elements linked to a specific cell type by finding significantly correlated
#' differences in gene expression and regulatory element methylation between
#' case and control input datasets, and identifying the top genes by number of
#' significant RE DNA methylation site links. It also includes many additional
#' tools to aid in the visualization and analysis of the results, including
#' plots displaying and comparing methylation and expression data and RE DNA
#' methylation site link counts, survival analysis, TF motif searching in the
#' vicinity of linked RE DNA methylation sites, custom TAD and peak overlap
#' analysis, and UCSC Genome Browser track file generation. A utility function
#' is also provided to download methylation, expression, and patient survival
#' data from The Cancer Genome Atlas (TCGA) for use in TENET or other analyses.
#' @docType package
#' @name TENET
#' @importFrom TENET.ExperimentHub exampleTENETTADRegions
#' @importFrom sesame openSesame
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
