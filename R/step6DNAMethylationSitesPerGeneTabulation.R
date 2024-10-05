## Internal functions used by step 6

## Internal function to perform tabulation on a given data quadrant
.tabulateQuadrant <- function(
    resultsDF, geneIDNameDF, substepDescription) {
    ## Check that there are results present from step 5.
    ## .ensureStepPresent checks that the run was performed, but there may
    ## still not be any links found in step 5. If there were not, warn the user
    ## and return NULL for this quadrant.
    if (nrow(resultsDF) == 0) {
        .warningNoCall(
            "No significant ", substepDescription, " were identified in ",
            "step 5. Please check and rerun the step3GetAnalysisZScores, ",
            "step4SelectMostSignificantLinksPerDNAMethylationSite, and/or ",
            "step5OptimizeLinks functions. No results will be generated for ",
            "this quadrant."
        )

        return(NULL)
    }

    ## Create a table dataset from the input data frame
    quadrantAllGeneDataset <- as.data.frame(
        table(resultsDF$geneID),
        stringsAsFactors = FALSE
    )

    ## Rename the columns
    colnames(quadrantAllGeneDataset) <- c("geneID", "count")

    ## Add the gene names to the dataset
    quadrantAllGeneDataset$geneName <- geneIDNameDF[
        quadrantAllGeneDataset$geneID, "geneName"
    ]

    ## Sort the data frame by decreasing number of links
    quadrantAllGeneDataset <- quadrantAllGeneDataset[
        order(quadrantAllGeneDataset$count, decreasing = TRUE),
    ]

    ## Set the row names to the gene IDs
    rownames(quadrantAllGeneDataset) <- quadrantAllGeneDataset$geneID

    ## Return the dataset
    return(quadrantAllGeneDataset)
}

## Main step 6 function

#' Tabulate the total number of RE DNA methylation sites linked to each of the
#' genes
#'
#' This function takes the final optimized regulatory element
#' RE DNA methylation site-gene links identified by the `step5OptimizeLinks`
#' function and tabulates the number of these links per gene. This tabulation is
#' done separately for both of the hyper- or hypomethylated G+ analysis
#' quadrants, as selected by the user.
#'
#' @param TENETMultiAssayExperiment Specify a MultiAssayExperiment object
#' containing expression and methylation SummarizedExperiment objects, such
#' as one created by the TCGADownloader function. This MultiAssayExperiment
#' object should also contain the results from the `step5OptimizeLinks`
#' function in its metadata.
#' @param geneAnnotationDataset Specify a gene annotation dataset which is
#' used to identify names for genes by their Ensembl IDs. The argument must be
#' either a GRanges object (such as one imported via `rtracklayer::import`) or a
#' path to a GFF3 or GTF file. Both GENCODE and Ensembl annotations are
#' supported. Other annotation datasets may work, but have not been tested.
#' See the "Input data" section of the vignette for information on the required
#' dataset format.
#' Specify NA to use the names for genes listed in the "geneName" column of the
#' elementMetadata of the rowRanges of the "expression" SummarizedExperiment
#' object within the TENETMultiAssayExperiment object. Defaults to NA.
#' @param hypermethGplusAnalysis Set to TRUE to calculate total links by
#' gene for hypermethylated RE DNA methylation sites with G+ links. Defaults to
#' TRUE.
#' @param hypomethGplusAnalysis Set to TRUE to calculate total links by
#' gene for hypomethylated RE DNA methylation sites with G+ links. Defaults to
#' TRUE.
#' @return Returns the MultiAssayExperiment object given as the
#' TENETMultiAssayExperiment argument with an additional list of data named
#' "step6DNAMethylationSitesPerGeneTabulation" in its metadata with the output
#' of this function, which includes data frames containing significant hyper-
#' and/or hypomethylated G+ link counts per gene after all TENET steps through
#' `step5OptimizeLinks` have been run.
#' @export
#'
#' @examplesIf interactive()
#' ## This example uses the example MultiAssayExperiment provided in the
#' ## TENET.ExperimentHub package to tabulate both hyper- and hypomethylated G+
#' ## RE DNA methylation site-gene links, using genes with names provided in the
#' ## provided MultiAssayExperiment object.
#'
#' ## Load the example TENET MultiAssayExperiment object
#' ## from the TENET.ExperimentHub package
#' exampleTENETMultiAssayExperiment <-
#'     TENET.ExperimentHub::exampleTENETMultiAssayExperiment()
#'
#' ## Calculate linked RE DNA methylation sites per gene
#' returnValue <- step6DNAMethylationSitesPerGeneTabulation(
#'     TENETMultiAssayExperiment = exampleTENETMultiAssayExperiment
#' )
#'
#' ## This example also uses the example MultiAssayExperiment provided
#' ## in the TENET.ExperimentHub package, but it only runs on hypomethylated
#' ## RE DNA methylation sites.
#'
#' ## Load the example TENET MultiAssayExperiment object
#' ## from the TENET.ExperimentHub package
#' exampleTENETMultiAssayExperiment <-
#'     TENET.ExperimentHub::exampleTENETMultiAssayExperiment()
#'
#' ## Calculate linked RE DNA methylation sites per gene
#' returnValue <- step6DNAMethylationSitesPerGeneTabulation(
#'     TENETMultiAssayExperiment = exampleTENETMultiAssayExperiment,
#'     hypermethGplusAnalysis = FALSE
#' )
step6DNAMethylationSitesPerGeneTabulation <- function(
    TENETMultiAssayExperiment,
    geneAnnotationDataset = NA,
    hypermethGplusAnalysis = TRUE,
    hypomethGplusAnalysis = TRUE) {
    ## Validate the analysis types and get a vector of the ones selected
    analysisTypes <- .validateAnalysisTypes(
        hypermethGplusAnalysis, hypomethGplusAnalysis
    )

    ## Return an error message if the input MultiAssayExperiment is invalid
    .validateMultiAssayExperiment(
        TENETMultiAssayExperiment,
        needGeneName = is.na(geneAnnotationDataset)
    )

    ## Create an empty list to hold the step 6 results
    step6ResultsList <- list()

    ## Get gene IDs and names from the MAE, or gene annotation dataset if
    ## provided
    geneIDNameDF <- .getGeneIDsAndNames(
        TENETMultiAssayExperiment, geneAnnotationDataset
    )

    ## Ensure that the relevant data are available from step 5, then run
    ## tabulation of RE DNA methylation site-gene links for the selected
    ## quadrants, adding the resulting data frame to the step 6 results list
    for (hyperHypo in analysisTypes) {
        substepName <- paste0(hyperHypo, "methGplusResults")
        substepDescription <- paste0(
            hyperHypo, "methylated G+ RE DNA methylation site-gene links"
        )

        .ensureStepPresent(
            TENETMultiAssayExperiment,
            "step5OptimizeLinks",
            substepName = substepName,
            substepDescription = substepDescription,
            substepParamDescription = paste0(
                hyperHypo, "methGplusAnalysis set to TRUE"
            )
        )

        step6ResultsList[[substepName]] <- .tabulateQuadrant(
            resultsDF = TENETMultiAssayExperiment@metadata$
                step5OptimizeLinks[[substepName]],
            geneIDNameDF = geneIDNameDF,
            substepDescription = substepDescription
        )
    }

    ## Add the step 6 results list to the MultiAssayExperiment object
    TENETMultiAssayExperiment@metadata$
        step6DNAMethylationSitesPerGeneTabulation <-
        step6ResultsList

    return(TENETMultiAssayExperiment)
}
