## Internal functions used by step 4

## Internal function to restrict the RE DNA methylation site-gene links by
## removing links whose number of linked genes, or multiple-testing-corrected
## p-value, exceeds the given maximum
.findSignificantLinksByMethSite <- function(
    DNAMethylationSiteID,
    quadrantResultsDF,
    linksPerREDNAMethylationSiteMaximum,
    multipleTestingPValue,
    step3Metadata) {
    ## Get the subset of links to this RE DNA methylation site
    methSiteLinks <- quadrantResultsDF[
        quadrantResultsDF$DNAMethylationSiteID %in% DNAMethylationSiteID,
    ]

    ## Depending on the directionality, isolate only the negative or positive
    ## Z-scores; Z-scores of the wrong directionality will never be
    ## significant. G- links have positive Z-scores and G+ links have
    ## negative Z-scores. In TENET, we only consider G+ links.
    methSiteLinks <- methSiteLinks[methSiteLinks$zScore < 0, ]

    ## Sort the Z-scores so the most negative Z-scores are listed first
    methSiteLinks <- methSiteLinks[order(methSiteLinks$zScore), ]

    ## Depending on the value of linksPerREDNAMethylationSiteMaximum, either
    ## get the most significant n links per RE DNA methylation site, or do a
    ## multiple testing correction on those links
    if (linksPerREDNAMethylationSiteMaximum %in%
        setdiff(stats::p.adjust.methods, "none")
    ) {
        ## Since the user wants to do a multiple testing correction based on
        ## the number of links to each unique RE DNA methylation site in the
        ## quadrant, first convert the Z-scores to p-values again
        methSiteLinks$zScorePValue <- stats::pnorm(methSiteLinks$zScore)

        ## Perform the multiple testing
        methSiteLinks$zScorePValueAdj <- stats::p.adjust(
            methSiteLinks$zScorePValue,
            method = linksPerREDNAMethylationSiteMaximum
        )

        ## Create a vector of the Z-scores with significant p-values after the
        ## multiple testing, with its names set to the names of the genes
        returnVector <- methSiteLinks[
            methSiteLinks$zScorePValueAdj < multipleTestingPValue,
            "zScore"
        ]
        names(returnVector) <- methSiteLinks[
            methSiteLinks$zScorePValueAdj < multipleTestingPValue,
            "geneID"
        ]
    } else {
        ## Restrict the dataset to only the significant Z-scores from step 3 if
        ## sparseResults was not set to TRUE. Otherwise, we will end up
        ## grabbing the top n Z-scores without regard for if they were
        ## significant to begin with.
        if (!step3Metadata$sparseResults) {
            ## Get the Z-score that corresponds to the significant p-value used
            ## in step 3
            significantZScore <- stats::qnorm(
                1 - step3Metadata$pValue
            )

            ## Restrict the Z-scores to only those with absolute values above
            ## the given p-value
            methSiteLinks <- methSiteLinks[
                (abs(methSiteLinks$zScore) > significantZScore),
            ]
        }

        ## Create a vector of all Z-scores with its names set to the names of
        ## the genes. If there are are more RE DNA methylation site-gene links
        ## for the given RE DNA methylation site than the maximum number
        ## specified, return the most significant Z-scores up to that number.
        ## If there are fewer, return all Z-scores.
        returnVector <- methSiteLinks$zScore[seq_len(
            min(nrow(methSiteLinks), linksPerREDNAMethylationSiteMaximum)
        )]
        names(returnVector) <- methSiteLinks$geneID[seq_len(
            min(nrow(methSiteLinks), linksPerREDNAMethylationSiteMaximum)
        )]
    }

    return(returnVector)
}

.restrictLinksPerMethSite <- function(
    MAE,
    hyperHypo,
    linksPerREDNAMethylationSiteMaximum,
    multipleTestingPValue,
    coreCount) {
    ## Define the name of the results category in the MAE
    resultsCategory <- paste0(hyperHypo, "methResults")

    ## Ensure that the relevant data are available from step 3
    .ensureStepPresent(
        MAE,
        "step3GetAnalysisZScores",
        substepName = resultsCategory,
        substepDescription = paste0(
            hyperHypo, "methylated RE DNA methylation sites"
        ),
        substepParamDescription = paste0(
            hyperHypo, "methAnalysis set to TRUE"
        )
    )

    ## Create a data frame which will contain the results for this function,
    ## starting with data on the gene-RE DNA methylation site links from step 3.
    ## Note: c() wraps the geneID part to handle an error where that
    ## part can result in a matrix instead of a single vector, so the c()
    ## collapses it back down if it occurs.
    quadrantResults <- data.frame(
        "geneID" = c(unlist(unname(mapply(
            rep,
            names(MAE@metadata$step3GetAnalysisZScores[[resultsCategory]]),
            lengths(
                MAE@metadata$step3GetAnalysisZScores[[resultsCategory]]
            )
        )))),
        ## Extract the methylation site ID from the name, which looks like
        ## <geneID>.<methSiteID>
        "DNAMethylationSiteID" = sub(
            "^.*?\\.",
            "",
            names(unlist(
                MAE@metadata$step3GetAnalysisZScores[[resultsCategory]]
            ))
        ),
        "zScore" = unname(
            unlist(MAE@metadata$step3GetAnalysisZScores[[resultsCategory]])
        ),
        stringsAsFactors = FALSE
    )

    ## Remove NA or NaN values
    quadrantResults <- quadrantResults[!is.na(quadrantResults$zScore), ]

    ## Order the data frame alphanumerically by the methylation site IDs
    quadrantResults <- quadrantResults[
        order(quadrantResults$DNAMethylationSiteID),
    ]

    ## Use the .findSignificantLinksByMethSite function to create a list of
    ## vectors, one per unique RE DNA methylation site, with the significant
    ## Z-scores for each gene linked to the RE DNA methylation site, and the
    ## names of those genes
    restrictedResults <- parallel::mclapply(
        unique(quadrantResults$DNAMethylationSiteID),
        .findSignificantLinksByMethSite,
        quadrantResultsDF = quadrantResults,
        linksPerREDNAMethylationSiteMaximum =
            linksPerREDNAMethylationSiteMaximum,
        multipleTestingPValue = multipleTestingPValue,
        step3Metadata =
            MAE@metadata$step3GetAnalysisZScores$metadata,
        mc.cores = coreCount
    )

    ## Name the vectors after their methylation site IDs
    names(restrictedResults) <- unique(quadrantResults$DNAMethylationSiteID)

    ## Return the Z-scores and gene names
    return(restrictedResults)
}

## Main step 4 function

#' Select the most significant RE DNA methylation site-gene links to each RE DNA
#' methylation site
#'
#' This function takes the calculated Z-scores for the hyper- or hypomethylated
#' G+ RE DNA methylation site-gene links and selects the most significant links
#' to each regulatory element DNA methylation site, either up to a number
#' specified by the user, or based on a significant p-value level set by the
#' user after multiple testing correction is performed on the Z-scores output
#' by the `step3GetAnalysisZScores` function per RE DNA methylation site in the
#' RE DNA methylation site-gene pairs.
#'
#' @param TENETMultiAssayExperiment Specify a MultiAssayExperiment object
#' containing expression and methylation SummarizedExperiment objects, such as
#' one created by the TCGADownloader function. This MultiAssayExperiment object
#' should also contain the results from the
#' `step2GetDifferentiallyMethylatedSites` and `step3GetAnalysisZScores`
#' functions in its metadata.
#' @param hypermethGplusAnalysis Set to TRUE to analyze hypermethylated G+ RE
#' DNA methylation site-gene links. Requires the hypermethAnalysis parameter to
#' have been set to TRUE in step 3.
#' @param hypomethGplusAnalysis Set to TRUE to analyze hypomethylated G+ RE DNA
#' methylation site-gene links. Requires the hypomethAnalysis parameter to have
#' been set to TRUE in step 3.
#' @param linksPerREDNAMethylationSiteMaximum This parameter must either be set
#' to an integer n greater than 0, in which case only the n most significant RE
#' DNA methylation site-gene link pairs from step 3 will be selected per RE DNA
#' methylation site, or a character string describing a multiple testing
#' correction method supported by `p.adjust` (see `?stats::p.adjust`) to perform
#' multiple testing correction on the Z-scores from step 3, using the
#' multipleTestingPValue argument to set a significant p-value cutoff.
#' **Note:** If multiple testing correction is performed, `sparseResults`
#' should have been set to FALSE in the `step3GetAnalysisZScores` function.
#' Defaults to 25 (maximum links per unique RE DNA methylation site).
#' @param multipleTestingPValue Cutoff for multiple testing corrected p-values.
#' This argument is only used if the `linksPerREDNAMethylationSiteMaximum`
#' argument is set to a multiple testing correction method. Defaults to 0.05.
#' @param coreCount Argument passed as the mc.cores argument to mclapply. See
#' `?parallel::mclapply` for more details. Defaults to 1.
#' @return Returns the MultiAssayExperiment object given as the
#' TENETMultiAssayExperiment argument with an additional list of data named
#' "step4SelectMostSignificantLinksPerDNAMethylationSite" in its metadata with
#' the output of this function, which includes the most significant selected
#' gene links to the hyper- or hypomethylated RE DNA methylation sites.
#' @export
#'
#' @examplesIf interactive()
#' ## This example uses the example MultiAssayExperiment provided in the
#' ## TENET.ExperimentHub package to identify the 25 most significant links
#' ## between both hyper- and hypomethylated enhancer DNA methylation sites and
#' ## all genes, using one CPU core to perform the analysis.
#'
#' ## Load the example TENET MultiAssayExperiment object
#' ## from the TENET.ExperimentHub package
#' exampleTENETMultiAssayExperiment <-
#'     TENET.ExperimentHub::exampleTENETMultiAssayExperiment()
#'
#' ## Perform the analysis
#' returnValue <- step4SelectMostSignificantLinksPerDNAMethylationSite(
#'     TENETMultiAssayExperiment = exampleTENETMultiAssayExperiment
#' )
#'
#' ## This example also uses the example MultiAssayExperiment but identifies
#' ## the most significant links between only hypomethylated enhancer DNA
#' ## methylation sites and all genes by performing Bonferroni multiple testing
#' ## correction using a significant p-value of 0.10, using 8 CPU cores to
#' ## perform the analysis. Note: running this code with the
#' ## exampleTENETMultiAssayExperiment will produce a warning message as
#' ## sparseResults was set to TRUE when the example dataset was generated.
#' ## However, this function will still run and is valid as an example.
#'
#' ## Load the example TENET MultiAssayExperiment object
#' ## from the TENET.ExperimentHub package
#' exampleTENETMultiAssayExperiment <-
#'     TENET.ExperimentHub::exampleTENETMultiAssayExperiment()
#'
#' ## Perform the analysis
#' returnValue <- step4SelectMostSignificantLinksPerDNAMethylationSite(
#'     TENETMultiAssayExperiment = exampleTENETMultiAssayExperiment,
#'     hypermethGplusAnalysis = FALSE,
#'     linksPerREDNAMethylationSiteMaximum = "bonferroni",
#'     multipleTestingPValue = 0.1,
#'     coreCount = 8
#' )
step4SelectMostSignificantLinksPerDNAMethylationSite <- function(
    TENETMultiAssayExperiment,
    hypermethGplusAnalysis = TRUE,
    hypomethGplusAnalysis = TRUE,
    linksPerREDNAMethylationSiteMaximum = 25,
    multipleTestingPValue = 0.05,
    coreCount = 1) {
    ## Validate the analysis types and get a vector of the ones selected
    analysisTypes <- .validateAnalysisTypes(
        hypermethGplusAnalysis, hypomethGplusAnalysis
    )

    ## Validate the linksPerREDNAMethylationSiteMaximum value or the multiple
    ## testing parameters
    if (linksPerREDNAMethylationSiteMaximum %in%
        setdiff(stats::p.adjust.methods, "none")
    ) {
        ## multipleTestingPValue must be greater than 0 and less than 1, as
        ## it is a p-value
        if (multipleTestingPValue <= 0 || multipleTestingPValue >= 1) {
            .stopNoCall(
                "The multipleTestingPValue argument must be a value ",
                "greater than 0 and less than 1, as it is the significant ",
                "p-value for Z-scores from step 3 after multiple testing ",
                "correction for the number of genes linked to each ",
                "RE DNA methylation site among the RE DNA methylation ",
                "site-gene link pairs."
            )
        }

        ## If the sparseResults argument was set to TRUE in step 3, that will
        ## affect how multiple testing correction is performed
        if (TENETMultiAssayExperiment@metadata$
            step3GetAnalysisZScores$metadata$sparseResults
        ) {
            warning(
                "The sparseResults argument was set to TRUE in the ",
                "step3GetAnalysisZScores function, so only significant ",
                "Z-scores with equivalent p-values below the pValue argument ",
                "for that function were saved. This will affect the results ",
                "of the multiple testing method selected. Consider re-running ",
                "the step3GetAnalysisZScores function with sparseResults ",
                "set to FALSE if you want to perform multiple testing ",
                "correction."
            )
        }
    } else {
        ## If linksPerREDNAMethylationSiteMaximum is not set to a method
        ## supported by stats::p.adjust, it must be greater than 0 and an
        ## integer
        if (!is.numeric(linksPerREDNAMethylationSiteMaximum) ||
            linksPerREDNAMethylationSiteMaximum <= 0 ||
            floor(linksPerREDNAMethylationSiteMaximum) !=
                linksPerREDNAMethylationSiteMaximum
        ) {
            .stopNoCall(
                "The linksPerREDNAMethylationSiteMaximum argument must either ",
                "be set to an integer n greater than 0, in which case only ",
                "the n most significant RE DNA methylation site-gene link ",
                "pairs from step 3 will be selected per RE DNA methylation ",
                "site, or a character string describing a multiple testing ",
                "correction method supported by p.adjust (see ",
                "?stats::p.adjust) to perform multiple testing correction on ",
                "the Z-scores from step 3, using the multipleTestingPValue ",
                "argument to set a significant p-value cutoff."
            )
        }
    }

    ## Return an error message if the input MultiAssayExperiment is invalid
    .validateMultiAssayExperiment(TENETMultiAssayExperiment)

    ## Ensure the output data from the step2GetDifferentiallyMethylatedSites
    ## function are present in the MultiAssayExperiment object
    .ensureStepPresent(
        TENETMultiAssayExperiment, "step2GetDifferentiallyMethylatedSites"
    )

    ## Ensure the output data from the step3GetAnalysisZScores function
    ## are present in the MultiAssayExperiment object
    .ensureStepPresent(
        TENETMultiAssayExperiment, "step3GetAnalysisZScores"
    )

    ## Create an empty list to store the results of this step 4 function
    TENETMultiAssayExperiment@metadata$
        step4SelectMostSignificantLinksPerDNAMethylationSite <- list()

    ## Perform the selection of links per RE DNA methylation site on the
    ## selected quadrants
    for (hyperHypo in analysisTypes) {
        TENETMultiAssayExperiment@metadata$
            step4SelectMostSignificantLinksPerDNAMethylationSite[[
            paste0(hyperHypo, "methGplusResults")
        ]] <- .restrictLinksPerMethSite(
            MAE = TENETMultiAssayExperiment,
            hyperHypo = hyperHypo,
            linksPerREDNAMethylationSiteMaximum =
                linksPerREDNAMethylationSiteMaximum,
            multipleTestingPValue = multipleTestingPValue,
            coreCount = coreCount
        )
    }

    return(TENETMultiAssayExperiment)
}
