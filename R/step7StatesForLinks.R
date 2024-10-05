## Internal functions used by step7StatesForLinks

## Internal function to calculate "personal links". This determines if a given
## case sample has a given RE DNA methylation site-gene link based on if it has
## expression significantly greater than the control sample mean and methylation
## greater than the hypermethylation cutoff, or expression significantly less
## than the control sample mean and methylation less than the hypomethylation
## cutoff, as determined by a Bonferroni-corrected 1-sided t-test with a p-value
## threshold of 0.05.
.linkEvaluator <- function(
    DNAMethylationSiteID,
    geneID,
    expressionDataControl,
    expressionDataCase,
    quadrantMethylationDataCase,
    hyperHypo,
    cutoff) {
    ## Extract methylation values from each case sample
    caseMethylationValuesForSite <- unlist(
        quadrantMethylationDataCase[DNAMethylationSiteID, ]
    )

    ## Extract gene expression values for each case sample
    caseExpressionValuesForGene <- unlist(
        expressionDataCase[geneID, ]
    )

    ## Get expression values from the control samples
    controlExpressionValuesForGene <- unlist(
        expressionDataControl[geneID, ]
    )

    ## Determine for which samples the mean expression is significantly lower
    ## or higher than the control mean expression
    pValuesComparedToControl <- numeric()

    for (i in caseExpressionValuesForGene) {
        pValuesComparedToControl <- c(
            pValuesComparedToControl,
            stats::t.test(
                controlExpressionValuesForGene,
                mu = i,
                alternative = ifelse(hyperHypo == "hyper", "greater", "less")
            )$p.value
        )
    }

    ## Bonferroni correct the p-values
    pValuesComparedToControlBonferroni <- stats::p.adjust(
        pValuesComparedToControl,
        method = "bonferroni",
        n = length(pValuesComparedToControl)
    )

    ## Create a vector that is TRUE if the given sample has expression for the
    ## gene significantly different from the control samples
    significantFromControl <- (pValuesComparedToControlBonferroni < 0.05)

    ## For each methylation sample, check if it is hyper- or hypomethylated
    ## for the given RE DNA methylation site
    if (hyperHypo == "hypo") {
        cutoffValues <- (caseMethylationValuesForSite < cutoff)
    } else if (hyperHypo == "hyper") {
        cutoffValues <- (caseMethylationValuesForSite > cutoff)
    }

    ## Create a final vector with 1 for where both conditions are met, and 0
    ## when they aren't both met
    finalOutput <- as.numeric(significantFromControl & cutoffValues)

    ## Assign sample names to the output values
    names(finalOutput) <- colnames(expressionDataCase)

    ## Return the final output
    return(finalOutput)
}

## Main step7StatesForLinks function

#' Identify which of the case samples harbor each of the identified regulatory
#' element DNA methylation site-gene links
#'
#' This function generates data frames with information for each of the case
#' samples and each RE DNA methylation site-gene link from both of the
#' hyper- or hypomethylated G+ analysis quadrants, as selected by the user, on
#' if a given sample is said to "harbor" each link, depending on if the
#' methylation of the given sample for the RE DNA methylation site in the link
#' is above or below the hyper- or hypomethylation cutoff defined in the
#' `step2GetDifferentiallyMethylatedSites` function and the expression of the
#' gene in the link is significantly greater than, or less than, the mean
#' expression in the control samples, as determined by a Bonferroni-corrected
#' 1-sided t-test with a p-value threshold of 0.05.
#'
#' @param TENETMultiAssayExperiment Specify a MultiAssayExperiment object
#' containing expression and methylation SummarizedExperiment objects, such as
#' one created by the TCGADownloader function. This MultiAssayExperiment object
#' should also contain the results from the `step5OptimizeLinks` functions in
#' its metadata.
#' @param hypermethGplusAnalysis Set to TRUE to analyze which case samples
#' harbor each hypermethylated G+ RE DNA methylation site-gene link. Defaults to
#' TRUE.
#' @param hypomethGplusAnalysis Set to TRUE to analyze which case samples
#' harbor each hypomethylated G+ RE DNA methylation site-gene link. Defaults to
#' TRUE.
#' @param coreCount Argument passed as the mc.cores argument to mcmapply. See
#' `?parallel::mcmapply` for more details. Defaults to 1.
#' @return Returns the MultiAssayExperiment object given as the
#' TENETMultiAssayExperiment argument with an additional list of information
#' named 'step7StatesForLinks' in its metadata with the output of this
#' function. This list is subdivided into hypermethGplus or hypomethGplus
#' results as selected by the user, which contain data frames for each of the
#' specified analysis types with case samples in the columns and each RE DNA
#' methylation site-gene link for that analysis type in the rows, with a 1
#' indicating the sample is positive for that link and a 0 indicating it is not.
#' NA values are shown for samples that lack methylation data for the site or
#' expression data for the gene.
#' @export
#'
#' @examplesIf interactive()
#' ## This example uses the example MultiAssayExperiment provided in the
#' ## TENET.ExperimentHub package to create data frames for link states between
#' ## genes and both hyper- and hypomethylated G+ RE DNA methylation sites. The
#' ## analysis will be performed using one CPU core.
#'
#' ## Load the example TENET MultiAssayExperiment object
#' ## from the TENET.ExperimentHub package
#' exampleTENETMultiAssayExperiment <-
#'     TENET.ExperimentHub::exampleTENETMultiAssayExperiment()
#'
#' ## Use the example dataset to do the states for links analysis
#' returnValue <- step7StatesForLinks(
#'     TENETMultiAssayExperiment = exampleTENETMultiAssayExperiment
#' )
#'
#' ## This example uses the example MultiAssayExperiment provided in the
#' ## TENET.ExperimentHub package to create data frames for link states between
#' ## genes and only hypomethylated G+ RE DNA methylation sites. The analysis
#' ## will be performed using 8 CPU cores.
#'
#' ## Load the example TENET MultiAssayExperiment object
#' ## from the TENET.ExperimentHub package
#' exampleTENETMultiAssayExperiment <-
#'     TENET.ExperimentHub::exampleTENETMultiAssayExperiment()
#'
#' ## Use the example dataset to do the states for links analysis
#' returnValue <- step7StatesForLinks(
#'     TENETMultiAssayExperiment = exampleTENETMultiAssayExperiment,
#'     hypermethGplusAnalysis = FALSE,
#'     coreCount = 8
#' )
step7StatesForLinks <- function(
    TENETMultiAssayExperiment,
    hypermethGplusAnalysis = TRUE,
    hypomethGplusAnalysis = TRUE,
    coreCount = 1) {
    ## Validate the analysis types and get a vector of the ones selected
    analysisTypes <- .validateAnalysisTypes(
        hypermethGplusAnalysis, hypomethGplusAnalysis
    )

    ## Return an error message if the input MultiAssayExperiment is invalid
    .validateMultiAssayExperiment(TENETMultiAssayExperiment)

    ## Ensure the output data from the step2GetDifferentiallyMethylatedSites
    ## function are present in the MultiAssayExperiment object
    .ensureStepPresent(
        TENETMultiAssayExperiment, "step2GetDifferentiallyMethylatedSites"
    )

    ## Get the methylation names of the case samples (control methylation data
    ## are not used here)
    methylationSampleNamesCase <- .getExpOrMetSamplesOfType(
        TENETMultiAssayExperiment,
        "methylation",
        "Case",
        namesOnly = TRUE
    )

    ## Get the control and case expression datasets
    expressionDataControl <- .getExpOrMetSamplesOfType(
        TENETMultiAssayExperiment,
        "expression",
        "Control"
    )

    expressionDataCase <- .getExpOrMetSamplesOfType(
        TENETMultiAssayExperiment,
        "expression",
        "Case"
    )

    ## Create an empty list to hold the results from this step 7 function
    resultsList <- list()

    ## Generate results for the selected analysis types
    for (hyperHypo in analysisTypes) {
        quadrantResultsName <- paste0(hyperHypo, "methGplusResults")

        ## Ensure the quadrant's results are present in step 5
        .ensureStepPresent(
            TENETMultiAssayExperiment,
            stepName = "step5OptimizeLinks",
            substepName = quadrantResultsName
        )

        ## Load the quadrant's significant links from step 5
        quadrantSigLinkZScores <- TENETMultiAssayExperiment@metadata$
            step5OptimizeLinks[[quadrantResultsName]]

        ## Get the quadrant's cutoff value from step 2
        quadrantCutoff <- unname(
            TENETMultiAssayExperiment@metadata$
                step2GetDifferentiallyMethylatedSites
            [[paste0(hyperHypo, "methCutoff")]]
        )

        ## Get the quadrant's RE DNA methylation site identities from step 2
        quadrantMethSites <- sort(
            TENETMultiAssayExperiment@metadata$
                step2GetDifferentiallyMethylatedSites$
                siteIdentitiesList[[paste0(hyperHypo, "methylatedSites")]]
        )

        ## Get the quadrant's case methylation dataset
        quadrantMethylationDataCase <- MultiAssayExperiment::assays(
            TENETMultiAssayExperiment@ExperimentList$methylation
        )[[1]][
            quadrantMethSites,
            methylationSampleNamesCase
        ]

        ## Calculate the sample states for the quadrant's links
        quadrantLinksDataset <- parallel::mcmapply(
            .linkEvaluator,
            DNAMethylationSiteID = quadrantSigLinkZScores$DNAMethylationSiteID,
            geneID = quadrantSigLinkZScores$geneID,
            MoreArgs = list(
                expressionDataControl = expressionDataControl,
                expressionDataCase = expressionDataCase,
                quadrantMethylationDataCase = quadrantMethylationDataCase,
                hyperHypo = hyperHypo,
                cutoff = quadrantCutoff
            ),
            mc.cores = coreCount
        )

        ## Change the column names to the names of the RE DNA methylation
        ## site-gene links
        colnames(quadrantLinksDataset) <- paste(
            quadrantSigLinkZScores$DNAMethylationSiteID,
            quadrantSigLinkZScores$geneID,
            sep = "_"
        )

        ## Transpose the dataset to make the case samples the columns
        quadrantLinksDatasetTransposed <- t(quadrantLinksDataset)

        ## Convert it to a data frame
        quadrantLinksDatasetTransposedDF <- as.data.frame(
            quadrantLinksDatasetTransposed,
            stringsAsFactors = FALSE
        )

        ## Add it to the results list
        resultsList[[quadrantResultsName]] <-
            quadrantLinksDatasetTransposedDF
    }

    ## Add the results list to the MultiAssayExperiment
    TENETMultiAssayExperiment@metadata$step7StatesForLinks <- resultsList

    return(TENETMultiAssayExperiment)
}
