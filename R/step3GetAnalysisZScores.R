## Internal functions used by step 3

## Internal RE DNA methylation site Z-score calculation function for use in the
## step3GetAnalysisZScores function
.calcZScoreForMethSite <- function(
    DNAMethylationSiteID,
    geneID,
    methDataOfInterest,
    expData,
    metToExpSampleConversion,
    zScoreCalculation) {
    ## Get the methylation just for the RE DNA methylation site of interest
    methSiteOfInterestMethylation <- methDataOfInterest[DNAMethylationSiteID, ]
    methSiteOfInterestMethylation <- methSiteOfInterestMethylation[
        !is.na(methSiteOfInterestMethylation)
    ]

    ## Take the input methylation dataset passed through from .zScoreCalc
    ## and split it into samples with higher and lower methylation
    siteMethylationDataLowMethNames <- names(
        methSiteOfInterestMethylation[methSiteOfInterestMethylation == TRUE]
    )
    siteMethylationDataHighMethNames <- names(
        methSiteOfInterestMethylation[methSiteOfInterestMethylation == FALSE]
    )

    ## Check if the RE DNA methylation site is viable to prevent a divide by 0
    ## issue. For non-viable methylation sites, return an NA value.
    if (length(siteMethylationDataLowMethNames) == 0) {
        return(NA)
    }

    ## Get expression for the gene of interest from the hyper- and
    ## hypomethylated samples from the given RE DNA methylation site
    geneExpressionDataLowMeth <- expData[
        geneID,
        unname(metToExpSampleConversion[siteMethylationDataLowMethNames])
    ]
    geneExpressionDataHighMeth <- expData[
        geneID,
        unname(metToExpSampleConversion[siteMethylationDataHighMethNames])
    ]

    ## Select a Z-score calculation method depending on the selected
    ## zScoreCalculation value
    if (zScoreCalculation == "twoSample") {
        ## Get the mean, variance, and number of non-NA expression samples
        ## from each group
        meanLowMeth <- mean(geneExpressionDataLowMeth, na.rm = TRUE)
        meanHighMeth <- mean(geneExpressionDataHighMeth, na.rm = TRUE)

        varianceLowMeth <- stats::var(geneExpressionDataLowMeth, na.rm = TRUE)
        varianceHighMeth <- stats::var(geneExpressionDataHighMeth, na.rm = TRUE)

        nonNALowCount <- length(
            geneExpressionDataLowMeth[!is.na(geneExpressionDataLowMeth)]
        )
        nonNAHighCount <- length(
            geneExpressionDataHighMeth[!is.na(geneExpressionDataHighMeth)]
        )

        ## Calculate the Z-score using the two-sample method
        zScore <- (meanHighMeth - meanLowMeth) / sqrt(
            (varianceHighMeth / nonNAHighCount) +
                (varianceLowMeth / nonNALowCount)
        )
    } else {
        ## Calculate the Z-score using the one-sample method
        zScore <- (
            mean(geneExpressionDataHighMeth, na.rm = TRUE) -
                mean(geneExpressionDataLowMeth, na.rm = TRUE)
        ) / sqrt(stats::var(geneExpressionDataLowMeth, na.rm = TRUE))
    }

    return(zScore)
}

## Internal function to calculate all RE DNA methylation site Z-scores for a
## given gene
.zScoreCalc <- function(
    geneID,
    methDataOfInterest,
    expData,
    metToExpSampleConversion,
    significantZScore,
    zScoreCalculation,
    sparseResults) {
    ## Use the above function to calculate Z-scores between all RE DNA
    ## methylation sites and the gene
    methSiteZScores <- vapply(
        rownames(methDataOfInterest),
        FUN = .calcZScoreForMethSite,
        FUN.VALUE = numeric(1),
        geneID = geneID,
        methDataOfInterest = methDataOfInterest,
        expData = expData,
        metToExpSampleConversion = metToExpSampleConversion,
        zScoreCalculation = zScoreCalculation
    )

    ## If the user has decided to only save sparse results, get only the
    ## significant Z-scores at both ends of the spectrum
    if (sparseResults) {
        ## Remove the infinite and NA values
        methSiteZScores <- methSiteZScores[is.finite(methSiteZScores)]

        ## Subset to only the significant Z-scores
        methSiteZScores <- methSiteZScores[
            abs(methSiteZScores) > significantZScore
        ]
    }

    ## Round the Z-scores to 4 digits to save memory and return them
    return(round(methSiteZScores, 4))
}

## Main step 3 function

#' Calculate Z-scores comparing the mean expression of each gene in the case
#' samples that are hyper- or hypomethylated for each RE DNA methylation site
#' chosen in step 2
#'
#' This function takes the identified hyper- and hypomethylated RE DNA
#' methylation sites from `step2GetDifferentiallyMethylatedSites` function and
#' calculates Z-scores comparing the mean expression of each gene in the case
#' samples that are hyper- or hypomethylated for each RE DNA methylation site,
#' according to hyper- and hypomethylation cutoffs set during the previous
#' function, to those that are not, across all hyper- or hypomethylated RE DNA
#' methylation sites, calculating Z-scores for each unique RE DNA methylation
#' site and gene combination, also known as a link.
#'
#' @param TENETMultiAssayExperiment Specify a MultiAssayExperiment object
#' containing expression and methylation SummarizedExperiment objects,
#' such as one created by the TCGADownloader function. This
#' MultiAssayExperiment object should also contain the results from
#' the `step2GetDifferentiallyMethylatedSites` function in the metadata of that
#' object.
#' @param hypermethAnalysis Set to TRUE to calculate Z-scores for
#' hypermethylated RE DNA methylation sites. Defaults to TRUE.
#' @param hypomethAnalysis Set to TRUE to calculate Z-scores for
#' hypomethylated RE DNA methylation sites. Defaults to TRUE.
#' @param includeControl Set to TRUE to include the control samples when
#' identifying hyper/hypomethylated groups and calculating Z-scores.
#' Defaults to FALSE.
#' @param TFOnly Set to TRUE to only consider genes that are accepted
#' transcription factors in "The Human Transcription Factors" by
#' Lambert et al. 2018 when calculating Z-scores. Defaults to TRUE.
#' @param zScoreCalculation Set to 'oneSample' to use a one-sample Z-score
#' calculation or 'twoSample' to use a two sample Z-score calculation. Note
#' that 'twoSample' tends to be much more lenient, and identifies many more
#' significant RE DNA methylation site-gene links. Defaults to 'oneSample'.
#' @param sparseResults Set to TRUE to save only the significant Z-scores
#' for RE DNA methylation site-gene links. **Note:** If multiple testing
#' correction will be performed in the subsequent
#' `step4SelectMostSignificantLinksPerDNAMethylationSite` function, this
#' argument should be set to FALSE. Defaults to TRUE.
#' @param pValue Specify the p-value below which Z-scores will be considered
#' significant during comparison of gene expression values between case samples
#' that are hyper- or hypomethylated and those that are not. If `sparseResults`
#' is set to TRUE, only significant Z-scores will be saved in the output
#' MultiAssayExperiment object. Defaults to 0.05.
#' @param coreCount Argument passed as the mc.cores argument to
#' mclapply. See `?parallel::mclapply` for more details. Defaults to 1.
#' @return Returns the MultiAssayExperiment object given as the
#' TENETMultiAssayExperiment argument with an additional list
#' of data named "step3GetAnalysisZScores" in its metadata with the output
#' of this function, which includes Z-scores for each TF or gene analyzed to
#' each identified hyper- and/or hypomethylated RE DNA methylation site.
#' @export
#'
#' @examplesIf interactive()
#' ## This example uses the example MultiAssayExperiment provided in the
#' ## TENET.ExperimentHub package to calculate one-sample Z-scores for links
#' ## between both the hypermethylated and hypomethylated RE DNA methylation
#' ## sites and the expression of transcription factor genes only, considering
#' ## only the expression/methylation of case samples. Only significant
#' ## Z-scores equivalent to p<0.05 will be saved to the
#' ## TENETMultiAssayExperiment object. The analysis will be performed using
#' ## one CPU core.
#'
#' ## Load the example TENET MultiAssayExperiment object
#' ## from the TENET.ExperimentHub package
#' exampleTENETMultiAssayExperiment <-
#'     TENET.ExperimentHub::exampleTENETMultiAssayExperiment()
#'
#' ## Calculate Z-scores for hyper- and hypomethylated RE DNA methylation sites
#' returnValue <- step3GetAnalysisZScores(
#'     TENETMultiAssayExperiment = exampleTENETMultiAssayExperiment
#' )
#'
#' ## This example also uses the example MultiAssayExperiment, but calculates
#' ## two-sample Z-scores for links between only hypomethylated RE DNA
#' ## methylation sites and all genes, considering the expression/methylation
#' ## of both case and control samples. For this analysis, all Z-scores are
#' ## saved to the TENETMultiAssayExperiment object (though it should be noted
#' ## this takes a large amount of memory). Only Z-scores with p-values less
#' ## than 0.1 will be considered significant. This analysis will be performed
#' ## using 8 CPU cores.
#'
#' ## Load the example TENET MultiAssayExperiment object
#' ## from the TENET.ExperimentHub package
#' exampleTENETMultiAssayExperiment <-
#'     TENET.ExperimentHub::exampleTENETMultiAssayExperiment()
#'
#' ## Calculate Z-scores for hyper- and hypomethylated RE DNA methylation sites
#' returnValue <- step3GetAnalysisZScores(
#'     TENETMultiAssayExperiment = exampleTENETMultiAssayExperiment,
#'     hypermethAnalysis = FALSE,
#'     includeControl = TRUE,
#'     TFOnly = FALSE,
#'     zScoreCalculation = "twoSample",
#'     sparseResults = FALSE,
#'     pValue = 0.1,
#'     coreCount = 8
#' )
step3GetAnalysisZScores <- function(
    TENETMultiAssayExperiment,
    hypermethAnalysis = TRUE,
    hypomethAnalysis = TRUE,
    includeControl = FALSE,
    TFOnly = TRUE,
    zScoreCalculation = "oneSample",
    sparseResults = TRUE,
    pValue = 0.05,
    coreCount = 1) {
    ## Return an error message if no analysis types have been selected
    .validateAnalysisTypes(hypermethAnalysis, hypomethAnalysis)

    ## Return an error message if an improper p-value has been selected
    if (pValue <= 0 || pValue >= 1) {
        .stopNoCall("The pValue parameter must be a value between 0 and 1.")
    }

    ## Check that the user has supplied a valid zScoreCalculation value
    if (!(zScoreCalculation %in% c("oneSample", "twoSample"))) {
        .stopNoCall(
            "zScoreCalculation is not properly defined. Please set it to ",
            'either "oneSample" or "twoSample".'
        )
    }

    ## Return an error message if the input MultiAssayExperiment is invalid
    .validateMultiAssayExperiment(TENETMultiAssayExperiment)

    ## Ensure the output data from the step2GetDifferentiallyMethylatedSites
    ## function are present in the MultiAssayExperiment object
    .ensureStepPresent(
        TENETMultiAssayExperiment, "step2GetDifferentiallyMethylatedSites"
    )

    ## Create an environment to store data from the TENET package
    TENETDataEnv <- new.env(parent = emptyenv())

    ## Get the methylation values that match with expression values
    ## using the mapping data. This assumes the methylation and expression
    ## values share a clinical data match within the mapping.
    metToExpSampleConversion <- .createMetToExpSampleConversionVector(
        TENETMultiAssayExperiment
    )

    ## Get the names of the control and case samples in the
    ## expression and methylation data
    expressionSampleNamesControl <- .getExpOrMetSamplesOfType(
        TENETMultiAssayExperiment,
        "expression",
        "Control",
        namesOnly = TRUE
    )

    expressionSampleNamesCase <- .getExpOrMetSamplesOfType(
        TENETMultiAssayExperiment,
        "expression",
        "Case",
        namesOnly = TRUE
    )

    methylationSampleNamesControl <- .getExpOrMetSamplesOfType(
        TENETMultiAssayExperiment,
        "methylation",
        "Control",
        namesOnly = TRUE
    )

    methylationSampleNamesCase <- .getExpOrMetSamplesOfType(
        TENETMultiAssayExperiment,
        "methylation",
        "Case",
        namesOnly = TRUE
    )

    ## If the user is using only TF genes, get the gene IDs of those which are
    ## present in the data and sort them. Otherwise, sort the names of all
    ## genes present.
    if (TFOnly) {
        ## Get the IDs of the human transcription factors
        utils::data(
            "humanTranscriptionFactorList",
            package = "TENET",
            envir = TENETDataEnv
        )
        humanTFs <- TENETDataEnv$humanTranscriptionFactorList

        ## Identify the TFs which are present in the "expression"
        ## SummarizedExperiment object
        genesOrdered <- sort(humanTFs[
            humanTFs %in% rownames(
                TENETMultiAssayExperiment@ExperimentList$expression@
                assays@data[[1]]
            )
        ])
    } else {
        genesOrdered <- sort(rownames(
            TENETMultiAssayExperiment@ExperimentList$expression@
            assays@data[[1]]
        ))
    }

    ## Define the samples of interest for the analysis, depending on if
    ## the user wants to include control samples or not
    if (includeControl) {
        expressionSamplesOfInterest <- c(
            expressionSampleNamesControl,
            expressionSampleNamesCase
        )
        methylationSamplesOfInterest <- c(
            methylationSampleNamesControl,
            methylationSampleNamesCase
        )
    } else {
        expressionSamplesOfInterest <- expressionSampleNamesCase
        methylationSamplesOfInterest <- methylationSampleNamesCase
    }

    ## Extract the expression samples of interest
    expData <- MultiAssayExperiment::assays(
        TENETMultiAssayExperiment@ExperimentList$expression
    )[[1]][
        ,
        expressionSamplesOfInterest
    ]

    ## Determine a significant Z-score from the input p-value
    significantZScore <- stats::qnorm(1 - pValue)

    ## Create a list in which to store the relevant results
    resultsList <- list()

    ## Generate results for the analysis types selected by the user
    if (hypermethAnalysis) {
        ## Load the cutoff values and RE DNA methylation site identities
        ## determined by the step2GetDifferentiallyMethylatedSites function
        hypermethCutoff <- TENETMultiAssayExperiment@metadata$
            step2GetDifferentiallyMethylatedSites$hypermethCutoff
        hypermethylatedSites <- sort(
            TENETMultiAssayExperiment@metadata$
                step2GetDifferentiallyMethylatedSites$
                siteIdentitiesList$hypermethylatedSites
        )

        ## Extract the hypermethylated samples of interest which meet the cutoff
        methDataHyper <- MultiAssayExperiment::assays(
            TENETMultiAssayExperiment@ExperimentList$methylation
        )[[1]][
            hypermethylatedSites,
            methylationSamplesOfInterest
        ] < hypermethCutoff

        ## Calculate the Z-scores for each gene and save them in resultsList
        hypermethZScoreResultsList <- parallel::mclapply(
            genesOrdered,
            FUN = .zScoreCalc,
            methDataOfInterest = methDataHyper,
            expData = expData,
            metToExpSampleConversion = metToExpSampleConversion,
            significantZScore = significantZScore,
            zScoreCalculation = zScoreCalculation,
            sparseResults = sparseResults,
            mc.cores = coreCount
        )
        names(hypermethZScoreResultsList) <- genesOrdered
        resultsList$hypermethResults <- hypermethZScoreResultsList
    }

    if (hypomethAnalysis) {
        ## Load the cutoff values and RE DNA methylation site identities
        ## determined by the step2GetDifferentiallyMethylatedSites function
        hypomethCutoff <- TENETMultiAssayExperiment@metadata$
            step2GetDifferentiallyMethylatedSites$hypomethCutoff
        hypomethylatedSites <- sort(
            TENETMultiAssayExperiment@metadata$
                step2GetDifferentiallyMethylatedSites$
                siteIdentitiesList$hypomethylatedSites
        )

        ## Extract the hypomethylated samples of interest which meet the cutoff
        methDataHypo <- MultiAssayExperiment::assays(
            TENETMultiAssayExperiment@ExperimentList$methylation
        )[[1]][
            hypomethylatedSites,
            methylationSamplesOfInterest
        ] < hypomethCutoff

        ## Calculate the Z-scores for each gene and save them in resultsList
        hypomethZScoreResultsList <- parallel::mclapply(
            genesOrdered,
            FUN = .zScoreCalc,
            methDataOfInterest = methDataHypo,
            expData = expData,
            metToExpSampleConversion = metToExpSampleConversion,
            significantZScore = significantZScore,
            zScoreCalculation = zScoreCalculation,
            sparseResults = sparseResults,
            mc.cores = coreCount
        )
        names(hypomethZScoreResultsList) <- genesOrdered
        resultsList$hypomethResults <- hypomethZScoreResultsList
    }

    ## Create a data frame noting the settings for the function,
    ## as they are relevant in step 4
    metadataDF <- data.frame(
        "sparseResults" = sparseResults,
        "pValue" = pValue,
        "zScoreCalculation" = zScoreCalculation,
        stringsAsFactors = FALSE
    )

    ## Add settings to the step 3 return for use in the step 4 function
    resultsList$metadata <- metadataDF

    ## Add the step 3 results list to the MultiAssayExperiment object
    TENETMultiAssayExperiment@metadata$step3GetAnalysisZScores <-
        resultsList

    return(TENETMultiAssayExperiment)
}
