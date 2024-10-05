## Internal functions used by step 5

## Internal function to calculate a Wilcoxon p-value for expression of a
## given gene linked to an RE DNA methylation site between the control and case
## samples
.getWilcoxonPValueForExpressionInControlVsCase <- function(
    geneID,
    methSiteID,
    cutoffInfo,
    expressionDataControl,
    expressionDataCase,
    metToExpSampleConversion) {
    ## Subset the cutoff info to just this RE DNA methylation site
    cutoffInfoForMethSite <- cutoffInfo[methSiteID, ]

    ## Remove any NA values
    cutoffInfoForMethSite <- cutoffInfoForMethSite[
        !is.na(cutoffInfoForMethSite)
    ]

    ## Get the names of the samples for the RE DNA methylation site that meet
    ## the cutoff
    namesOfSamplesMeetingCutoff <- names(
        cutoffInfoForMethSite[cutoffInfoForMethSite == TRUE]
    )

    ## Get the expression for the gene in the control samples
    controlExpression <- expressionDataControl[geneID, ]

    ## Get expression values for only the expressed control samples
    expressedControlSamples <- controlExpression[controlExpression != 0]

    ## Get the expression for the gene in the case samples
    caseExpression <- expressionDataCase[geneID, ]

    ## Get names of the case samples that have expression for the gene
    expressedCaseNames <- names(caseExpression[caseExpression != 0])

    ## Use the converter to get the names of the expression samples that
    ## correspond with the identified methylation samples
    expressedCaseNamesMeetingCutoff <- unname(
        metToExpSampleConversion[namesOfSamplesMeetingCutoff]
    )

    ## Get samples that have expression of the gene *and* that match the
    ## methylation criterion described by cutoffInfo
    matchedSamples <- intersect(
        expressedCaseNames, expressedCaseNamesMeetingCutoff
    )

    ## Get expression values for only the expressed case samples
    expressedCaseSamples <- caseExpression[matchedSamples]

    ## If there are samples expressing it in both datasets, return the
    ## Wilcoxon p-value comparing their differential expression. Otherwise,
    ## return NA
    if (
        length(expressedControlSamples) > 0 &&
            length(expressedCaseSamples) > 0
    ) {
        wilcoxonPValue <- stats::wilcox.test(
            expressedControlSamples,
            expressedCaseSamples,
            exact = FALSE
        )$p.value
    } else {
        wilcoxonPValue <- NA
    }

    ## Return the final Wilcoxon p-value
    return(wilcoxonPValue)
}

## Internal function to get only the hypo/hypermethylated samples for an RE DNA
## methylation site
.getSamplesMeetingCutoff <- function(methSiteID, cutoffInfo) {
    ## Subset the cutoff info to just this RE DNA methylation site
    cutoffInfoForMethSite <- cutoffInfo[methSiteID, ]

    ## Remove any NA values
    cutoffInfoForMethSite <- cutoffInfoForMethSite[
        !is.na(cutoffInfoForMethSite)
    ]

    ## Return the samples for the RE DNA methylation site that meet the cutoff
    return(cutoffInfoForMethSite[cutoffInfoForMethSite == TRUE])
}

## Internal function to calculate the number of hypo/hypermethylated samples
## for an RE DNA methylation site
.getNumOfSamplesMeetingCutoff <- function(methSiteID, cutoffInfo) {
    return(length(.getSamplesMeetingCutoff(methSiteID, cutoffInfo)))
}

## Internal function to get the mean expression of the gene of interest
## in the hyper/hypomethylated samples with expression of the gene
.getMeanExpressionInSamplesMeetingCutoffWithExpression <- function(
    geneID,
    methSiteID,
    cutoffInfo,
    expressionDataCase,
    metToExpSampleConversion) {
    ## Get the samples for the RE DNA methylation site that meet the cutoff
    samplesMeetingCutoff <- .getSamplesMeetingCutoff(methSiteID, cutoffInfo)

    ## Get the names of the samples for the RE DNA methylation site that meet
    ## the cutoff
    namesOfSamplesMeetingCutoff <- names(samplesMeetingCutoff)

    ## Get the expression for the gene in the case samples
    caseExpression <- expressionDataCase[geneID, ]

    ## Get names of the case samples that have expression for the gene
    expressedCaseNames <- names(caseExpression[caseExpression != 0])

    ## Use the conversion table to get the names of the expression
    ## samples that correspond with the identified methylation samples
    expressedCaseNamesMeetingCutoff <- unname(
        metToExpSampleConversion[namesOfSamplesMeetingCutoff]
    )

    ## Get samples that have expression of the gene *and* that match the
    ## methylation criterion described by cutoffInfo
    matchedSamples <- intersect(
        expressedCaseNames, expressedCaseNamesMeetingCutoff
    )

    ## Get expression values for only the expressed case samples
    expressedCaseSamples <- caseExpression[matchedSamples]

    ## Return the mean expression in these samples
    return(mean(expressedCaseSamples, na.rm = TRUE))
}

## Internal function calculating the number of hypo/hypermethylated samples
## greater/less than the case expression mean
.getNumSamplesMeetingCutoffWithExpressionGTOrLTCaseMean <- function(
    geneID,
    methSiteID,
    cutoffInfo,
    expressionDataCase,
    metToExpSampleConversion,
    hyperHypo) {
    ## Get the samples for the RE DNA methylation site that meet the cutoff
    samplesMeetingCutoff <- .getSamplesMeetingCutoff(methSiteID, cutoffInfo)

    ## Get the names of the samples for the RE DNA methylation site that meet
    ## the cutoff
    namesOfSamplesMeetingCutoff <- names(samplesMeetingCutoff)

    ## Get the expression for the gene in the case samples
    caseExpression <- expressionDataCase[geneID, ]

    ## Calculate the mean of these samples
    meanExpressionCase <- mean(caseExpression, na.rm = TRUE)

    ## Get names of the case samples that have expression for the gene
    expressedCaseNames <- names(caseExpression[caseExpression != 0])

    ## Use the conversion table to get the names of the expression samples
    ## that correspond with the identified methylation samples
    expressedCaseNamesMeetingCutoff <- unname(
        metToExpSampleConversion[namesOfSamplesMeetingCutoff]
    )

    ## Get samples that have expression of the gene *and* that match the
    ## methylation criterion described by cutoffInfo
    matchedSamples <- intersect(
        expressedCaseNames, expressedCaseNamesMeetingCutoff
    )

    ## Get expression values for only the expressed case samples
    expressedCaseSamples <- caseExpression[matchedSamples]

    if (hyperHypo == "hyper") {
        ## Return the number of hypermethylated samples less than the mean
        return(
            sum(expressedCaseSamples < meanExpressionCase, na.rm = TRUE)
        )
    } else {
        ## Return the number of hypomethylated samples greater than the mean
        return(
            sum(expressedCaseSamples > meanExpressionCase, na.rm = TRUE)
        )
    }
}

## Internal function to get the min/max methylation value for
## hyper/hypomethylated RE DNA methylation sites
.getSiteMinOrMaxMethValue <- function(
    methSiteID,
    cutoffInfo,
    methylationDataCase,
    hyperHypo) {
    ## Get the samples for the RE DNA methylation site that meet the cutoff
    samplesMeetingCutoff <- .getSamplesMeetingCutoff(methSiteID, cutoffInfo)

    ## Get the names of the samples for the RE DNA methylation site that meet
    ## the cutoff
    namesOfSamplesMeetingCutoff <- names(samplesMeetingCutoff)

    ## Return the min/max methylation value for these hyper/hypomethylated
    ## RE DNA methylation sites
    minMaxFunc <- ifelse(hyperHypo == "hyper", max, min)
    return(minMaxFunc(
        methylationDataCase[
            methSiteID,
            match(
                namesOfSamplesMeetingCutoff, colnames(methylationDataCase)
            )
        ],
        na.rm = TRUE
    ))
}

## Internal function to perform link optimization on a given data quadrant
.optimizeQuadrantLinks <- function(
    MAE,
    hyperHypo,
    expressionDataControl,
    expressionDataCase,
    methylationSampleNamesControl,
    methylationSampleNamesCase,
    minCaseCount,
    stringency = NA,
    coreCount,
    expressionPvalCutoff,
    metToExpSampleConversion) {
    ## Define variable names based on the selected methylation type
    caseLengthName <- paste0(hyperHypo, "methCaseLength")
    lowHighCaseLengthName <- paste0(
        hyperHypo, "meth", ifelse(hyperHypo == "hyper", "Low", "High"),
        "erCaseLength"
    )
    meanLowHighName <- paste0(
        "meanExpressionControl", ifelse(hyperHypo == "hyper", "High", "Low"),
        "ExpressionCase"
    )
    minMaxName <- paste0(
        ifelse(hyperHypo == "hyper", "max", "min"), "MethCaseValue"
    )
    meanHyperHypoName <- paste0(
        "mean", tools::toTitleCase(hyperHypo), "methCase"
    )
    resultsCategory <- paste0(hyperHypo, "methResults")
    resultsName <- paste0(hyperHypo, "methGplusResults")

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

    ## Ensure that the relevant data are available from step 4
    .ensureStepPresent(
        MAE,
        "step4SelectMostSignificantLinksPerDNAMethylationSite",
        substepName = resultsName,
        substepDescription = paste0(
            hyperHypo, "methylated G+ RE DNA methylation site-gene links"
        ),
        substepParamDescription = paste0(
            hyperHypo, "methGplusAnalysis set to TRUE"
        )
    )

    ## Get the quadrant's cutoff value from step 2
    quadrantCutoff <- unname(
        MAE@metadata$step2GetDifferentiallyMethylatedSites[[
            paste0(hyperHypo, "methCutoff")
        ]]
    )

    ## If stringency values aren't set, and the analysis types relevant to
    ## them are selected to be performed, set the stringency cutoffs to be
    ## the hyper/hypomethylation cutoffs determined or set in step 2
    if (is.na(stringency)) {
        stringency <- quadrantCutoff
    }

    ## Get the quadrant's RE DNA methylation site identities from step 2
    quadrantMethSites <- sort(
        MAE@metadata$step2GetDifferentiallyMethylatedSites$siteIdentitiesList[[
            paste0(hyperHypo, "methylatedSites")
        ]]
    )

    ## Create a data frame which will contain the results for this function,
    ## starting with data on the gene-RE DNA methylation site links from steps
    ## 3 and 4.
    step4QuadrantResults <- MAE@metadata$
        step4SelectMostSignificantLinksPerDNAMethylationSite[[resultsName]]

    ## Note: c() wraps the DNAMethylationSiteID part to handle a rare error
    ## where that part can result in a matrix instead of a single vector, so
    ## the c() collapses it back down if it occurs.
    quadrantResults <- data.frame(
        "geneID" = sub(
            "^.*?\\.",
            "",
            names(unlist(step4QuadrantResults))
        ),
        "DNAMethylationSiteID" = c(unlist(unname(mapply(
            rep,
            names(step4QuadrantResults),
            lengths(step4QuadrantResults)
        )))),
        "zScore" = unname(unlist(step4QuadrantResults)),
        "pValue" = stats::pnorm(unname(unlist(step4QuadrantResults))),
        stringsAsFactors = FALSE
    )

    ## Sort the data frame by gene ID then by methylation site ID
    quadrantResults <- quadrantResults[
        with(quadrantResults, order(geneID, DNAMethylationSiteID)),
    ]

    ## Get the case methylation dataset
    quadrantMethylationDataCase <- MultiAssayExperiment::assays(
        MAE@ExperimentList$methylation
    )[[1]][
        quadrantMethSites,
        methylationSampleNamesCase
    ]

    ## Create empty columns with NA values
    quadrantResults$meanExpressionControl <- NA
    quadrantResults$meanExpressionCase <- NA
    quadrantResults$wilcoxonPValExpressionControlVsCase <- NA
    quadrantResults[[caseLengthName]] <- NA
    quadrantResults[[meanHyperHypoName]] <- NA
    quadrantResults[[lowHighCaseLengthName]] <- NA
    quadrantResults[[minMaxName]] <- NA
    quadrantResults[[meanLowHighName]] <- NA

    if (hyperHypo == "hyper") {
        ## For each RE DNA methylation site, note which of the samples are above
        ## the hypermethylation cutoff
        cutoffInfo <- (quadrantMethylationDataCase > quadrantCutoff)
    } else {
        ## For each RE DNA methylation site, note which of the samples are below
        ## the hypomethylation cutoff
        cutoffInfo <- (quadrantMethylationDataCase < quadrantCutoff)
    }

    ## Calculate mean expression in case vs control samples
    quadrantResults$meanExpressionControl <- apply(
        expressionDataControl[quadrantResults$geneID, ], 1, mean,
        na.rm = TRUE
    )

    quadrantResults$meanExpressionCase <- apply(
        expressionDataCase[quadrantResults$geneID, ], 1, mean,
        na.rm = TRUE
    )

    quadrantResults$wilcoxonPValExpressionControlVsCase <- unname(
        parallel::mcmapply(
            .getWilcoxonPValueForExpressionInControlVsCase,
            geneID = quadrantResults$geneID,
            methSiteID = quadrantResults$DNAMethylationSiteID,
            MoreArgs = list(
                "cutoffInfo" = cutoffInfo,
                "expressionDataControl" = expressionDataControl,
                "expressionDataCase" = expressionDataCase,
                "metToExpSampleConversion" = metToExpSampleConversion
            ),
            mc.cores = coreCount
        )
    )

    quadrantResults[[caseLengthName]] <- unname(
        parallel::mcmapply(
            .getNumOfSamplesMeetingCutoff,
            methSiteID = quadrantResults$DNAMethylationSiteID,
            MoreArgs = list("cutoffInfo" = cutoffInfo),
            mc.cores = coreCount
        )
    )

    quadrantResults[[meanHyperHypoName]] <- parallel::mcmapply(
        .getMeanExpressionInSamplesMeetingCutoffWithExpression,
        geneID = quadrantResults$geneID,
        methSiteID = quadrantResults$DNAMethylationSiteID,
        MoreArgs = list(
            "cutoffInfo" = cutoffInfo,
            "expressionDataCase" = expressionDataCase,
            "metToExpSampleConversion" = metToExpSampleConversion
        ),
        mc.cores = coreCount
    )

    quadrantResults[[lowHighCaseLengthName]] <- parallel::mcmapply(
        .getNumSamplesMeetingCutoffWithExpressionGTOrLTCaseMean,
        geneID = quadrantResults$geneID,
        methSiteID = quadrantResults$DNAMethylationSiteID,
        MoreArgs = list(
            "cutoffInfo" = cutoffInfo,
            "expressionDataCase" = expressionDataCase,
            "metToExpSampleConversion" = metToExpSampleConversion,
            "hyperHypo" = hyperHypo
        ),
        mc.cores = coreCount
    )

    quadrantResults[[minMaxName]] <- parallel::mcmapply(
        .getSiteMinOrMaxMethValue,
        methSiteID = quadrantResults$DNAMethylationSiteID,
        MoreArgs = list(
            "cutoffInfo" = cutoffInfo,
            "methylationDataCase" = quadrantMethylationDataCase,
            "hyperHypo" = hyperHypo
        ),
        mc.cores = coreCount
    )

    if (hyperHypo == "hyper") {
        quadrantResults[[meanLowHighName]] <- (
            quadrantResults$meanExpressionControl >
                quadrantResults$meanExpressionCase
        )
    } else {
        quadrantResults[[meanLowHighName]] <- (
            quadrantResults$meanExpressionControl <
                quadrantResults$meanExpressionCase
        )
    }

    ## Perform BH multiple testing correction on the Wilcoxon p-values
    quadrantResults$wilcoxonPValExpressionControlVsCaseAdjusted <-
        stats::p.adjust(
            quadrantResults$wilcoxonPValExpressionControlVsCase, "BH"
        )

    ## Select the final links to carry over into post-hoc analysis
    quadrantResults <- quadrantResults[
        which(
            quadrantResults$meanExpressionControl != 0 &
                quadrantResults$meanExpressionCase != 0 &
                quadrantResults[[meanLowHighName]]
        ),
    ]

    quadrantResults <- quadrantResults[
        which(
            quadrantResults$
                wilcoxonPValExpressionControlVsCaseAdjusted <
                expressionPvalCutoff &
                quadrantResults[[caseLengthName]] >= minCaseCount &
                quadrantResults[[lowHighCaseLengthName]] >=
                    minCaseCount &
                .ifelseNoIterate(
                    hyperHypo == "hyper",
                    quadrantResults[[minMaxName]] > stringency,
                    quadrantResults[[minMaxName]] < stringency
                )
        ),
    ]

    ## Return the final data frame
    return(quadrantResults)
}

## Main step 5 function

#' Find final RE DNA methylation site-gene links using various optimization
#' metrics
#'
#' This function takes the most significant hyper- or hypomethylated G+ RE DNA
#' methylation site-gene links selected in step 4, and selects optimized links
#' based on the relative expression of the given gene in hyper- or
#' hypomethylated case samples compared to control samples, using a Wilcoxon
#' test to check that the hyper- or hypomethylated samples for that given RE
#' DNA methylation site-gene link also show appropriately higher/lower
#' expression of the linked gene in a number of case samples greater than or
#' equal to the `minCaseCount` number specified in the
#' `step2GetDifferentiallyMethylatedSites` function that have maximum/minimum
#' methylation above/below the `hyperStringency`/`hypoStringency` cutoff values
#' selected.
#'
#' @param TENETMultiAssayExperiment Specify a MultiAssayExperiment object
#' containing expression and methylation SummarizedExperiment objects, such
#' as one created by the TCGADownloader function. This MultiAssayExperiment
#' object should also contain the results from the
#' `step2GetDifferentiallyMethylatedSites`, `step3GetAnalysisZScores`, and
#' `step4SelectMostSignificantLinksPerDNAMethylationSite` functions in its
#' metadata.
#' @param hypermethGplusAnalysis Set to TRUE to optimize hypermethylated
#' G+ links. Requires hypermethAnalysis from step 4 to have been set to TRUE.
#' @param hypomethGplusAnalysis Set to TRUE to optimize hypomethylated
#' G+ links. Requires hypomethAnalysis from step 4 to have been set to TRUE.
#' @param expressionPvalCutoff Cutoff for Benjamini-Hochberg corrected
#' Wilcoxon p-values used during comparison of gene expression values between
#' hyper/hypomethylated case and control samples. Defaults to 0.05.
#' @param hyperStringency Specify a number from 0 to 1 to be the beta-value
#' cutoff to optimize for hypermethylated links with methylation values above
#' the cutoff if a more/less selective cutoff is desired. Defaults to the
#' hypermethCutoff value specified in step2GetDifferentiallyMethylatedSites.
#' @param hypoStringency Specify a number from 0 to 1 to be the beta-value
#' cutoff to optimize for hypomethylated links with methylation values below
#' the cutoff if a more/less selective cutoff is desired. Defaults to the
#' hypomethCutoff value specified in step2GetDifferentiallyMethylatedSites.
#' @param coreCount Argument passed as the mc.cores argument to mcmapply.
#' See `?parallel::mcmapply` for more details. Defaults to 1.
#' @return Returns the MultiAssayExperiment object given as the
#' TENETMultiAssayExperiment argument with an additional list of data named
#' "step5OptimizeLinks" in its metadata with the output of this function, which
#' includes the gene and RE DNA methylation site IDs of each optimized gene-RE
#' DNA methylation site link, the Z-scores and corresponding p-values
#' calculated in steps 3-4 for those links, and the various optimization
#' metrics performed by this step.
#' @export
#'
#' @examplesIf interactive()
#' ## This example uses the example MultiAssayExperiment provided in the
#' ## TENET.ExperimentHub package to optimize links between the hypermethylated
#' ## and hypomethylated RE DNA methylation sites and the expression of
#' ## transcription factor genes. For this analysis, the default p-value cutoff
#' ## of 0.05 is used, with `hyperStringency` and `hypoStringency` values set
#' ## to the originally calculated hyper- and hypomethylation cutoffs. The
#' ## analysis uses one CPU core.
#'
#' ## Load the example TENET MultiAssayExperiment object
#' ## from the TENET.ExperimentHub package
#' exampleTENETMultiAssayExperiment <-
#'     TENET.ExperimentHub::exampleTENETMultiAssayExperiment()
#'
#' ## Perform the link optimization
#' returnValue <- step5OptimizeLinks(
#'     TENETMultiAssayExperiment = exampleTENETMultiAssayExperiment
#' )
#'
#' ## This example also uses the example MultiAssayExperiment provided in the
#' ## TENET.ExperimentHub package, but it only runs on hypomethylated RE DNA
#' ## methylation sites and uses a p-value cutoff of 0.01. The `hypoStringency`
#' ## value is set to 0.3 specifically rather than using the originally set
#' ## hypomethylation cutoff. Eight CPU cores are used to perform the analysis.
#'
#' ## Load the example TENET MultiAssayExperiment object
#' ## from the TENET.ExperimentHub package
#' exampleTENETMultiAssayExperiment <-
#'     TENET.ExperimentHub::exampleTENETMultiAssayExperiment()
#'
#' ## Perform the link optimization
#' returnValue <- step5OptimizeLinks(
#'     TENETMultiAssayExperiment = exampleTENETMultiAssayExperiment,
#'     hypermethGplusAnalysis = FALSE,
#'     expressionPvalCutoff = 0.01,
#'     hypoStringency = 0.3,
#'     coreCount = 8
#' )
step5OptimizeLinks <- function(
    TENETMultiAssayExperiment,
    hypermethGplusAnalysis = TRUE,
    hypomethGplusAnalysis = TRUE,
    expressionPvalCutoff = 0.05,
    hyperStringency = NA,
    hypoStringency = NA,
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

    ## Ensure the output data from the step3GetAnalysisZScores function
    ## are present in the MultiAssayExperiment object
    .ensureStepPresent(
        TENETMultiAssayExperiment, "step3GetAnalysisZScores"
    )

    ## Ensure the output data from the
    ## step4SelectMostSignificantLinksPerDNAMethylationSite function are present
    ## in the MultiAssayExperiment object
    .ensureStepPresent(
        TENETMultiAssayExperiment,
        "step4SelectMostSignificantLinksPerDNAMethylationSite"
    )

    ## Create an empty list to hold the step 5 results
    step5ResultsList <- list()

    ## Get the methylation values that match with expression values
    ## using the mapping data. This assumes the methylation and expression
    ## values share a clinical data match within the mapping.
    metToExpSampleConversion <- .createMetToExpSampleConversionVector(
        TENETMultiAssayExperiment
    )

    ## Get the names of the control and case samples in the methylation data
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

    ## Get the expression data for the control and case samples
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

    ## Get the minCaseCount value that was used in step 2
    minCaseCount <- unname(
        TENETMultiAssayExperiment@metadata$
            step2GetDifferentiallyMethylatedSites$minCaseCount
    )

    ## Depending on the analysis types selected, run optimization for RE DNA
    ## methylation sites of that type and add the resulting data frame to the
    ## step 5 metadata
    for (hyperHypo in analysisTypes) {
        step5ResultsList[[
            paste0(hyperHypo, "methGplusResults")
        ]] <- .optimizeQuadrantLinks(
            MAE = TENETMultiAssayExperiment,
            hyperHypo = hyperHypo,
            expressionDataControl = expressionDataControl,
            expressionDataCase = expressionDataCase,
            methylationSampleNamesControl = methylationSampleNamesControl,
            methylationSampleNamesCase = methylationSampleNamesCase,
            minCaseCount = minCaseCount,
            stringency = get(paste0(hyperHypo, "Stringency")),
            coreCount = coreCount,
            expressionPvalCutoff = expressionPvalCutoff,
            metToExpSampleConversion = metToExpSampleConversion
        )
    }

    ## Add the step 5 results list to the MultiAssayExperiment object
    TENETMultiAssayExperiment@metadata$step5OptimizeLinks <- step5ResultsList

    return(TENETMultiAssayExperiment)
}
