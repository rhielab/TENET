#' Generate boxplots comparing the methylation level of the specified RE DNA
#' methylation sites in case and control samples
#'
#' This function takes a vector of RE DNA methylation sites specified by the
#' user and generates boxplots displaying the methylation level of each of these
#' DNA methylation sites in the case compared to control samples, along with the
#' results of a Student's t-test comparing the methylation level between these
#' two groups.
#'
#' @param TENETMultiAssayExperiment Specify a MultiAssayExperiment object
#' containing a methylation SummarizedExperiment object, such as one created by
#' the TCGADownloader function.
#' @param DNAMethylationSites Supply a vector of RE DNA methylation site IDs
#' for which boxplots with the methylation of those RE DNA methylation sites
#' will be created.
#' @param coreCount Argument passed as the mc.cores argument to mcmapply. See
#' `?parallel::mcmapply` for more details. Defaults to 1.
#' @return Returns the MultiAssayExperiment object given as the
#' TENETMultiAssayExperiment argument with an additional list of data named
#' 'step7SelectedDNAMethylationSitesCaseVsControlBoxplots' in its metadata with
#' the output of this function, which contains boxplots showing the methylation
#' of the RE DNA methylation sites of interest in the case and control samples,
#' with Student's t-test p-values and the ID of the RE DNA methylation site in
#' the title.
#' @export
#'
#' @examplesIf interactive()
#' ## This example uses the example MultiAssayExperiment provided in the
#' ## TENET.ExperimentHub package to generate boxplots for several selected
#' ## RE DNA methylation sites, using one CPU core.
#'
#' ## Load the example TENET MultiAssayExperiment object
#' ## from the TENET.ExperimentHub package
#' exampleTENETMultiAssayExperiment <-
#'     TENET.ExperimentHub::exampleTENETMultiAssayExperiment()
#'
#' ## Use the example dataset to create RE DNA methylation site case vs. control
#' ## boxplots
#' returnValue <- step7SelectedDNAMethylationSitesCaseVsControlBoxplots(
#'     TENETMultiAssayExperiment = exampleTENETMultiAssayExperiment,
#'     DNAMethylationSites = c("cg03095778", "cg24011501", "cg12989041"),
#'     coreCount = 1
#' )
step7SelectedDNAMethylationSitesCaseVsControlBoxplots <- function(
    TENETMultiAssayExperiment,
    DNAMethylationSites,
    coreCount = 1) {
    ## Return an error message if the input MultiAssayExperiment is invalid
    .validateMultiAssayExperiment(TENETMultiAssayExperiment)

    if (missing(DNAMethylationSites)) {
        .stopNoCall("The DNAMethylationSites parameter must be specified.")
    }

    ## Get the methylation data
    methylationData <- .getExpOrMetSamplesOfType(
        TENETMultiAssayExperiment,
        "methylation"
    )

    ## Exclude RE DNA methylation sites that are not present in the methylation
    ## data
    methSiteList <- .excludeMissingMethylationSites(
        DNAMethylationSites, methylationData
    )

    ## Get the names and types of the methylation samples in the data
    sampleNames <- .getExpOrMetSamplesOfType(
        TENETMultiAssayExperiment,
        "methylation",
        namesOnly = TRUE
    )

    sampleTypes <- TENETMultiAssayExperiment@sampleMap[
        match(sampleNames, TENETMultiAssayExperiment@sampleMap$colname),
        "sampleType"
    ]

    ## Create a data frame containing group info for the case and control
    ## samples
    groupInfo <- data.frame(
        group = sampleNames,
        cluster = sampleTypes,
        stringsAsFactors = FALSE
    )

    ## Generate the plots for all RE DNA methylation sites of interest and add
    ## them to the MultiAssayExperiment
    TENETMultiAssayExperiment@metadata$
        step7SelectedDNAMethylationSitesCaseVsControlBoxplots <-
        parallel::mcmapply(
            FUN = .quadrantBoxplotFunction,
            geneOrMethSiteID = methSiteList,
            MoreArgs = list(
                expOrMet = "methylation",
                expOrMetData = methylationData,
                groupInfo = groupInfo
            ),
            mc.cores = coreCount,
            USE.NAMES = FALSE
        )

    return(TENETMultiAssayExperiment)
}
