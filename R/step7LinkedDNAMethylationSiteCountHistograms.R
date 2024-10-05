## Internal functions used by step7LinkedDNAMethylationSiteCountHistograms

## Internal function to generate a histogram for the given quadrant
.generateQuadrantHistogram <- function(
    hyperHypo, geneOrTF, TENETMultiAssayExperiment, quadrantResultsName) {
    ## Ensure the quadrant's results are present in step 6
    .ensureStepPresent(
        TENETMultiAssayExperiment,
        stepName = "step6DNAMethylationSitesPerGeneTabulation",
        substepName = quadrantResultsName
    )

    ## Load the quadrant's link counts from step 6 (NA means all genes/TFs)
    quadrantAllGeneOrTFFreq <- .getQuadrantTopGenesOrTFs(
        TENETMultiAssayExperiment, geneOrTF, hyperHypo, NA
    )
    if (.isSingleNA(quadrantAllGeneOrTFFreq$geneID)) {
        return(NA)
    }

    ## Create a histogram that shows the number of genes/TFs with a given number
    ## of the quadrant's links linked to them
    linksDescription <- paste0(
        "linked RE DNA methylation sites per ",
        ifelse(geneOrTF == "TF", "TF ", ""),
        "gene"
    )

    ## Using get here because R CMD check does not understand lazy evaluation
    quadrantHistogram <- ggplot2::ggplot(
        quadrantAllGeneOrTFFreq, ggplot2::aes(x = get("count"))
    ) +
        ggplot2::geom_histogram(
            color = "black",
            fill = "red",
            binwidth = ceiling(max(quadrantAllGeneOrTFFreq$count) / 200)
        ) +
        ggplot2::ggtitle(paste0(
            "Histogram of ", hyperHypo, "methylated G+ ", linksDescription
        )) +
        ggplot2::xlab(paste("Number of", linksDescription)) +
        ggplot2::ylab("Frequency") +
        ggplot2::theme_bw() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(hjust = 0.5, size = 20),
            panel.border = ggplot2::element_rect(
                color = "black", fill = NA, linewidth = 1
            ),
            plot.background = ggplot2::element_rect(fill = "white"),
            axis.title.x = ggplot2::element_text(size = 20),
            axis.title.y = ggplot2::element_text(size = 20),
            axis.text.x = ggplot2::element_text(size = 18, color = "black"),
            axis.text.y = ggplot2::element_text(size = 16, color = "black"),
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank()
        )

    ## To display this plot correctly, the plot area should be at least 10
    ## inches wide by 7 inches tall
    quadrantHistogram <- TENETSavedSizePlot(
        quadrantHistogram,
        width = 10, height = 7
    )

    ## Return the plot as a ggplot2 object
    return(quadrantHistogram)
}

## Main step7LinkedDNAMethylationSiteCountHistograms function

#' Create histograms displaying the number of genes or transcription factors
#' linked to a given number of RE DNA methylation sites
#'
#' This function generates histograms displaying the number of all genes,
#' or transcription factor (TF) genes only, with links to a given number of
#' regulatory element DNA methylation sites from both of the hyper- or
#' hypomethylated G+ analysis quadrants, as selected by the user.
#'
#' @param TENETMultiAssayExperiment Specify a MultiAssayExperiment object
#' containing expression and methylation SummarizedExperiment objects, such as
#' one created by the TCGADownloader function. This MultiAssayExperiment object
#' should also contain the results from the
#' step6DNAMethylationSitesPerGeneTabulation function in its metadata.
#' @param hypermethGplusAnalysis Set to TRUE to create histograms of genes/TFs
#' linked to hypermethylated RE DNA methylation sites with G+ links. Defaults
#' to TRUE.
#' @param hypomethGplusAnalysis Set to TRUE to create histograms of genes/TFs
#' linked to hypomethylated RE DNA methylation sites with G+ links. Defaults to
#' TRUE.
#' @return Returns the MultiAssayExperiment object given as the
#' TENETMultiAssayExperiment argument with an additional list of information
#' named 'step7LinkedDNAMethylationSiteCountHistograms' in its metadata with
#' the output of this function. This list is subdivided into hypermethGplus or
#' hypomethGplus results as selected by the user. Each of these contains a
#' histogram displaying the number of all genes, or TFs only, linked to a given
#' number of RE DNA methylation sites in each of the selected analysis
#' quadrants.
#' @export
#'
#' @examplesIf interactive()
#' ## This example uses the example MultiAssayExperiment provided in the
#' ## TENET.ExperimentHub package to create histograms displaying counts of
#' ## genes/TFs by the number of linked hyper- or hypomethylated G+ RE DNA
#' ## methylation sites.
#'
#' ## Load the example TENET MultiAssayExperiment object
#' ## from the TENET.ExperimentHub package
#' exampleTENETMultiAssayExperiment <-
#'     TENET.ExperimentHub::exampleTENETMultiAssayExperiment()
#'
#' ## Use the example dataset to create the RE DNA methylation site count
#' ## histograms
#' returnValue <- step7LinkedDNAMethylationSiteCountHistograms(
#'     TENETMultiAssayExperiment = exampleTENETMultiAssayExperiment
#' )
#'
#' ## This example also uses the example MultiAssayExperiment provided in the
#' ## TENET.ExperimentHub package, but only creates a histogram of counts of
#' ## genes/TFs by the number of linked hypomethylated G+ RE DNA methylation
#' ## sites.
#'
#' ## Load the example TENET MultiAssayExperiment object
#' ## from the TENET.ExperimentHub package
#' exampleTENETMultiAssayExperiment <-
#'     TENET.ExperimentHub::exampleTENETMultiAssayExperiment()
#'
#' ## Use the example dataset to create the RE DNA methylation site count
#' ## histograms
#' returnValue <- step7LinkedDNAMethylationSiteCountHistograms(
#'     TENETMultiAssayExperiment = exampleTENETMultiAssayExperiment,
#'     hypermethGplusAnalysis = FALSE
#' )
step7LinkedDNAMethylationSiteCountHistograms <- function(
    TENETMultiAssayExperiment,
    hypermethGplusAnalysis = TRUE,
    hypomethGplusAnalysis = TRUE) {
    ## Validate the analysis types and get a vector of the ones selected
    analysisTypes <- .validateAnalysisTypes(
        hypermethGplusAnalysis, hypomethGplusAnalysis
    )

    ## Return an error message if the input MultiAssayExperiment is invalid
    .validateMultiAssayExperiment(TENETMultiAssayExperiment)

    ## Create an empty list to hold the results from this step 7 function
    resultsList <- list()

    ## Generate histograms for the selected analysis types
    for (hyperHypo in analysisTypes) {
        quadrantResultsName <- paste0(hyperHypo, "methGplusResults")
        resultsList[[quadrantResultsName]] <- list()

        for (geneOrTF in c("Gene", "TF")) {
            resultsList[[quadrantResultsName]][[
                paste0("top", geneOrTF, "s")
            ]] <- .generateQuadrantHistogram(
                hyperHypo,
                geneOrTF,
                TENETMultiAssayExperiment,
                quadrantResultsName
            )
        }
    }

    ## Add the results list to the MultiAssayExperiment
    TENETMultiAssayExperiment@metadata$
        step7LinkedDNAMethylationSiteCountHistograms <- resultsList

    return(TENETMultiAssayExperiment)
}
