#' Create boxplots comparing the expression level of the top
#' genes/transcription factors in case and control samples
#'
#' This function takes the top genes/transcription factors (TFs) for each
#' analysis type by number of linked RE DNA methylation sites identified by the
#' `step6DNAMethylationSitesPerGeneTabulation` function, up to the number
#' specified by the user, and generates boxplots displaying the expression
#' level of each of these genes in the case compared to control samples, along
#' with the results of a Student's t-test comparing the expression level
#' between these two groups.
#'
#' @param TENETMultiAssayExperiment Specify a MultiAssayExperiment object
#' containing expression and methylation SummarizedExperiment objects, such as
#' one created by the TCGADownloader function. This MultiAssayExperiment object
#' should also contain the results from the `step5OptimizeLinks` and
#' `step6DNAMethylationSitesPerGeneTabulation` functions in its metadata.
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
#' @param hypermethGplusAnalysis Set to TRUE to create expression boxplots for
#' the top genes/TFs with the most hypermethylated RE DNA methylation sites
#' with G+ links. Defaults to TRUE.
#' @param hypomethGplusAnalysis Set to TRUE to create expresion boxplots for
#' the top genes/TFs with the most hypomethylated RE DNA methylation sites with
#' G+ links. Defaults to TRUE.
#' @param topGeneNumber Specify the number of top genes/TFs, based on the most
#' linked RE DNA methylation sites of a given analysis type, for which to
#' generate expression boxplots. Defaults to 10.
#' @param coreCount Argument passed as the mc.cores argument to mcmapply. See
#' `?parallel::mcmapply` for more details. Defaults to 1.
#' @return Returns the MultiAssayExperiment object given as the
#' TENETMultiAssayExperiment argument with an additional list of information
#' named 'step7TopGenesCaseVsControlExpressionBoxplots' in its metadata
#' with the output of this function. This list is subdivided into
#' hypermethGplus or hypomethGplus results as selected by the user, which are
#' further subdivided into lists with plots for the top overall genes, and for
#' top TF genes only. These contain boxplots showing the expression of the genes
#' of interest in the case and control samples, with Student's t-test p-values
#' and the name and ID of the gene in the title.
#' @export
#'
#' @examplesIf interactive()
#' ## This example uses the example MultiAssayExperiment provided in the
#' ## TENET.ExperimentHub package to create expression boxplots in case vs.
#' ## control samples for the top 10 genes/TFs, by number of linked hyper- or
#' ## hypomethylated RE DNA methylation sites. Gene names will be retrieved
#' ## from the rowRanges of the 'expression' SummarizedExperiment object in the
#' ## example MultiAssayExperiment. The analysis will be performed using one
#' ## CPU core.
#'
#' ## Load the example TENET MultiAssayExperiment object
#' ## from the TENET.ExperimentHub package
#' exampleTENETMultiAssayExperiment <-
#'     TENET.ExperimentHub::exampleTENETMultiAssayExperiment()
#'
#' ## Use the example dataset to create expression case vs. control boxplots
#' returnValue <- step7TopGenesCaseVsControlExpressionBoxplots(
#'     TENETMultiAssayExperiment = exampleTENETMultiAssayExperiment
#' )
#'
#' ## This example uses the example MultiAssayExperiment provided in the
#' ## TENET.ExperimentHub package to create expression boxplots in case vs.
#' ## control samples for the top 5 genes/TFs, by number of linked
#' ## hypomethylated RE DNA methylation sites only. Gene names will be
#' ## retrieved from the rowRanges of the 'expression' SummarizedExperiment
#' ## object in the example MultiAssayExperiment. The analysis will be
#' ## performed using 8 CPU cores.
#'
#' ## Load the example TENET MultiAssayExperiment object
#' ## from the TENET.ExperimentHub package
#' exampleTENETMultiAssayExperiment <-
#'     TENET.ExperimentHub::exampleTENETMultiAssayExperiment()
#'
#' ## Use the example dataset to create expression case vs. control boxplots
#' returnValue <- step7TopGenesCaseVsControlExpressionBoxplots(
#'     TENETMultiAssayExperiment = exampleTENETMultiAssayExperiment,
#'     hypermethGplusAnalysis = FALSE,
#'     topGeneNumber = 5,
#'     coreCount = 8
#' )
step7TopGenesCaseVsControlExpressionBoxplots <- function(
    TENETMultiAssayExperiment,
    geneAnnotationDataset = NA,
    hypermethGplusAnalysis = TRUE,
    hypomethGplusAnalysis = TRUE,
    topGeneNumber = 10,
    coreCount = 1) {
    ## Validate the analysis types and get a vector of the ones selected
    analysisTypes <- .validateAnalysisTypes(
        hypermethGplusAnalysis, hypomethGplusAnalysis
    )

    ## Return an error message if the input MultiAssayExperiment is invalid
    .validateMultiAssayExperiment(
        TENETMultiAssayExperiment,
        needGeneName = is.na(geneAnnotationDataset)
    )

    ## Get gene IDs and names from the MAE, or gene annotation dataset if
    ## provided
    geneIDNameDF <- .getGeneIDsAndNames(
        TENETMultiAssayExperiment, geneAnnotationDataset
    )

    ## Get the expression data
    expressionData <- .getExpOrMetSamplesOfType(
        TENETMultiAssayExperiment,
        "expression"
    )

    ## Get the names and types of the samples in the data
    sampleNames <- .getExpOrMetSamplesOfType(
        TENETMultiAssayExperiment,
        "expression",
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

    ## Create an empty list to hold the results from this step 7 function
    resultsList <- list()

    for (hyperHypo in analysisTypes) {
        quadrantResultsName <- paste0(hyperHypo, "methGplusResults")

        for (geneOrTF in c("Gene", "TF")) {
            ## Get the IDs of the top genes/TFs in this quadrant
            topQuadrantGeneOrTFIDs <- .getQuadrantTopGenesOrTFs(
                TENETMultiAssayExperiment, geneOrTF, hyperHypo,
                topGeneNumber
            )$geneID
            if (.isSingleNA(topQuadrantGeneOrTFIDs)) {
                resultsList[[
                    quadrantResultsName
                ]][[
                    paste0("top", geneOrTF, "s")
                ]] <- NA
                next
            }

            ## Generate the plots for the genes of interest and add them to the
            ## results list
            resultsList[[
                quadrantResultsName
            ]][[
                paste0("top", geneOrTF, "s")
            ]] <- parallel::mcmapply(
                FUN = .quadrantBoxplotFunction,
                geneOrMethSiteID = topQuadrantGeneOrTFIDs,
                MoreArgs = list(
                    expOrMet = "expression",
                    expOrMetData = expressionData,
                    geneIDNameDF = geneIDNameDF,
                    groupInfo = groupInfo
                ),
                mc.cores = coreCount,
                USE.NAMES = FALSE
            )
        }
    }

    ## Add the results list to the MultiAssayExperiment
    TENETMultiAssayExperiment@metadata$
        step7TopGenesCaseVsControlExpressionBoxplots <- resultsList

    return(TENETMultiAssayExperiment)
}
