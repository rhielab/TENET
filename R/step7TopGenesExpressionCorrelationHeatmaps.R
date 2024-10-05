## Internal functions used by step7TopGenesExpressionCorrelationHeatmaps

## Internal function to plot expression correlation heatmaps for the given
## quadrant
.generateQuadrantExpressionHeatmapsAndCorMatrix <- function(
    TENETMultiAssayExperiment,
    expressionDataCase,
    hyperHypo,
    geneOrTF,
    topGeneNumber,
    geneIDNameDF,
    quadrantResultsName) {
    ## Get the IDs of the top genes/TFs in this quadrant
    topQuadrantGeneOrTFIDs <- .getQuadrantTopGenesOrTFs(
        TENETMultiAssayExperiment, geneOrTF, hyperHypo,
        topGeneNumber
    )$geneID
    if (.isSingleNA(topQuadrantGeneOrTFIDs)) {
        return(NA)
    }

    ## Convert the gene IDs into gene names to display in the expression
    ## correlation heatmap and table
    topQuadrantGeneOrTFNames <- geneIDNameDF[
        topQuadrantGeneOrTFIDs, "geneName"
    ]
    names(topQuadrantGeneOrTFNames) <- topQuadrantGeneOrTFIDs

    ## Create a subset of the original case expression data to include just
    ## the genes of interest. We will also need to transpose the rows and
    ## columns so the genes are correlated, rather than the samples.
    topQuadrantGeneOrTFExpression <- t(
        expressionDataCase[topQuadrantGeneOrTFIDs, ]
    )

    ## Perform the correlation analysis to correlate the expression of each
    ## gene onto each other gene (including itself) in the case data and
    ## produce a correlation matrix
    topQuadrantGeneOrTFExpressionCorrelationMatrix <- stats::cor(
        topQuadrantGeneOrTFExpression
    )

    ## Create dendrograms that can be applied to both the rows and columns
    ## (since they are mirrored)
    topQuadrantGeneOrTFDist <- .step7HeatmapDistf(
        topQuadrantGeneOrTFExpressionCorrelationMatrix
    )

    topQuadrantGeneOrTFClust <- .step7HeatmapClustf(
        topQuadrantGeneOrTFDist
    )

    topQuadrantGeneOrTFDend <- stats::as.dendrogram(
        topQuadrantGeneOrTFClust
    )

    ## Create empty row and column labels (since we like the organization of
    ## heatmap.3 with the row and column labels, but we don't want to actually
    ## label them with anything for this analysis).
    topQuadrantGeneOrTFRowColorLabels <- rbind(rep(
        "white", nrow(topQuadrantGeneOrTFExpressionCorrelationMatrix)
    ))
    rownames(topQuadrantGeneOrTFRowColorLabels) <- ""

    topQuadrantGeneOrTFColColorLabels <- cbind(rep(
        "white", ncol(topQuadrantGeneOrTFExpressionCorrelationMatrix)
    ))
    colnames(topQuadrantGeneOrTFColColorLabels) <- ""

    ## Create the plot for the top genes/TFs of interest
    .newInvisibleRecordablePlot(width = 10, height = 10.5)
    .TENETInternalHeatmap3(
        x = topQuadrantGeneOrTFExpressionCorrelationMatrix,
        Rowv = topQuadrantGeneOrTFDend,
        Colv = topQuadrantGeneOrTFDend,
        RowSideColors = topQuadrantGeneOrTFRowColorLabels,
        ColSideColors = topQuadrantGeneOrTFColColorLabels,
        dendrogram = "col",
        labRow = NA,
        lmat = rbind(c(0, 0, 5), c(0, 0, 2), c(4, 1, 3), c(0, 0, 6)),
        lwid = c(0.02, 0.02, 2),
        lhei = c(0.6, 0.03, 2, 0.15),
        margins = c(6, 2),
        col = grDevices::colorRampPalette(c("blue", "white", "red"))(100),
        trace = "none",
        key = FALSE,
        main = NULL,
        ylab = paste0(
            "Top ", geneOrTF, "s linked to ",
            hyperHypo, "methylated G+ RE DNA methylation sites"
        ),
        xlab = NULL
    )

    ## Save the plot to an object
    quadrantHeatmap <- .recordTENETSavedSizePlot()

    ## Close the plot
    grDevices::dev.off()

    ## Create a copy of the expression correlation matrix in the same
    ## orientation as the heatmap sort
    sortedQuadrantExpressionCorrelationMatrix <-
        data.frame(
            topQuadrantGeneOrTFExpressionCorrelationMatrix[
                rev(labels(topQuadrantGeneOrTFDend)),
                labels(topQuadrantGeneOrTFDend)
            ]
        )

    ## Add gene names to the file and set them as the first column
    sortedQuadrantExpressionCorrelationMatrix$
        geneNames <- topQuadrantGeneOrTFNames[rownames(
        sortedQuadrantExpressionCorrelationMatrix
    )]

    sortedQuadrantExpressionCorrelationMatrix <-
        sortedQuadrantExpressionCorrelationMatrix[
            ,
            c(ncol(
                sortedQuadrantExpressionCorrelationMatrix
            ), seq_len(ncol(
                sortedQuadrantExpressionCorrelationMatrix
            ) - 1))
        ]

    ## Return a combined results list containing the heatmap and matrix
    return(list(
        "correlationMatrix" = sortedQuadrantExpressionCorrelationMatrix,
        "heatmap" = quadrantHeatmap
    ))
}

## Main step7TopGenesExpressionCorrelationHeatmaps function

#' Generate mirrored heatmaps displaying the correlation of the expression
#' values of the top genes/TFs
#'
#' This function takes the top genes/TFs for each analysis type by number of
#' linked RE DNA methylation sites identified by the
#' `step6DNAMethylationSitesPerGeneTabulation` function, up to the number
#' specified by the user, and generates heatmaps showing the correlation
#' R-values between the expression of the top genes/TFs, as well as data frames
#' with the R-values.
#'
#' @param TENETMultiAssayExperiment Specify a MultiAssayExperiment object
#' containing expression and methylation SummarizedExperiment objects, such as
#' one created by the TCGADownloader function. This MultiAssayExperiment object
#' should also contain the results from the
#' `step6DNAMethylationSitesPerGeneTabulation` function in its metadata.
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
#' @param hypermethGplusAnalysis Set to TRUE to create heatmaps and tables
#' showing expression correlation values for the top genes/TFs with the most
#' hypermethylated RE DNA methylation sites with G+ links. Defaults to TRUE.
#' @param hypomethGplusAnalysis Set to TRUE to create heatmaps and tables
#' showing expression correlation values for the top genes/TFs with the most
# ; hypomethylated RE DNA methylation sites with G+ links. Defaults to TRUE.
#' @param topGeneNumber Specify the number of top genes/TFs, based on the most
#' linked RE DNA methylation sites of a given analysis type, for which to
#' generate expression correlation heatmaps and tables. Defaults to 10.
#' @return Returns the MultiAssayExperiment object given as the
#' TENETMultiAssayExperiment argument with an additional list of information
#' named 'step7TopGenesExpressionCorrelationHeatmaps' in its metadata with
#' the output of this function. This list is subdivided into hypermethGplus or
#' hypomethGplus results as selected by the user, which are further subdivided
#' into lists with data for the top overall genes and for top TF genes only.
#' These contain the mirrored heatmaps displaying the expression correlation
#' values for the expression of top genes/TFs as well as data frames with names
#' and correlation values for each of those genes/TFs.
#' @export
#'
#' @examplesIf interactive()
#' ## This example uses the example MultiAssayExperiment provided in the
#' ## TENET.ExperimentHub package to create correlation heatmaps, and
#' ## corresponding data frames, for the top 10 genes/TFs by number of linked
#' ## hyper- or hypomethylated RE DNA methylation sites. Gene names will be
#' ## retrieved from the rowRanges of the 'expression' SummarizedExperiment
#' ## object in the example MultiAssayExperiment.
#'
#' ## Load the example TENET MultiAssayExperiment object
#' ## from the TENET.ExperimentHub package
#' exampleTENETMultiAssayExperiment <-
#'     TENET.ExperimentHub::exampleTENETMultiAssayExperiment()
#'
#' ## Use the example dataset to create expression correlation heatmaps
#' returnValue <- step7TopGenesExpressionCorrelationHeatmaps(
#'     TENETMultiAssayExperiment = exampleTENETMultiAssayExperiment
#' )
#'
#' ## This example uses the example MultiAssayExperiment provided in the
#' ## TENET.ExperimentHub package to create correlation heatmaps, and
#' ## corresponding data frames, for the top 5 genes/TFs by number of linked
#' ## hypomethylated RE DNA methylation sites only. Gene names will be retrieved
#' ## from the rowRanges of the 'expression' SummarizedExperiment object in the
#' ## example MultiAssayExperiment.
#'
#' ## Load the example TENET MultiAssayExperiment object
#' ## from the TENET.ExperimentHub package
#' exampleTENETMultiAssayExperiment <-
#'     TENET.ExperimentHub::exampleTENETMultiAssayExperiment()
#'
#' ## Use the example dataset to create expression correlation heatmaps
#' returnValue <- step7TopGenesExpressionCorrelationHeatmaps(
#'     TENETMultiAssayExperiment = exampleTENETMultiAssayExperiment,
#'     hypermethGplusAnalysis = FALSE,
#'     topGeneNumber = 5
#' )
step7TopGenesExpressionCorrelationHeatmaps <- function(
    TENETMultiAssayExperiment,
    geneAnnotationDataset = NA,
    hypermethGplusAnalysis = TRUE,
    hypomethGplusAnalysis = TRUE,
    topGeneNumber = 10) {
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

    ## Get the case expression dataset
    expressionDataCase <- .getExpOrMetSamplesOfType(
        TENETMultiAssayExperiment,
        "expression",
        "Case"
    )

    ## Create an empty list to hold the results from this step 7 function
    resultsList <- list()

    ## Generate heatmaps for the selected analysis types
    for (hyperHypo in analysisTypes) {
        quadrantResultsName <- paste0(hyperHypo, "methGplusResults")
        resultsList[[quadrantResultsName]] <- list()

        for (geneOrTF in c("Gene", "TF")) {
            resultsList[[quadrantResultsName]][[
                paste0("top", geneOrTF, "s")
            ]] <- .generateQuadrantExpressionHeatmapsAndCorMatrix(
                TENETMultiAssayExperiment,
                expressionDataCase,
                hyperHypo,
                geneOrTF,
                topGeneNumber,
                geneIDNameDF,
                quadrantResultsName
            )
        }
    }

    ## Add the results list to the MultiAssayExperiment
    TENETMultiAssayExperiment@metadata$
        step7TopGenesExpressionCorrelationHeatmaps <- resultsList

    return(TENETMultiAssayExperiment)
}
