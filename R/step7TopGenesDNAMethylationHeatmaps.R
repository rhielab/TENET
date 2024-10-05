## Internal functions used by step7TopGenesDNAMethylationHeatmaps

## Internal function to scale expression values across samples on a proportional
## scale from 1 to 200, ignoring 0 values
.rescaleZeroIgnored <- function(x) {
    ## Get the minimum and maximum of non-zero, non-NA values
    nonZeroVec <- x[x != 0]
    nonZeroVec <- nonZeroVec[!is.na(nonZeroVec)]
    usedMin <- min(nonZeroVec)
    usedMax <- max(nonZeroVec)

    place <- ((x - usedMin) / (usedMax - usedMin)) * 200

    returnValue <- ifelse(place <= 0, 1, place)

    return(ceiling(returnValue))
}

## Internal function to plot a methylation heatmap for the given quadrant
.plotQuadrantMethylationHeatmap <- function(
    TENETMultiAssayExperiment,
    expressionDataCase,
    quadrantMethylationDataCase,
    hyperHypo,
    geneOrTF,
    topGeneNumber,
    geneIDNameDF,
    quadrantSigLinkZScores,
    quadrantResultsName,
    metToExpSampleConversion) {
    ## Get the IDs of the top genes/TFs in this quadrant
    topQuadrantGeneOrTFIDs <- .getQuadrantTopGenesOrTFs(
        TENETMultiAssayExperiment, geneOrTF, hyperHypo, topGeneNumber
    )$geneID
    if (.isSingleNA(topQuadrantGeneOrTFIDs)) {
        return(NA)
    }

    ## Get the names for each of the top genes/TFs
    topQuadrantGeneOrTFNames <- geneIDNameDF[topQuadrantGeneOrTFIDs, "geneName"]

    ## Get the list of RE DNA methylation sites linked to the significant
    ## genes/TFs
    topGeneOrTFLinkedMethSites <- unique(
        quadrantSigLinkZScores[
            quadrantSigLinkZScores$geneID %in% topQuadrantGeneOrTFIDs,
            "DNAMethylationSiteID"
        ]
    )

    ## Create methylation datasets from the case samples of interest
    topGeneOrTFSitesMethylation <- as.matrix(
        quadrantMethylationDataCase[topGeneOrTFLinkedMethSites, ]
    )

    ## Create a subset of the original case expression data to include just
    ## the genes of interest. We will also need to transpose the rows and
    ## columns so the genes are correlated, rather than the samples.
    topQuadrantGeneOrTFExpression <- t(
        expressionDataCase[topQuadrantGeneOrTFIDs, ]
    )

    ## Rescale the expression values of the genes of interest
    topQuadrantGeneOrTFExpressionRescaled <- apply(
        topQuadrantGeneOrTFExpression, 2, .rescaleZeroIgnored
    )

    ## Get a vector of 200 heatmap colors ranging from dark blue to dark red
    heatmapColors <- matlab::jet.colors(200)

    ## Convert the rescaled values to colors
    topQuadrantGeneOrTFExpressionRescaledColor <- apply(
        topQuadrantGeneOrTFExpressionRescaled,
        2,
        function(x) {
            heatmapColors[x]
        }
    )
    rownames(topQuadrantGeneOrTFExpressionRescaledColor) <-
        rownames(topQuadrantGeneOrTFExpressionRescaled)

    ## Create empty row labels (since we like the organization of
    ## heatmap.3 with the row and column labels, but we don't want to actually
    ## label them with anything for this analysis).
    topQuadrantGeneOrTFRowColorLabels <- rbind(
        rep("white", nrow(topGeneOrTFSitesMethylation))
    )
    rownames(topQuadrantGeneOrTFRowColorLabels) <- ""

    ## Create column labels for the heatmap
    topQuadrantGeneOrTFColColorLabels <-
        topQuadrantGeneOrTFExpressionRescaledColor[
            metToExpSampleConversion[colnames(topGeneOrTFSitesMethylation)],
            rev(colnames(topQuadrantGeneOrTFExpressionRescaledColor))
        ]
    colnames(topQuadrantGeneOrTFColColorLabels) <- rev(
        topQuadrantGeneOrTFNames
    )

    ## Create a natural clustering of the heatmap
    topQuadrantGeneOrTFColDist <- .step7HeatmapDistf(
        t(topGeneOrTFSitesMethylation)
    )
    topQuadrantGeneOrTFColClust <- .step7HeatmapClustf(
        topQuadrantGeneOrTFColDist
    )
    topQuadrantGeneOrTFColDend <- stats::as.dendrogram(
        topQuadrantGeneOrTFColClust
    )

    ## Create a row clustering for the heatmap
    topQuadrantGeneOrTFRowDist <- .step7HeatmapDistf(
        topGeneOrTFSitesMethylation
    )
    topQuadrantGeneOrTFRowClust <- .step7HeatmapClustf(
        topQuadrantGeneOrTFRowDist
    )
    topQuadrantGeneOrTFRowDend <- stats::as.dendrogram(
        topQuadrantGeneOrTFRowClust
    )

    ## Create the plot for the top genes/TFs of interest
    .newInvisibleRecordablePlot(width = 20, height = 14)
    .TENETInternalHeatmap3(
        x = topGeneOrTFSitesMethylation,
        Rowv = topQuadrantGeneOrTFRowDend,
        Colv = topQuadrantGeneOrTFColDend,
        RowSideColors = topQuadrantGeneOrTFRowColorLabels,
        ColSideColors = topQuadrantGeneOrTFColColorLabels,
        dendrogram = "col",
        labCol = NA,
        labRow = NA,
        lmat = rbind(c(0, 0, 5), c(0, 0, 2), c(4, 1, 3), c(0, 0, 6)),
        lwid = c(0.25, 0.02, 2),
        lhei = c(0.4, 0.6, 2, 0.001),
        margins = c(2, 2),
        col = heatmapColors,
        trace = "none",
        key = FALSE,
        main = NULL,
        ylab = paste0(
            "Regulatory element DNA methylation sites linked to top ",
            ifelse(geneOrTF == "TF", "transcription factors", "genes")
        ),
        xlab = "Case samples"
    )

    ## Save the plot to an object
    quadrantHeatmap <- .recordTENETSavedSizePlot()

    ## Close the plot
    grDevices::dev.off()

    ## Return the heatmap as an object
    return(quadrantHeatmap)
}

## Main step7TopGenesDNAMethylationHeatmaps function

#' Generate heatmaps displaying the methylation level of all RE DNA methylation
#' sites linked to the top genes/transcription factors, along with the
#' expression of those genes in the column headers, in the case samples within
#' the supplied MultiAssayExperiment object
#'
#' This function takes the top genes/transcription factors (TFs) for each
#' analysis type by number of linked RE DNA methylation sites identified by the
#' `step6DNAMethylationSitesPerGeneTabulation` function, up to the number
#' specified by the user, and generates heatmaps displaying the methylation
#' level of the unique RE DNA methylation sites linked to any of those genes,
#' along with the expression of those genes in the case samples only.
#'
#' @param TENETMultiAssayExperiment Specify a MultiAssayExperiment object
#' containing expression and methylation SummarizedExperiment objects, such as
#' one created by the TCGADownloader function. This MultiAssayExperiment object
#' should also contain the results from the
#' step2GetDifferentiallyMethylatedSites, `step5OptimizeLinks`, and
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
#' @param hypermethGplusAnalysis Set to TRUE to create heatmaps showing DNA
#' methylation levels of RE DNA methylation sites linked to the top genes/TFs
#' with the most hypermethylated RE DNA methylation sites with G+ links.
#' Defaults to TRUE.
#' @param hypomethGplusAnalysis Set to TRUE to create heatmaps showing DNA
#' methylation levels of RE DNA methylation sites linked to the top genes/TFs
#' with the most hypomethylated RE DNA methylation sites with G+ links.
#' Defaults to TRUE.
#' @param topGeneNumber Specify the number of top genes/TFs, based on the most
#' linked RE DNA methylation sites of a given analysis type, for which to
#' generate heatmaps with their linked RE DNA methylation sites' methylation
#' levels. Defaults to 10.
#' @return Returns the MultiAssayExperiment object given as the
#' TENETMultiAssayExperiment argument with an additional list of information
#' named 'step7TopGenesDNAMethylationHeatmaps' in its metadata with the output
#' of this function. This list is subdivided into hypermethGplus or
#' hypomethGplus results as selected by the user, which are further subdivided
#' into lists with data for the top overall genes and for top TF genes only.
#' Each of these contains a single heatmap, with the expression of the top
#' genes/TFs in the column headers and the methylation of their unique linked
#' RE DNA methylation sites in the body of the heatmaps.
#' @export
#'
#' @examplesIf interactive()
#' ## This example uses the example MultiAssayExperiment provided in the
#' ## TENET.ExperimentHub package to create heatmaps for the top 10 genes/TFs,
#' ## by number of linked hyper- or hypomethylated RE DNA methylation sites, as
#' ## well as the unique RE DNA methylation sites linked to those 10 genes.
#' ## Gene names will be retrieved from the rowRanges of the 'expression'
#' ## SummarizedExperiment object in the example MultiAssayExperiment.
#'
#' ## Load the example TENET MultiAssayExperiment object
#' ## from the TENET.ExperimentHub package
#' exampleTENETMultiAssayExperiment <-
#'     TENET.ExperimentHub::exampleTENETMultiAssayExperiment()
#'
#' ## Use the example dataset to create methylation heatmaps
#' returnValue <- step7TopGenesDNAMethylationHeatmaps(
#'     TENETMultiAssayExperiment = exampleTENETMultiAssayExperiment
#' )
#'
#' ## This example uses the example MultiAssayExperiment provided in the
#' ## TENET.ExperimentHub package to create heatmaps for only the top 5
#' ## genes/TFs, by number of linked hypomethylated RE DNA methylation sites
#' ## only, as well as the unique RE DNA methylation sites linked to those 5
#' ## genes/TFs. Gene names will be retrieved from the rowRanges of the
#' ## 'expression' SummarizedExperiment object in the example
#' ## MultiAssayExperiment.
#'
#' ## Load the example TENET MultiAssayExperiment object
#' ## from the TENET.ExperimentHub package
#' exampleTENETMultiAssayExperiment <-
#'     TENET.ExperimentHub::exampleTENETMultiAssayExperiment()
#'
#' ## Use the example dataset to create methylation heatmaps
#' returnValue <- step7TopGenesDNAMethylationHeatmaps(
#'     TENETMultiAssayExperiment = exampleTENETMultiAssayExperiment,
#'     hypermethGplusAnalysis = FALSE,
#'     topGeneNumber = 5
#' )
step7TopGenesDNAMethylationHeatmaps <- function(
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

    ## Ensure the output data from the step2GetDifferentiallyMethylatedSites
    ## function are present in the MultiAssayExperiment object
    .ensureStepPresent(
        TENETMultiAssayExperiment, "step2GetDifferentiallyMethylatedSites"
    )

    ## Get gene IDs and names from the MAE, or gene annotation dataset if
    ## provided
    geneIDNameDF <- .getGeneIDsAndNames(
        TENETMultiAssayExperiment, geneAnnotationDataset
    )

    ## Get the expression/methylation names of the case samples
    expressionSampleNamesCase <- .getExpOrMetSamplesOfType(
        TENETMultiAssayExperiment,
        "expression",
        "Case",
        namesOnly = TRUE
    )

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

    ## Get the methylation values that match with expression values
    ## using the mapping data. This assumes the methylation and expression
    ## values share a clinical data match within the mapping.
    metToExpSampleConversion <- .createMetToExpSampleConversionVector(
        TENETMultiAssayExperiment
    )

    ## Create an empty list to hold the results from this step 7 function
    resultsList <- list()

    ## Generate heatmaps for the selected analysis types
    for (hyperHypo in analysisTypes) {
        quadrantResultsName <- paste0(hyperHypo, "methGplusResults")

        ## Ensure the quadrant's results are present in step 5
        .ensureStepPresent(
            TENETMultiAssayExperiment,
            stepName = "step5OptimizeLinks",
            substepName = quadrantResultsName
        )

        resultsList[[quadrantResultsName]] <- list()

        ## Load the quadrant's significant links from step 5
        quadrantSigLinkZScores <- TENETMultiAssayExperiment@metadata$
            step5OptimizeLinks[[quadrantResultsName]]

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

        for (geneOrTF in c("Gene", "TF")) {
            ## Plot the heatmap and add it to the results list
            resultsList[[quadrantResultsName]][[
                paste0("top", geneOrTF, "s")
            ]] <- .plotQuadrantMethylationHeatmap(
                TENETMultiAssayExperiment,
                expressionDataCase,
                quadrantMethylationDataCase,
                hyperHypo,
                geneOrTF,
                topGeneNumber,
                geneIDNameDF,
                quadrantSigLinkZScores,
                quadrantResultsName,
                metToExpSampleConversion
            )
        }
    }

    ## Add the results list to the MultiAssayExperiment
    TENETMultiAssayExperiment@metadata$
        step7TopGenesDNAMethylationHeatmaps <- resultsList

    return(TENETMultiAssayExperiment)
}
