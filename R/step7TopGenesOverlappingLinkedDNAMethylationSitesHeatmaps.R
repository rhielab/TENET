## Internal functions used by
## step7TopGenesOverlappingLinkedDNAMethylationSitesHeatmaps

## Internal clustering function (unique to this heatmap function)
.step7HeatmapDistfBinary <- function(d) {
    stats::dist(d, method = "binary")
}

## Internal function to plot overlapping linked RE DNA methylation site heatmaps
## for the given quadrant
.plotQuadrantLinkedDNAMethylationSitesHeatmaps <- function(
    TENETMultiAssayExperiment,
    hyperHypo,
    geneOrTF,
    topGeneNumber,
    geneIDNameDF,
    quadrantSigLinkZScores,
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
    topQuadrantGeneOrTFNames <- geneIDNameDF[topQuadrantGeneOrTFIDs, "geneName"]
    names(topQuadrantGeneOrTFNames) <- topQuadrantGeneOrTFIDs

    ## Create a list of vectors, with each vector containing the names of the
    ## RE DNA methylation sites linked to each TF/gene using a for loop, by
    ## first creating an empty list, then looping through each of the gene IDs
    ## to grab RE DNA methylation sites linked to each and add them to the list

    ## For each gene, list the RE DNA methylation sites linked to it
    topQuadrantGeneOrTFLinkedMethSiteList <- list()

    for (i in seq_along(topQuadrantGeneOrTFIDs)) {
        topQuadrantGeneOrTFLinkedMethSiteList[[i]] <-
            quadrantSigLinkZScores[
                quadrantSigLinkZScores$geneID ==
                    topQuadrantGeneOrTFIDs[i],
                "DNAMethylationSiteID"
            ]
    }

    ## Then identify all unique RE DNA methylation sites linked to at least one
    ## of the top genes/TFs
    topQuadrantGeneOrTFLinkedUniqueMethSites <- unique(
        c(unlist(topQuadrantGeneOrTFLinkedMethSiteList))
    )

    ## For each of the genes/TFs of interest, we need to take the unique RE DNA
    ## methylation sites and identify which of them are linked to each of the
    ## genes/TFs of interest, creating a list of vectors of 1s for TRUE and 0s
    ## for FALSE
    topQuadrantGeneOrTFUniqueMethSitesLinkedList <- list()

    for (i in seq_along(topQuadrantGeneOrTFLinkedMethSiteList)) {
        ## For each unique DNA methylation site, identify if it is linked to a
        ## given gene
        TFVector <- topQuadrantGeneOrTFLinkedUniqueMethSites %in%
            topQuadrantGeneOrTFLinkedMethSiteList[[i]]

        ## Convert the TRUE and FALSEs into 1 and 0s and add them to the list
        topQuadrantGeneOrTFUniqueMethSitesLinkedList[[i]] <- ifelse(
            TFVector == TRUE,
            1,
            0
        )
    }

    ## Take the uniqueMethSitesLinkedList and collapse it into a data frame with
    ## the genes in the rows and the RE DNA methylation sites in the columns
    topQuadrantGeneOrTFUniqueMethSitesLinkedDF <- as.data.frame(
        do.call(rbind, topQuadrantGeneOrTFUniqueMethSitesLinkedList)
    )

    ## Add the gene names to the rows, and the methylation site IDs as the
    ## column names
    rownames(topQuadrantGeneOrTFUniqueMethSitesLinkedDF) <-
        topQuadrantGeneOrTFNames

    colnames(topQuadrantGeneOrTFUniqueMethSitesLinkedDF) <-
        topQuadrantGeneOrTFLinkedUniqueMethSites

    ## Create a column dendrogram
    topQuadrantGeneOrTFColDist <- .step7HeatmapDistfBinary(
        t(topQuadrantGeneOrTFUniqueMethSitesLinkedDF)
    )
    topQuadrantGeneOrTFColClust <- .step7HeatmapClustf(
        topQuadrantGeneOrTFColDist
    )
    topQuadrantGeneOrTFColDend <- stats::as.dendrogram(
        topQuadrantGeneOrTFColClust
    )

    ## Create a row dendrogram
    topQuadrantGeneOrTFRowDist <- .step7HeatmapDistfBinary(
        topQuadrantGeneOrTFUniqueMethSitesLinkedDF
    )
    topQuadrantGeneOrTFRowClust <- .step7HeatmapClustf(
        topQuadrantGeneOrTFRowDist
    )
    topQuadrantGeneOrTFRowDend <- stats::as.dendrogram(
        topQuadrantGeneOrTFRowClust
    )

    ## Create blank row labels for the heatmap
    topQuadrantGeneOrTFColColorLabels <- cbind(
        rep("white", ncol(topQuadrantGeneOrTFUniqueMethSitesLinkedDF))
    )
    colnames(topQuadrantGeneOrTFColColorLabels) <- ""

    ## Create blank row labels for the heatmap
    topQuadrantGeneOrTFRowColorLabels <- rbind(
        rep("white", nrow(topQuadrantGeneOrTFUniqueMethSitesLinkedDF))
    )
    rownames(topQuadrantGeneOrTFRowColorLabels) <- ""

    ## Create the plot for the top genes/TFs of interest
    .newInvisibleRecordablePlot(width = 20, height = 14)
    .TENETInternalHeatmap3(
        x = topQuadrantGeneOrTFUniqueMethSitesLinkedDF,
        Rowv = topQuadrantGeneOrTFRowDend,
        Colv = topQuadrantGeneOrTFColDend,
        RowSideColors = topQuadrantGeneOrTFRowColorLabels,
        ColSideColors = topQuadrantGeneOrTFColColorLabels,
        dendrogram = "col",
        labCol = NA,
        labRow = rownames(topQuadrantGeneOrTFUniqueMethSitesLinkedDF),
        lmat = rbind(c(0, 0, 5), c(0, 0, 2), c(4, 1, 3), c(0, 0, 6)),
        lwid = c(0.02, 0.02, 2),
        lhei = c(0.6, 0.03, 2, 0.15),
        margins = c(2, 6),
        col = grDevices::colorRampPalette(c("white", "black"))(60),
        trace = "none",
        key = FALSE,
        main = NULL,
        xlab = paste0(
            "Unique RE DNA methylation sites linked to top ",
            hyperHypo, "methylated G+ ", geneOrTF, "s"
        )
    )

    ## Save the plot to an object
    quadrantHeatmap <- .recordTENETSavedSizePlot()

    ## Close the plot
    grDevices::dev.off()

    ## Also save the table indicating the RE DNA methylation sites and genes/TFs
    ## and which are linked to each other (1=link, 0=no link). It will also be
    ## formatted to fit the ordering displayed in the heatmap.
    topQuadrantGeneOrTFUniqueMethSitesLinkedOrderedDF <-
        topQuadrantGeneOrTFUniqueMethSitesLinkedDF[
            rev(labels(topQuadrantGeneOrTFRowDend)),
            labels(topQuadrantGeneOrTFColDend)
        ]

    ## Return a combined results list containing the heatmap and link table
    return(list(
        "linkTable" = topQuadrantGeneOrTFUniqueMethSitesLinkedOrderedDF,
        "heatmap" = quadrantHeatmap
    ))
}

## Main step7TopGenesOverlappingLinkedDNAMethylationSitesHeatmaps function

#' Generate binary heatmaps displaying which of the top genes/transcription
#' factors share links with each of the unique regulatory element DNA
#' methylation sites linked to at least one top gene/TF
#'
#' This function takes the top genes/TFs for each analysis type by number of
#' linked RE DNA methylation sites identified by the
#' `step6DNAMethylationSitesPerGeneTabulation` function, up to the number
#' specified by the user, and identifies the unique RE DNA methylation sites
#' linked to them, then generates two-color binary heatmaps displaying which of
#' the top genes/TFs the RE DNA methylation sites are linked to, as well as
#' data frames with that information.
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
#' @param hypermethGplusAnalysis Set to TRUE to create heatmaps and tables
#' showing the linked RE DNA methylation sites for the top genes/TFs with the
#' most hypermethylated RE DNA methylation sites with G+ links. Defaults to
#' TRUE.
#' @param hypomethGplusAnalysis Set to TRUE to create heatmaps and tables
#' showing the linked RE DNA methylation sites for the top genes/TFs with the
#' most hypomethylated RE DNA methylation sites with G+ links. Defaults to TRUE.
#' @param topGeneNumber Specify the number of top genes/TFs, based on the most
#' linked RE DNA methylation sites of a given analysis type, for which to
#' generate linked RE DNA methylation site heatmaps and tables. Defaults to 10.
#' @return Returns the MultiAssayExperiment object given as the
#' TENETMultiAssayExperiment argument with an additional list of information
#' named 'step7TopGenesOverlappingLinkedDNAMethylationSitesHeatmaps' in its
#' metadata with the output of this function. This list is subdivided into
#' hypermethGplus or hypomethGplus results as selected by the user, which are
#' further subdivided into lists with data for the top overall genes and for
#' top TF genes only. These contain the binary heatmaps displaying the unique
#' RE DNA methylation sites linked to the top genes/TFs in the columns and each
#' of the top genes/TFs in the rows, with black indicating the given RE DNA
#' methylation site is linked to the given gene/TF. Data frames are also
#' included with this same data, with 1s indicating an RE DNA methylation site
#' is linked to the gene/TF.
#' @export
#'
#' @examplesIf interactive()
#' ## This example uses the example MultiAssayExperiment provided in the
#' ## TENET.ExperimentHub package to create overlap heatmaps, and corresponding
#' ## data frames, for the top 10 genes/TFs by number of linked hyper- or
#' ## hypomethylated RE DNA methylation sites. Gene names will be retrieved
#' ## from the rowRanges of the 'expression' SummarizedExperiment object in the
#' ## example MultiAssayExperiment.
#'
#' ## Load the example TENET MultiAssayExperiment object
#' ## from the TENET.ExperimentHub package
#' exampleTENETMultiAssayExperiment <-
#'     TENET.ExperimentHub::exampleTENETMultiAssayExperiment()
#'
#' ## Use the example dataset to create overlap heatmaps
#' returnValue <- step7TopGenesOverlappingLinkedDNAMethylationSitesHeatmaps(
#'     TENETMultiAssayExperiment = exampleTENETMultiAssayExperiment
#' )
#'
#' ## This example uses the example MultiAssayExperiment provided in the
#' ## TENET.ExperimentHub package to create overlap heatmaps, and corresponding
#' ## data frames, for the top 5 genes/TFs by number of linked hypomethylated
#' ## RE DNA methylation sites only. Gene names will be retrieved from the
#' ## rowRanges of the 'expression' SummarizedExperiment object in the example
#' ## MultiAssayExperiment.
#'
#' ## Load the example TENET MultiAssayExperiment object
#' ## from the TENET.ExperimentHub package
#' exampleTENETMultiAssayExperiment <-
#'     TENET.ExperimentHub::exampleTENETMultiAssayExperiment()
#'
#' ## Use the example dataset to create overlap heatmaps
#' returnValue <- step7TopGenesOverlappingLinkedDNAMethylationSitesHeatmaps(
#'     TENETMultiAssayExperiment = exampleTENETMultiAssayExperiment,
#'     hypermethGplusAnalysis = FALSE,
#'     topGeneNumber = 5
#' )
step7TopGenesOverlappingLinkedDNAMethylationSitesHeatmaps <- function(
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

    ## Create an empty list to hold the results from this step 7 function
    resultsList <- list()

    ## Generate heatmaps for the selected analysis types
    for (hyperHypo in analysisTypes) {
        quadrantResultsName <- paste0(hyperHypo, "methGplusResults")
        resultsList[[quadrantResultsName]] <- list()

        ## Ensure the quadrant's results are present in step 5
        .ensureStepPresent(
            TENETMultiAssayExperiment,
            stepName = "step5OptimizeLinks",
            substepName = quadrantResultsName
        )

        ## Load the quadrant's significant links from step 5
        quadrantSigLinkZScores <- TENETMultiAssayExperiment@metadata$
            step5OptimizeLinks[[quadrantResultsName]]

        for (geneOrTF in c("Gene", "TF")) {
            resultsList[[quadrantResultsName]][[
                paste0("top", geneOrTF, "s")
            ]] <- .plotQuadrantLinkedDNAMethylationSitesHeatmaps(
                TENETMultiAssayExperiment,
                hyperHypo,
                geneOrTF,
                topGeneNumber,
                geneIDNameDF,
                quadrantSigLinkZScores,
                quadrantResultsName
            )
        }
    }

    ## Add the results list to the MultiAssayExperiment
    TENETMultiAssayExperiment@metadata$
        step7TopGenesOverlappingLinkedDNAMethylationSitesHeatmaps <- resultsList

    return(TENETMultiAssayExperiment)
}
