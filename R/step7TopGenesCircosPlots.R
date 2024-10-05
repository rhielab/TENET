## Internal functions used by step7TopGenesCircosPlots

## Internal function to generate and save a Circos plot for a given gene
## to each of its linked RE DNA methylation sites
.geneCircosFunction <- function(
    geneID,
    geneIDNameDF,
    quadrantSigLinkZScores,
    methSiteIDdf) {
    ## Get the gene name corresponding to the gene ID
    geneName <- geneIDNameDF[
        geneIDNameDF$geneID == geneID, "geneName"
    ]

    ## Get a list of all the RE DNA methylation sites linked to the given gene
    linkedDNAMethylationSites <- unique(
        quadrantSigLinkZScores[
            quadrantSigLinkZScores$geneID == geneID,
            "DNAMethylationSiteID"
        ]
    )

    ## Create lists of the gene chromosome and start position repeated a number
    ## of times equal to the number of linked RE DNA methylation sites
    circosStartChromosome <- rep(
        as.character(geneIDNameDF[
            geneIDNameDF$geneID == geneID, "chromosome"
        ]),
        length(linkedDNAMethylationSites)
    )

    circosStartPosition <- rep(
        geneIDNameDF[geneIDNameDF$geneID == geneID, "TSS"],
        length(linkedDNAMethylationSites)
    )

    ## List the chromosome and start for RE DNA methylation sites linked to the
    ## given gene as end points on the Circos plot
    circosEndChromosome <- methSiteIDdf[linkedDNAMethylationSites, "chromosome"]
    circosEndPosition <- methSiteIDdf[linkedDNAMethylationSites, "start"]

    ## Create the RE DNA methylation sites and genes table for RCircos
    RCircosMethSites <- data.frame(
        chr = circosEndChromosome,
        start = (circosEndPosition - 1),
        end = (circosEndPosition + 1),
        name = linkedDNAMethylationSites
    )

    RCircosGenes <- data.frame(
        chr = circosStartChromosome,
        start = (circosStartPosition - 1),
        end = (circosStartPosition + 1),
        name = paste0(geneName, " (", geneID, ")")
    )

    ## Combine the RE DNA methylation site and gene data frames
    methSitesAndGenes <- rbind(RCircosMethSites, RCircosGenes)

    ## Get unique entries for each gene
    methSitesAndGenes <- unique(methSitesAndGenes)

    ## Create a new links data frame with info from the previously
    ## separated data frames
    links <- RCircosGenes
    links$V1 <- RCircosMethSites$chr
    links$V2 <- RCircosMethSites$start
    links$V3 <- RCircosMethSites$end
    links$V4 <- RCircosMethSites$name
    colnames(links) <- c(
        "geneChromosome",
        "geneChromStart",
        "geneChromEnd",
        "geneNameAndID",
        "methSiteChromosome",
        "methSiteChromStart",
        "methSiteChromEnd",
        "methSiteID"
    )

    ## Set the parameters
    params <- RCircos::RCircos.Get.Plot.Parameters()
    params$base.per.unit <- 4500
    RCircos::RCircos.Reset.Plot.Parameters(params)

    ## Use the RCircos package to create the Circos plot
    .newInvisibleRecordablePlot(width = 8, height = 8)
    RCircos::RCircos.Set.Plot.Area()
    RCircos::RCircos.Chromosome.Ideogram.Plot()
    RCircos::RCircos.Gene.Connector.Plot(
        methSitesAndGenes,
        track.num = 1,
        side = "in"
    )

    linkData <- links[, c(1, 2, 3, 5, 6, 7)]
    linkData$PlotColor <- "red"
    RCircos::RCircos.Link.Plot(
        linkData,
        track.num = 2,
        by.chromosome = FALSE
    )

    circosPlot <- .recordTENETSavedSizePlot()
    grDevices::dev.off()

    return(circosPlot)
}

## Main step7TopGenesCircosPlots function

#' Generate Circos plots displaying the links between top identified genes and
#' each of the RE DNA methylation sites linked to them
#'
#' This function takes the top genes/TFs by number of linked regulatory element
#' DNA methylation sites identified by the
#' `step6DNAMethylationSitesPerGeneTabulation` function, up to the number
#' specified by the user, and generates Circos plots for each gene showing the
#' genomic links between each gene and each RE DNA methylation site linked to
#' the gene for the analysis types specified.
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
#' @param DNAMethylationArray Specify the name of a DNA methylation probe array
#' supported by the sesameData package (see
#' `?sesameData::sesameData_getManifestGRanges`). If an array is specified, RE
#' DNA methylation sites and their locations in that array's manifest are
#' cross-referenced with RE DNA methylation site IDs included in the rownames
#' of the methylation dataset provided in the "methylation"
#' SummarizedExperiment object within the TENETMultiAssayExperiment object, and
#' only those overlapping will be considered for analysis. If set to NA, all RE
#' DNA methylation sites with locations listed in the rowRanges of the
#' "methylation" SummarizedExperiment object are used. Defaults to NA.
#' @param hypermethGplusAnalysis Set to TRUE to create Circos plots displaying
#' links between the top genes/TFs by most hypermethylated RE DNA methylation
#' sites with G+ links and their linked RE DNA methylation sites of that type.
#' Defaults to TRUE.
#' @param hypomethGplusAnalysis Set to TRUE to create Circos plots displaying
#' links between the top genes/TFs by most hypomethylated RE DNA methylation
#' sites with G+ links and their linked RE DNA methylation sites of that type.
#' Defaults to TRUE.
#' @param topGeneNumber Specify the number of top genes/TFs, based on the most
#' linked RE DNA methylation sites of a given analysis type, for which to
#' generate Circos plots showing the links between the genes and each of their
#' linked RE DNA methylation sites. Defaults to 10.
#' @param coreCount Argument passed as the mc.cores argument to mclapply. See
#' `?parallel::mclapply` for more details. Defaults to 1.
#' @return Returns the MultiAssayExperiment object given as the
#' TENETMultiAssayExperiment argument with an additional list of information
#' named 'step7TopGenesCircosPlots' in its metadata with the output of this
#' function. This list is subdivided into hypermethGplus or hypomethGplus
#' results as selected by the user, which are further subdivided into lists
#' with data for the top overall genes and for top TF genes only. Each of these
#' lists contain Circos plots for each of the top genes/TFs visualizing the
#' links between the genes and their linked RE DNA methylation sites for the
#' selected analysis types.
#' @export
#'
#' @examplesIf interactive()
#' ## This example uses the example MultiAssayExperiment provided in the
#' ## TENET.ExperimentHub package to create Circos plots for the top 10
#' ## genes/TFs, by number of linked hyper- or hypomethylated RE DNA
#' ## methylation sites. Gene names and locations and RE DNA methylation site
#' ## locations will be retrieved from the rowRanges of the 'expression' and
#' ## 'methylation' SummarizedExperiment objects in the example
#' ## MultiAssayExperiment. The analysis will be performed using one CPU core.
#'
#' ## Load the example TENET MultiAssayExperiment object
#' ## from the TENET.ExperimentHub package
#' exampleTENETMultiAssayExperiment <-
#'     TENET.ExperimentHub::exampleTENETMultiAssayExperiment()
#'
#' ## Use the example dataset to create Circos plots
#' returnValue <- step7TopGenesCircosPlots(
#'     TENETMultiAssayExperiment = exampleTENETMultiAssayExperiment
#' )
#'
#' ## This example uses the example MultiAssayExperiment provided in the
#' ## TENET.ExperimentHub package to create Circos plots for the top 5
#' ## genes/TFs, by number of linked hypomethylated RE DNA methylation sites
#' ## only. Gene names and locations will be retrieved from the rowRanges of
#' ## the 'expression' and 'methylation' SummarizedExperiment objects in the
#' ## example MultiAssayExperiment. RE DNA methylation site IDs and locations
#' ## will be retrieved from the HM450 array via the sesameData package. Eight
#' ## CPU cores will be used to create plots.
#'
#' ## Load the example TENET MultiAssayExperiment object
#' ## from the TENET.ExperimentHub package
#' exampleTENETMultiAssayExperiment <-
#'     TENET.ExperimentHub::exampleTENETMultiAssayExperiment()
#'
#' ## Use the example dataset to create Circos plots
#' returnValue <- step7TopGenesCircosPlots(
#'     TENETMultiAssayExperiment = exampleTENETMultiAssayExperiment,
#'     DNAMethylationArray = "HM450",
#'     hypermethGplusAnalysis = FALSE,
#'     topGeneNumber = 5,
#'     coreCount = 8
#' )
step7TopGenesCircosPlots <- function(
    TENETMultiAssayExperiment,
    DNAMethylationArray = NA,
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

    ## Get methylation site IDs and names from the MAE, or methylation array if
    ## provided
    methSiteIDdf <- .getMethSiteIDsAndLocations(
        TENETMultiAssayExperiment, DNAMethylationArray
    )

    ## Create an environment to store data from the RCircos package
    circosDataEnv <- new.env(parent = emptyenv())

    ## Load the hg38 ideogram dataset
    utils::data(
        "UCSC.HG38.Human.CytoBandIdeogram",
        package = "RCircos",
        envir = circosDataEnv
    )

    ## Load the RCircos environment. RCircos explicitly reads it from globalenv,
    ## see https://github.com/cran/RCircos/blob/c46b2d54eb6123f5a365b4cc6139b4225c23de08/R/RCircosMain.R#L201
    oldRCircos.Env <- NA
    if (exists("RCircos.Env", envir = globalenv())) {
        oldRCircos.Env <- get("RCircos.Env", envir = globalenv())
    }
    assign("RCircos.Env", RCircos::RCircos.Env, envir = globalenv())

    ## Configure the RCircos settings
    RCircos::RCircos.Set.Core.Components(
        cyto.info = circosDataEnv$UCSC.HG38.Human.CytoBandIdeogram,
        chr.exclude = NULL,
        tracks.inside = 2,
        tracks.outside = 0
    )

    ## Create an empty list to hold the results from this step 7 function
    resultsList <- list()

    ## Generate Circos plots for the selected analysis types
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
            ## Get the IDs of the top genes/TFs in this quadrant
            topQuadrantGeneOrTFIDs <- .getQuadrantTopGenesOrTFs(
                TENETMultiAssayExperiment, geneOrTF, hyperHypo,
                topGeneNumber
            )$geneID
            if (.isSingleNA(topQuadrantGeneOrTFIDs)) {
                resultsList[[quadrantResultsName]][[
                    paste0("top", geneOrTF, "s")
                ]] <- NA
                next
            }

            ## Generate the plots for the genes of interest and add them to the
            ## results list
            resultsList[[quadrantResultsName]][[
                paste0("top", geneOrTF, "s")
            ]] <-
                parallel::mclapply(
                    topQuadrantGeneOrTFIDs,
                    .geneCircosFunction,
                    geneIDNameDF = geneIDNameDF,
                    quadrantSigLinkZScores = quadrantSigLinkZScores,
                    methSiteIDdf = methSiteIDdf,
                    mc.cores = coreCount
                )

            names(
                resultsList[[quadrantResultsName]][[
                    paste0("top", geneOrTF, "s")
                ]]
            ) <- topQuadrantGeneOrTFIDs
        }
    }

    ## Restore the old global RCircos environment, if there was one
    if (is.na(oldRCircos.Env)) {
        rm("RCircos.Env", envir = globalenv())
    } else {
        assign("RCircos.Env", oldRCircos.Env, envir = globalenv())
    }

    ## Add the results list to the MultiAssayExperiment
    TENETMultiAssayExperiment@metadata$
        step7TopGenesCircosPlots <- resultsList

    return(TENETMultiAssayExperiment)
}
