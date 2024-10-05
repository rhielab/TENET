#' Create bed-formatted interact files which can be loaded on the UCSC
#' Genome Browser to display links between top genes and transcription factors
#' and their linked RE DNA methylation sites
#'
#' This function takes the top genes/transcription factors (TFs) by number of
#' linked RE DNA methylation sites identified by the
#' `step6DNAMethylationSitesPerGeneTabulation` function, up to the number
#' specified by the user, and generates bed-formatted interact files (see
#' <https://genome.ucsc.edu/goldenPath/help/interact.html>) that can be
#' uploaded to the UCSC Genome Browser (<https://genome.ucsc.edu>) to visualize
#' the links between each of the top specified genes/TFs and the RE DNA
#' methylation sites linked to them for both of the hyper- or hypomethylated G+
#' analysis quadrants, as selected by the user.
#'
#' @param TENETMultiAssayExperiment Specify a MultiAssayExperiment object
#' containing expression and methylation SummarizedExperiment objects, such as
#' one created by the TCGADownloader function. This MultiAssayExperiment object
#' should also contain the results from the `step5OptimizeLinks` and
#' `step6DNAMethylationSitesPerGeneTabulation` functions in its metadata.
#' @param outputDirectory Specify the path to the output directory in which to
#' save the .bed files created by this function. The directory will be created
#' if it does not exist.
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
#'
#'
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
#' @param hypermethGplusAnalysis Set to TRUE to create interact files showing
#' links between the top genes/TFs by most RE hypermethylated RE DNA
#' methylation sites with G+ links, and these linked RE DNA methylation sites.
#' Defaults to TRUE.
#' @param hypomethGplusAnalysis Set to TRUE to create interact files showing
#' links between the top genes/TFs by most hypomethylated RE DNA methylation
#' sites with G+ links, and these linked RE DNA methylation sites. Defaults to
#' TRUE.
#' @param topGeneNumber Specify the number of top genes/TFs, based on the most
#' linked RE DNA methylation sites of a given analysis type, for which to
#' generate interact files showing the links between those genes and each of
#' their linked RE DNA methylation sites. Defaults to 10.
#' @return Outputs .bed formatted interact files to upload to the UCSC Genome
#' Browser to the specified output directory. These files display the
#' interactions between the top genes/TFs and their linked RE DNA methylation
#' sites for the given analysis types. Returns a list of lists named after each
#' selected analysis type, each containing the file paths to the created .bed
#' files for top genes and top TFs for that analysis type.
#' @export
#'
#' @examplesIf interactive()
#' ## This example uses the example MultiAssayExperiment provided in the
#' ## TENET.ExperimentHub package to create UCSC Genome Browser interact files
#' ## for the top 10 genes/TFs by number of linked hyper- or hypomethylated RE
#' ## DNA methylation sites. The interact files for the top genes/TFs will be
#' ## saved in the user's working directory. Gene names and locations, and the
#' ## locations of RE DNA methylation sites, will be retrieved from the
#' ## rowRanges of the 'expression' and 'methylation' SummarizedExperiment
#' ## objects in the example MultiAssayExperiment.
#'
#' ## Load the example TENET MultiAssayExperiment object
#' ## from the TENET.ExperimentHub package
#' exampleTENETMultiAssayExperiment <-
#'     TENET.ExperimentHub::exampleTENETMultiAssayExperiment()
#'
#' ## Use the example dataset to create and save the UCSC Genome Browser
#' ## interact files
#' filePaths <- step7TopGenesUCSCBedFiles(
#'     TENETMultiAssayExperiment = exampleTENETMultiAssayExperiment,
#'     outputDirectory = "."
#' )
#'
#' ## Get the path to the bed file for the top TFs by number of
#' ## hypomethylated G+ RE DNA methylation sites
#' filePaths$hypoGplus$topTFs
#'
#' ## This example uses the example MultiAssayExperiment provided in the
#' ## TENET.ExperimentHub package to create UCSC Genome Browser interact files
#' ## for the top 5 genes/TFs by number of linked hypomethylated RE DNA
#' ## methylation sites only. The interact files for the top genes/TFs will be
#' ## saved in the user's working directory. Gene names and locations will be
#' ## retrieved from the rowRanges of the 'expression' and 'methylation'
#' ## SummarizedExperiment objects in the example MultiAssayExperiment. RE DNA
#' ## methylation site IDs and locations will be retrieved from the HM450 array
#' ## via the sesameData package.
#'
#' ## Load the example TENET MultiAssayExperiment object
#' ## from the TENET.ExperimentHub package
#' exampleTENETMultiAssayExperiment <-
#'     TENET.ExperimentHub::exampleTENETMultiAssayExperiment()
#'
#' ## Use the example dataset to create and save the UCSC Genome Browser
#' ## interact files
#' filePaths <- step7TopGenesUCSCBedFiles(
#'     TENETMultiAssayExperiment = exampleTENETMultiAssayExperiment,
#'     outputDirectory = ".",
#'     DNAMethylationArray = "HM450",
#'     hypermethGplusAnalysis = FALSE,
#'     topGeneNumber = 5
#' )
#'
#' ## Get the path to the bed file for the top TFs by number of
#' ## hypomethylated G+ RE DNA methylation sites.
#' ## Note: Since we performed analyses only using TFs in the step 3 function,
#' ## the top genes are all TFs, so topTFs will be NA here, and topGenes
#' ## should be used instead.
#' filePaths$hypoGplus$topGenes
step7TopGenesUCSCBedFiles <- function(
    TENETMultiAssayExperiment,
    outputDirectory,
    geneAnnotationDataset = NA,
    DNAMethylationArray = NA,
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

    if (missing(outputDirectory)) {
        .stopNoCall("The outputDirectory parameter must be specified.")
    }

    ## Trim extra slashes from the output directory. file.path does it on
    ## Windows, but not on other OSes.
    outputDirectory <- trimws(
        file.path(outputDirectory),
        which = "right", whitespace = "/"
    )

    ## Create the output directory if it doesn't exist
    R.utils::mkdirs(outputDirectory)

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

    ## Create an empty list to hold the results from this step 7 function
    resultsList <- list()

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

        for (geneOrTF in c("Gene", "TF")) {
            ## Get the IDs of the top genes/TFs in this quadrant
            topQuadrantGeneOrTFIDs <- .getQuadrantTopGenesOrTFs(
                TENETMultiAssayExperiment, geneOrTF, hyperHypo,
                topGeneNumber
            )$geneID
            if (.isSingleNA(topQuadrantGeneOrTFIDs)) {
                resultsList[[
                    paste0(hyperHypo, "Gplus")
                ]][[
                    paste0("top", geneOrTF, "s")
                ]] <- NA
                next
            }

            ## Create data frames listing the top genes/TFs in one column
            ## (repeated), and all the RE DNA methylation sites linked to them
            ## in a second column

            ## Create empty data frames
            topQuadrantGenesOrTFsIntersectDF <- data.frame(
                "geneID" = character(),
                "geneName" = character(),
                "DNAMethylationSiteID" = character(),
                stringsAsFactors = FALSE
            )

            ## Add the RE DNA methylation sites linked to each of the top
            ## genes/TFs to the data frame
            for (geneID in topQuadrantGeneOrTFIDs) {
                ## Get the gene name corresponding to the gene ID
                geneName <- geneIDNameDF[geneID, "geneName"]

                ## Get a list of all the RE DNA methylation sites linked to the
                ## given gene
                linkedDNAMethylationSites <- unique(
                    quadrantSigLinkZScores[
                        quadrantSigLinkZScores$geneID == geneID,
                        "DNAMethylationSiteID"
                    ]
                )

                ## Create a temporary data frame with the 3 vectors of info
                tempDF <- data.frame(
                    "geneID" = geneID,
                    "geneName" = geneName,
                    "DNAMethylationSiteID" = linkedDNAMethylationSites,
                    stringsAsFactors = FALSE
                )

                ## Bind the tempDF to the IntersectDF
                topQuadrantGenesOrTFsIntersectDF <- rbind(
                    topQuadrantGenesOrTFsIntersectDF,
                    tempDF
                )
            }

            ## For each of the genes in the intersect DFs, get the chromosome,
            ## "start", and "end" of the gene. Start and end will both be the
            ## TSS.
            topQuadrantGenesOrTFsIntersectDF$geneChr <- geneIDNameDF[
                topQuadrantGenesOrTFsIntersectDF$geneID, "chromosome"
            ]

            topQuadrantGenesOrTFsIntersectDF$geneStart <- geneIDNameDF[
                topQuadrantGenesOrTFsIntersectDF$geneID, "TSS"
            ]

            topQuadrantGenesOrTFsIntersectDF$geneEnd <- geneIDNameDF[
                topQuadrantGenesOrTFsIntersectDF$geneID, "TSS"
            ]

            ## For each of the RE DNA methylation sites in the intersect DFs,
            ## get the chromosome, "start", and "end" of the RE DNA methylation
            ## site
            topQuadrantGenesOrTFsIntersectDF$siteChr <- methSiteIDdf[
                topQuadrantGenesOrTFsIntersectDF$DNAMethylationSiteID,
                "chromosome"
            ]

            topQuadrantGenesOrTFsIntersectDF$siteStart <- methSiteIDdf[
                topQuadrantGenesOrTFsIntersectDF$DNAMethylationSiteID, "start"
            ]

            topQuadrantGenesOrTFsIntersectDF$siteEnd <- methSiteIDdf[
                topQuadrantGenesOrTFsIntersectDF$DNAMethylationSiteID, "end"
            ]

            ## Use the grDevices::rainbow color function to set up a gradient of
            ## colors equal to the number of genes analyzed and assign each
            ## gene a color
            rainbowNumericColorGrad <- grDevices::rainbow(topGeneNumber)
            names(rainbowNumericColorGrad) <- topQuadrantGeneOrTFIDs

            ## Assemble the output file for the top genes/TFs.
            ## Note: this code cannot be changed to use rtracklayer because
            ## UCSC interact files are not standard bed or bedGraph files, and
            ## rtracklayer does not support adding custom columns when exporting
            ## a bed or bedGraph file.
            topQuadrantGenesOrTFsOutputDF <- data.frame(
                "chrom" = topQuadrantGenesOrTFsIntersectDF$geneChr,
                "chromStart" = ifelse(
                    topQuadrantGenesOrTFsIntersectDF$geneChr ==
                        topQuadrantGenesOrTFsIntersectDF$siteChr,
                    ifelse(
                        (topQuadrantGenesOrTFsIntersectDF$siteStart <
                            topQuadrantGenesOrTFsIntersectDF$geneStart
                        ),
                        topQuadrantGenesOrTFsIntersectDF$siteStart,
                        topQuadrantGenesOrTFsIntersectDF$geneStart
                    ),
                    topQuadrantGenesOrTFsIntersectDF$geneStart - 1
                ),
                "chromEnd" = ifelse(
                    topQuadrantGenesOrTFsIntersectDF$geneChr ==
                        topQuadrantGenesOrTFsIntersectDF$siteChr,
                    ifelse(
                        (topQuadrantGenesOrTFsIntersectDF$siteStart <
                            topQuadrantGenesOrTFsIntersectDF$geneStart
                        ),
                        topQuadrantGenesOrTFsIntersectDF$geneStart - 1,
                        topQuadrantGenesOrTFsIntersectDF$siteStart
                    ),
                    topQuadrantGenesOrTFsIntersectDF$geneEnd
                ),
                "name" = paste(
                    topQuadrantGenesOrTFsIntersectDF$geneID,
                    topQuadrantGenesOrTFsIntersectDF$DNAMethylationSiteID,
                    "link",
                    sep = "_"
                ),
                "score" = 0,
                "value" = 0,
                "exp" = ".",
                "color" = rainbowNumericColorGrad[
                    topQuadrantGenesOrTFsIntersectDF$geneID
                ],
                "sourceChrom" = topQuadrantGenesOrTFsIntersectDF$geneChr,
                "sourceStart" = (
                    topQuadrantGenesOrTFsIntersectDF$geneStart - 1
                ),
                "sourceEnd" = topQuadrantGenesOrTFsIntersectDF$geneEnd,
                "sourceName" = topQuadrantGenesOrTFsIntersectDF$geneName,
                "sourceStrand" = ".",
                "targetChrom" = topQuadrantGenesOrTFsIntersectDF$siteChr,
                "targetStart" = topQuadrantGenesOrTFsIntersectDF$
                    siteStart,
                "targetEnd" = topQuadrantGenesOrTFsIntersectDF$siteEnd,
                "targetName" = topQuadrantGenesOrTFsIntersectDF$
                    DNAMethylationSiteID,
                "targetStrand" = ".",
                stringsAsFactors = FALSE
            )

            ## Create text for the header line
            topQuadrantGenesOrTFsHeaderText <- paste0(
                "track type=interact name=\"TENET",
                tools::toTitleCase(hyperHypo),
                "G+Interactions\" description=\"TENET top ",
                ifelse(geneOrTF == "Gene", "gene", "TF"),
                " to RE DNA methylation site links\""
            )

            ## Create a file name for the output bed file. UCSC documentation
            ## uses an .inter.bed extension.
            topQuadrantGenesOrTFsBedFileName <- file.path(
                outputDirectory,
                paste0(
                    "top",
                    tools::toTitleCase(hyperHypo),
                    "Gplus",
                    geneOrTF,
                    "ToDNAMethylationSiteLinks.inter.bed"
                )
            )

            ## Add the header line to the new bed file
            write(
                topQuadrantGenesOrTFsHeaderText,
                file = topQuadrantGenesOrTFsBedFileName
            )

            ## Write the info to the file
            utils::write.table(
                topQuadrantGenesOrTFsOutputDF,
                file = topQuadrantGenesOrTFsBedFileName,
                append = TRUE,
                row.names = FALSE,
                col.names = FALSE,
                quote = FALSE
            )

            ## Add the file name to the results list
            resultsList[[paste0(hyperHypo, "Gplus")]][[
                paste0("top", geneOrTF, "s")
            ]] <- topQuadrantGenesOrTFsBedFileName
        }
    }

    return(resultsList)
}
