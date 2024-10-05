## Internal functions used by step7LinkedDNAMethylationSitesMotifSearching

## Internal function to get the DNA string of the motif matches
.internalMotifGrabber <- function(istart, iend, DNAString) {
    return(as.character(DNAString[istart:iend]))
}

## Internal function to find the motif occurrences in the specified vicinity of
## each linked RE DNA methylation site. Note: All instances of + 1 and -1 are
## necessary and correct.
.findMotifSurroundingMethSite <- function(
    DNAMethylationSiteID,
    dfOrCountOutput,
    genome,
    motifPWM,
    quadrantMethSiteDatasetLinkedToGene,
    matchPWMMinScore,
    distanceFromREDNAMethylationSites) {
    ## Get the methylation site ID, chromosome, start, and end locations
    chr <- as.character(quadrantMethSiteDatasetLinkedToGene[
        DNAMethylationSiteID, 2
    ])
    start <- as.numeric(quadrantMethSiteDatasetLinkedToGene[
        DNAMethylationSiteID, 3
    ])
    end <- as.numeric(quadrantMethSiteDatasetLinkedToGene[
        DNAMethylationSiteID, 4
    ])

    ## Subtract and add the buffer specified by the user with the
    ## distanceFromREDNAMethylationSites argument from the start and end values
    startBuffer <- start - distanceFromREDNAMethylationSites + 1
    endBuffer <- end + distanceFromREDNAMethylationSites

    ## Get the DNA string sequence for that segment from the human genome
    DNAString <- genome[[chr]][startBuffer:endBuffer]

    ## Do the motif search
    search <- Biostrings::matchPWM(motifPWM, DNAString, matchPWMMinScore)

    if (dfOrCountOutput == "df") {
        ## Return a data frame with the information about the found motifs

        ## Get the start and end locations where the motif was found
        starts <- IRanges::start(search@ranges)
        ends <- IRanges::end(search@ranges)

        ## Get the motif matches
        motifMatches <- mapply(
            .internalMotifGrabber,
            starts,
            ends,
            MoreArgs = list("DNAString" = DNAString)
        )

        ## Create a data frame with results and return it. The start and end
        ## position calculations look strange but are correct, since matchPWM
        ## is being run only on the DNA string of the region surrounding the RE
        ## DNA methylation site, so we need to add the start of the region
        ## surrounding the RE DNA methylation site to the position of the motif
        ## in the string to get the location of the motif on the whole
        ## chromosome.
        if (length(as.character(search)) == 0) {
            ## No results are found, so try creating an empty data frame
            ## with the columns of interest, but 0 rows:
            emptyDF <- data.frame(matrix(ncol = 5, nrow = 0))
            colnames(emptyDF) <- c(
                "DNAMethylationSiteID",
                "motifChromosome",
                "motifStartPosition",
                "motifEndPosition",
                "motifDNASequence"
            )
            return(emptyDF)
        } else {
            return(data.frame(
                "DNAMethylationSiteID" = DNAMethylationSiteID,
                "motifChromosome" = chr,
                "motifStartPosition" = (starts + startBuffer - 1),
                "motifEndPosition" = (ends + startBuffer - 1),
                "motifDNASequence" = motifMatches,
                stringsAsFactors = FALSE
            ))
        }
    } else if (dfOrCountOutput == "count") {
        ## Return the count of motifs found in the vicinity of the RE DNA
        ## methylation site
        if (length(as.character(search)) == 0) {
            return(0)
        } else {
            return(length(search))
        }
    }
}

## Main step7LinkedDNAMethylationSitesMotifSearching function

#' Perform motif searching for transcription factor motifs in the vicinity of
#' RE DNA methylation sites linked to the specified transcription factors
#'
#' This function takes a named list of transcription factors (TFs) and their
#' binding motifs, and for each of the TFs, it identifies if the specified
#' motif for that TF is found within a user-specified distance to each RE DNA
#' methylation site linked to that TF from both of the hyper- or hypomethylated
#' G+ analysis quadrants, as selected by the user.
#'
#' @param TENETMultiAssayExperiment Specify a MultiAssayExperiment object
#' containing expression and methylation SummarizedExperiment objects, such as
#' one created by the TCGADownloader function. This MultiAssayExperiment object
#' should also contain the results from the `step5OptimizeLinks` function
#' in its metadata.
#' @param TFMotifList Specify a named list mapping transcription factor gene
#' names and/or IDs to their respective motif position weight matrix (PWM). The
#' PWMs should be in the form of a 4xN matrix. If a PWM is not available, the
#' user can select one by querying the MotifDb package database. See the
#' example below for more infornation. **Note:** The name of each PWM in the
#' list must match the gene name or Ensembl ID of a gene in the
#' TENETMultiAssayExperiment with RE DNA methylation sites linked to it for the
#' specified analysis types.
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
#' @param hypermethGplusAnalysis Set to TRUE to do motif searching in the
#' vicinity of hypermethylated RE DNA methylation sites with G+ links to the TF
#' of interest. Defaults to TRUE.
#' @param hypomethGplusAnalysis Set to TRUE to do motif searching in the
#' vicinity of hypomethylated RE DNA methylation sites with G+ links to the TF
#' of interest. Defaults to TRUE.
#' @param distanceFromREDNAMethylationSites Specify a positive integer in base
#' pairs to be the distance from the linked RE DNA methylation sites within
#' which motif searching will be performed. Defaults to 100.
#' @param matchPWMMinScore Specify the `min.score` argument passed to the
#' matchPWM function for motif searching. See `?Biostrings::matchPWM` for more
#' details. Defaults to "75%".
#' @param coreCount Argument passed as the mc.cores argument to mclapply. See
#' `?parallel::mclapply` for more details. Defaults to 1.
#' @return Returns the MultiAssayExperiment object given as the
#' TENETMultiAssayExperiment argument with an additional list of information
#' named 'step7LinkedDNAMethylationSitesMotifSearching' in its metadata with
#' the output of this function. This list is subdivided into hypermethGplus or
#' hypomethGplus results as selected by the user. Each contains a list for each
#' specified TF gene with a seqLogo of the gene's motif PWM, a table of the
#' PWM, and two additional tables with information on motif occurrences found
#' in the vicinity of linked RE DNA methylation sites. The first table contains
#' information about each motif found in the vicinity of the linked RE DNA
#' methylation sites for the TF gene at the specified matchPWMMinScore. The
#' second lists each RE DNA methylation site linked to the TF gene for the
#' analysis types specified, and the number of motif occurrences found in the
#' vicinity of each of those RE DNA methylation sites.
#' @export
#'
#' @examplesIf interactive()
#' ## These examples require the MotifDb package
#' if (requireNamespace("MotifDb", quietly = TRUE)) {
#'     ## Show available motifs for example TF FOXA1
#'     names(MotifDb::query(MotifDb::MotifDb, "FOXA1"))
#'
#'     ## The seqLogos for all input motifs will also be included in the output
#'     ## of this function. Alternatively, individual motifs can be visualized
#'     ## with the seqLogo function from the seqLogo package.
#'     seqLogo::seqLogo(MotifDb::query(MotifDb::MotifDb, "FOXA1")[[3]])
#'
#'     ## Once PWMs have been selected for use, a list containing them must be
#'     ## created
#'     exampleTFMotifList <- list(
#'         "FOXA1" = MotifDb::query(MotifDb::MotifDb, "FOXA1")[[3]],
#'         "ESR1" = MotifDb::query(MotifDb::MotifDb, "ESR1")[[4]]
#'     )
#'
#'     ## This example uses the example MultiAssayExperiment provided in the
#'     ## TENET.ExperimentHub package to perform motif overlapping for all
#'     ## hyper- and hypomethylated RE DNA methylation sites linked to the
#'     ## FOXA1 and ESR1 genes using the motifs for each gene specified in the
#'     ## exampleTFMotifList. Gene names and locations, and the locations of RE
#'     ## DNA methylation sites, will be retrieved from the rowRanges of the
#'     ## 'expression' and 'methylation' SummarizedExperiment objects in the
#'     ## example MultiAssayExperiment. Regions within 100 bp of linked RE DNA
#'     ## methylation sites will be checked for motifs, and a similarity
#'     ## threshold of 75% will be used to identify motifs. The analysis will
#'     ## be performed using one CPU core.
#'
#'     ## Load the example TENET MultiAssayExperiment object
#'     ## from the TENET.ExperimentHub package
#'     exampleTENETMultiAssayExperiment <-
#'         TENET.ExperimentHub::exampleTENETMultiAssayExperiment()
#'
#'     ## Use the example dataset to do the motif searching
#'     returnValue <- step7LinkedDNAMethylationSitesMotifSearching(
#'         TENETMultiAssayExperiment = exampleTENETMultiAssayExperiment,
#'         TFMotifList = exampleTFMotifList
#'     )
#'
#'     ## This example uses the example MultiAssayExperiment provided in the
#'     ## TENET.ExperimentHub package to perform motif overlapping for only
#'     ## hypomethylated RE DNA methylation sites linked to the FOXA1 and ESR1
#'     ## genes using the motifs for each gene specified in the
#'     ## exampleTFMotifList. Gene names and locations, and the locations of RE
#'     ## DNA methylation sites, will be retrieved from the rowRanges of the
#'     ## 'expression' and 'methylation' SummarizedExperiment objects in the
#'     ## example MultiAssayExperiment. Regions within 50 bp of linked RE DNA
#'     ## methylation sites will be checked for motifs, and a similarity
#'     ## threshold of 80% will be used to identify motifs. The analysis will
#'     ## be performed using 8 CPU cores.
#'
#'     ## Load the example TENET MultiAssayExperiment object
#'     ## from the TENET.ExperimentHub package
#'     exampleTENETMultiAssayExperiment <-
#'         TENET.ExperimentHub::exampleTENETMultiAssayExperiment()
#'
#'     ## Use the example dataset to do the motif searching
#'     returnValue <- step7LinkedDNAMethylationSitesMotifSearching(
#'         TENETMultiAssayExperiment = exampleTENETMultiAssayExperiment,
#'         TFMotifList = exampleTFMotifList,
#'         hypermethGplusAnalysis = FALSE,
#'         distanceFromREDNAMethylationSites = 50,
#'         matchPWMMinScore = "80%",
#'         coreCount = 8
#'     )
#' } else {
#'     message(
#'         "These examples require the MotifDb package. Run ",
#'         'BiocManager::install("MotifDb") and rerun these examples.'
#'     )
#' }
step7LinkedDNAMethylationSitesMotifSearching <- function(
    TENETMultiAssayExperiment,
    TFMotifList,
    geneAnnotationDataset = NA,
    DNAMethylationArray = NA,
    hypermethGplusAnalysis = TRUE,
    hypomethGplusAnalysis = TRUE,
    distanceFromREDNAMethylationSites = 100,
    matchPWMMinScore = "75%",
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

    if (missing(TFMotifList)) {
        .stopNoCall("The TFMotifList parameter must be specified.")
    }

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

    ## Get the human genome
    genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38

    ## Create an empty list to hold the results from this step 7 function
    resultsList <- list()

    for (i in seq_along(TFMotifList)) {
        motifPWM <- TFMotifList[[i]]
        TFGene <- names(TFMotifList[i]) ## One bracket to keep names
        if (is.null(TFGene)) {
            ## There are no names
            .stopNoCall(
                "The TFMotifList argument must be a named list mapping TF ",
                "gene names and/or IDs to their respective motif PWM."
            )
        }

        ## Select a PWM to use in the analysis. Check if the user has
        ## happened to supply their own PWM.
        if (!is.matrix(motifPWM) || nrow(motifPWM) != 4) {
            ## The input PWM is not in the format of a matrix. Cancel the
            ## function and alert the user.
            .stopNoCall(
                "Input PWM for gene ", TFGene, " is not a 4xN matrix. ",
                "Please specify a 4xN matrix for this argument."
            )
        }

        ## Create a seqLogo of the motif PWM and add it to the results list
        .newInvisibleRecordablePlot(height = 3.5, width = 6)
        seqLogo::seqLogo(motifPWM)
        thisSeqLogo <- .recordTENETSavedSizePlot()

        ## Close the plot
        grDevices::dev.off()

        ## If the user has input a gene name, get the gene ID for it
        if (TFGene %in% geneIDNameDF$geneName) {
            geneID <- geneIDNameDF[
                which(geneIDNameDF$geneName == TFGene),
                "geneID"
            ]
        } else {
            geneID <- TFGene

            TFGene <- geneIDNameDF[
                which(geneIDNameDF$geneID == TFGene),
                "geneName"
            ]
        }

        ## For each of the analysis types that is selected, assemble the data
        ## frames of the complete motif occurrences and counts

        ## Check the step 5 results and find all the RE DNA methylation sites
        ## linked to the gene in each of the analysis quadrants
        for (hyperHypo in analysisTypes) {
            quadrantResultsName <- paste0(hyperHypo, "methGplusResults")

            ## Save the seqLogo for the TF's motif
            resultsList[[quadrantResultsName]][[TFGene]]$motifSeqLogo <-
                thisSeqLogo

            ## Ensure the quadrant's results are present in step 5
            .ensureStepPresent(
                TENETMultiAssayExperiment,
                stepName = "step5OptimizeLinks",
                substepName = quadrantResultsName
            )

            ## Load the quadrant's significant links from step 5
            quadrantSigLinkZScores <- TENETMultiAssayExperiment@metadata$
                step5OptimizeLinks[[quadrantResultsName]]

            ## Identify the RE DNA methylation sites linked to the specific gene
            quadrantMethSitesLinkedToGene <- quadrantSigLinkZScores[
                quadrantSigLinkZScores$geneID == geneID,
                "DNAMethylationSiteID"
            ]

            ## If there are sites which aren't found in the specified
            ## methSiteIDdf, return a warning with those RE DNA methylation
            ## sites
            quadrantMethSitesLinkedToGeneNotFound <-
                quadrantMethSitesLinkedToGene[
                    !(quadrantMethSitesLinkedToGene %in%
                        methSiteIDdf$DNAMethylationSiteID
                    )
                ]

            if (length(quadrantMethSitesLinkedToGeneNotFound) > 0) {
                .warningNoCall(
                    "Genomic locations for the following RE DNA methylation ",
                    "sites linked to ", TFGene, " were not identified: ",
                    paste(
                        quadrantMethSitesLinkedToGeneNotFound,
                        collapse = ", "
                    ),
                    ". Please check the specified DNAMethylationArray ",
                    "argument. These methylation sites have been excluded ",
                    "from this analysis."
                )
            }

            ## Get the RE DNA methylation sites which are found in the
            ## specified methSiteIDdf
            quadrantMethSitesLinkedToGene <- quadrantMethSitesLinkedToGene[
                (quadrantMethSitesLinkedToGene %in%
                    methSiteIDdf$DNAMethylationSiteID
                )
            ]

            ## Create a data frame with the RE DNA methylation sites as the
            ## first column
            quadrantMethSiteDatasetLinkedToGene <- data.frame(
                "DNAMethylationSiteID" = sort(quadrantMethSitesLinkedToGene),
                stringsAsFactors = FALSE
            )

            ## Add the RE DNA methylation site location info
            quadrantMethSiteDatasetLinkedToGene$seqnames <- methSiteIDdf[
                quadrantMethSiteDatasetLinkedToGene$DNAMethylationSiteID,
                "chromosome"
            ]

            quadrantMethSiteDatasetLinkedToGene$start <- methSiteIDdf[
                quadrantMethSiteDatasetLinkedToGene$DNAMethylationSiteID,
                "start"
            ]

            quadrantMethSiteDatasetLinkedToGene$end <- methSiteIDdf[
                quadrantMethSiteDatasetLinkedToGene$DNAMethylationSiteID,
                "end"
            ]

            ## Set the rownames to be the methylation site IDs
            rownames(quadrantMethSiteDatasetLinkedToGene) <-
                quadrantMethSiteDatasetLinkedToGene$DNAMethylationSiteID

            ## Get the list of the RE DNA methylation site info DFs
            listOfMethSiteInfo <- parallel::mclapply(
                rownames(quadrantMethSiteDatasetLinkedToGene),
                .findMotifSurroundingMethSite,
                "df",
                genome,
                motifPWM,
                quadrantMethSiteDatasetLinkedToGene,
                matchPWMMinScore,
                distanceFromREDNAMethylationSites,
                mc.cores = coreCount
            )

            ## Rebind them into a data frame
            quadrantMethSiteMotifsDF <- do.call(rbind, listOfMethSiteInfo)

            if (is.null(quadrantMethSiteMotifsDF)) {
                ## If there are no motifs found, recreate an empty data frame
                ## that column
                quadrantMethSiteMotifsDF <- data.frame(
                    "DNAMethylationSiteID" = NULL,
                    "seqnames" = NULL,
                    "start" = NULL,
                    "end" = NULL,
                    stringsAsFactors = FALSE
                )
            }

            ## Get a list of the counts of the number of motifs found per DNA
            ## methylation site
            countOfMotifsPerMethSiteList <- parallel::mclapply(
                rownames(quadrantMethSiteDatasetLinkedToGene),
                .findMotifSurroundingMethSite,
                "count",
                genome,
                motifPWM,
                quadrantMethSiteDatasetLinkedToGene,
                matchPWMMinScore,
                distanceFromREDNAMethylationSites,
                mc.cores = coreCount
            )

            ## Convert the list into a vector
            countOfMotifsPerMethSite <- unlist(c(
                countOfMotifsPerMethSiteList
            ))

            ## Convert the list into a data frame with the RE DNA methylation
            ## sites as the rownames
            quadrantMotifCountPerMethSiteDF <- data.frame(
                "DNAMethylationSiteID" = rownames(
                    quadrantMethSiteDatasetLinkedToGene
                ),
                "motifCount" = countOfMotifsPerMethSite,
                stringsAsFactors = FALSE
            )

            ## Add the quadrant's results to the list
            resultsList[[quadrantResultsName]][[TFGene]]$
                DNAMethylationSiteMotifOccurrences <- quadrantMethSiteMotifsDF
            resultsList[[quadrantResultsName]][[TFGene]]$
                totalMotifOccurrencesPerREDNAMethylationSite <-
                quadrantMotifCountPerMethSiteDF
        }
    }

    ## Add the results list to the MultiAssayExperiment
    TENETMultiAssayExperiment@metadata$
        step7LinkedDNAMethylationSitesMotifSearching <- resultsList

    return(TENETMultiAssayExperiment)
}
